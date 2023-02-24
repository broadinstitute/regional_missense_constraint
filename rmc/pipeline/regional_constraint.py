import argparse
import logging

import hail as hl

from gnomad.utils.slack import slack_notifications
from gnomad.utils.file_utils import file_exists

from rmc.resources.basics import (
    LOGGING_PATH,
    SINGLE_BREAK_TEMP_PATH,
    TEMP_PATH_WITH_FAST_DEL,
    TEMP_PATH_WITH_SLOW_DEL,
)
from rmc.resources.gnomad import (
    constraint_ht,
    filtered_exomes,
    processed_exomes,
    prop_obs_coverage,
)
from rmc.resources.reference_data import filtered_context, gene_model
from rmc.resources.rmc import (
    CURRENT_FREEZE,
    constraint_prep,
    merged_search_ht_path,
    P_VALUE,
    rmc_results,
    simul_search_round_bucket_path,
    single_search_round_ht_path,
)
from rmc.resources.resource_utils import MISSENSE
from rmc.slack_creds import slack_token
from rmc.utils.constraint import (
    add_obs_annotation,
    calculate_exp_per_transcript,
    calculate_observed,
    check_break_search_round_nums,
    create_no_breaks_he,
    get_max_chisq_per_group,
    GROUPINGS,
    merge_rmc_hts,
    process_sections,
)
from rmc.utils.generic import (
    filter_context_using_gnomad,
    filter_to_region_type,
    generate_models,
    get_constraint_transcripts,
    get_coverage_correction_expr,
    keep_criteria,
    process_context_ht,
    process_vep,
)


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_missense_constraint")
logger.setLevel(logging.INFO)


def main(args):
    """Call functions from `constraint.py` to calculate regional missense constraint."""
    try:
        if args.pre_process_data:
            hl.init(
                log="/RMC_pre_process.log",
                tmp_dir=TEMP_PATH_WITH_FAST_DEL,
                quiet=args.quiet,
            )
            # TODO: Add code to create annotations necessary for constraint_flag_expr and filter transcripts prior to running constraint
            logger.warning("Code currently only processes b37 data!")

            logger.info("Preprocessing reference fasta (context) HT...")
            context_ht = process_context_ht()

            logger.info(
                "Filtering context HT to all covered sites not found or rare in gnomAD exomes"
            )
            context_ht = filter_context_using_gnomad(
                context_ht, "exomes", filter_context_using_cov=True
            )
            context_ht.write(filtered_context.path, overwrite=args.overwrite)

            logger.info(
                "Filtering gnomAD exomes HT to missense variants in canonical transcripts only..."
            )
            exome_ht = processed_exomes.ht()
            exome_ht = process_vep(exome_ht, filter_csq=True, csq=MISSENSE)

            # Move nested annotations into top level annotations
            exome_ht = exome_ht.select(
                ac=exome_ht.freq[0].AC,
                af=exome_ht.freq[0].AF,
                filters=exome_ht.filters,
                exome_coverage=exome_ht.coverage.exomes.median,
                transcript_consequences=exome_ht.transcript_consequences,
            )
            exome_ht = exome_ht.filter(
                keep_criteria(
                    ac_expr=exome_ht.ac,
                    af_expr=exome_ht.af,
                    filters_expr=exome_ht.filters,
                    cov_expr=exome_ht.exome_coverage,
                )
            )
            exome_ht.write(filtered_exomes.path, overwrite=args.overwrite)

            logger.info("Done preprocessing files")

        if args.prep_for_constraint:
            hl.init(
                log="/RMC_prep_for_constraint.log",
                tmp_dir=TEMP_PATH_WITH_FAST_DEL,
                quiet=args.quiet,
            )
            logger.info("Reading in exome HT...")
            exome_ht = filtered_exomes.ht()

            logger.info("Reading in context HT...")
            context_ht = filtered_context.ht()

            logger.info("Building plateau and coverage models...")
            coverage_ht = prop_obs_coverage.ht()
            coverage_x_ht = hl.read_table(
                prop_obs_coverage.path.replace(".ht", "_x.ht")
            )
            coverage_y_ht = hl.read_table(
                prop_obs_coverage.path.replace(".ht", "_y.ht")
            )

            (
                coverage_model,
                plateau_models,
                plateau_x_models,
                plateau_y_models,
            ) = generate_models(
                coverage_ht,
                coverage_x_ht,
                coverage_y_ht,
            )

            context_ht = context_ht.annotate_globals(
                plateau_models=plateau_models,
                plateau_x_models=plateau_x_models,
                plateau_y_models=plateau_y_models,
                coverage_model=coverage_model,
            )
            context_ht = context_ht.annotate(
                coverage_correction=get_coverage_correction_expr(
                    context_ht.exome_coverage,
                    coverage_model,
                    args.high_cov_cutoff,
                )
            )

            # TODO: Make this a separate section of code (with its own argument)
            if not args.skip_calc_oe:
                logger.info(
                    "Adding coverage correction to mutation rate probabilities..."
                )
                context_ht = context_ht.annotate(
                    raw_mu_snp=context_ht.mu_snp,
                    mu_snp=context_ht.mu_snp
                    * get_coverage_correction_expr(
                        context_ht.exome_coverage, context_ht.coverage_model
                    ),
                )

                logger.info(
                    "Creating autosomes + chrX PAR, chrX non-PAR, and chrY non-PAR HT versions..."
                )
                context_auto_ht = filter_to_region_type(context_ht, "autosomes")
                context_x_ht = filter_to_region_type(context_ht, "chrX")
                context_y_ht = filter_to_region_type(context_ht, "chrY")

                logger.info("Calculating expected values per transcript...")
                exp_ht = calculate_exp_per_transcript(
                    context_auto_ht,
                    locus_type="autosomes",
                    groupings=GROUPINGS,
                )
                exp_x_ht = calculate_exp_per_transcript(
                    context_x_ht,
                    locus_type="X",
                    groupings=GROUPINGS,
                )
                exp_y_ht = calculate_exp_per_transcript(
                    context_y_ht,
                    locus_type="Y",
                    groupings=GROUPINGS,
                )
                exp_ht = exp_ht.union(exp_x_ht).union(exp_y_ht)

                logger.info(
                    "Fixing expected values for genes that span PAR and nonPAR regions..."
                )
                # Adding a sum here to make sure that genes like XG that span PAR/nonPAR regions
                # have correct total expected values
                exp_ht = exp_ht.group_by(transcript=exp_ht.transcript).aggregate(
                    section_exp=hl.agg.sum(exp_ht.expected),
                    section_mu=hl.agg.sum(exp_ht.mu_agg),
                )
                # TODO: Write exp HT here

                logger.info(
                    "Aggregating total observed variant counts per transcript..."
                )
                obs_ht = calculate_observed(exome_ht)
                # TODO: Write obs HT here

            else:
                logger.warning(
                    "Using observed and expected values calculated using LoF pipeline..."
                )
                exp_ht = (
                    constraint_ht.ht()
                    .key_by("transcript")
                    .select("obs_mis", "exp_mis", "oe_mis")
                )

                # Filter to canonical transcripts only and rename fields
                exp_ht = exp_ht.filter(exp_ht.canonical)
                exp_ht = exp_ht.transmute(
                    observed=exp_ht.obs_mis,
                    expected=exp_ht.exp_mis,
                    mu_agg=exp_ht.mu_mis,
                )
                obs_ht = exp_ht

            logger.info(
                "Annotating context HT with number of observed and expected variants per site..."
            )
            # Add observed variants to context HT
            context_ht = add_obs_annotation(context_ht, filter_csq=True)

            logger.info(
                "Collecting by key to run constraint per base and not per base-allele..."
            )
            # Context HT is keyed by locus and allele, which means there is one row for every possible missense variant
            # This means that any locus could be present up to three times (once for each possible missense)
            # Collect by key here to ensure all loci are unique
            context_ht = context_ht.key_by("locus", "transcript").collect_by_key()
            context_ht = context_ht.annotate(
                # Collect the mutation rate probabilities at each locus
                mu_snp=hl.sum(context_ht.values.mu_snp),
                # Collect the observed counts for each locus
                # (this includes counts for each possible missense at the locus)
                observed=hl.sum(context_ht.values.observed),
                # Take just the first coverage value, since the locus should have the same coverage across the possible variants
                coverage=context_ht.values.exome_coverage[0],
            )
            # NOTE: v2 constraint_prep HT has total and cumulative values annotated
            # (HT as written before we decided to move these annotations to within
            # `process_sections`)
            # TODO: Update get_subsection_exprs to pull expected values from exp_ht
            context_ht = context_ht.drop("values")
            context_ht = context_ht.write(
                constraint_prep.path, overwrite=args.overwrite
            )

        if args.search_for_single_break:
            hl.init(
                log=f"/round{args.search_num}_single_break_search.log",
                tmp_dir=TEMP_PATH_WITH_FAST_DEL,
                quiet=args.quiet,
            )
            chisq_threshold = hl.eval(hl.qchisqtail(P_VALUE, 1))
            if args.p_value:
                chisq_threshold = hl.eval(hl.qchisqtail(args.p_value, 1))
            run_single_search = True

            logger.info(
                "Searching for transcripts or transcript subsections with a single significant break..."
            )
            if args.search_num == 1:

                all_loci_chisq_ht_path = f"{SINGLE_BREAK_TEMP_PATH}/all_loci_chisq.ht"
                if file_exists(all_loci_chisq_ht_path) and not args.save_chisq_ht:
                    ht = hl.read_table(all_loci_chisq_ht_path)
                    ht = get_max_chisq_per_group(ht, "section", "chisq", args.freeze)
                    ht = ht.annotate(
                        is_break=(
                            (ht.chisq == ht.section_max_chisq)
                            & (ht.chisq >= args.chisq_threshold)
                        )
                    )
                    run_single_search = False

                else:
                    # Read constraint_prep resource HT if this is the first search
                    ht = constraint_prep.ht()

                    logger.info(
                        "Adding section annotation before searching for first break..."
                    )
                    # Add transcript start and stop positions from browser HT
                    transcript_ht = gene_model.ht().select("start", "stop")
                    ht = ht.annotate(**transcript_ht[ht.transcript])
                    ht = ht.annotate(
                        section=hl.format("%s_%s_%s", ht.transcript, ht.start, ht.stop)
                    ).drop("start", "stop")
                    ht = ht.key_by("locus", "section").drop("transcript")

            else:
                # Read in merged single and simultaneous breaks results HT
                ht = hl.read_table(
                    merged_search_ht_path(
                        search_num=args.search_num - 1,
                        freeze=args.freeze,
                    ),
                )

            if run_single_search:
                logger.info(
                    "Calculating nulls, alts, and chi square values and checkpointing..."
                )
                ht = process_sections(
                    ht=ht,
                    search_num=args.search_num,
                    freeze=args.freeze,
                    chisq_threshold=chisq_threshold,
                    save_chisq_ht=args.save_chisq_ht,
                )
                ht = ht.checkpoint(
                    f"{TEMP_PATH_WITH_FAST_DEL}/freeze{args.freeze}_round{args.search_num}_temp.ht",
                    overwrite=True,
                )

            logger.info(
                "Extracting breakpoints found in round %i...",
                args.search_num,
            )
            breakpoint_ht = ht.filter(ht.is_break)
            breakpoint_ht = breakpoint_ht.annotate_globals(
                chisq_threshold=chisq_threshold
            )
            breakpoint_ht = breakpoint_ht.key_by("section")
            breakpoint_ht = breakpoint_ht.checkpoint(
                single_search_round_ht_path(
                    search_num=args.search_num,
                    is_break_found=True,
                    is_breakpoint_only=True,
                    freeze=args.freeze,
                ),
                overwrite=args.overwrite,
            )

            logger.info(
                "Filtering to transcripts or transcript subsections with breaks and checkpointing..."
            )
            ht = ht.annotate(breakpoint=breakpoint_ht[ht.section].locus.position)
            # Possible checkpoint here if necessary
            logger.info(
                "Splitting at breakpoints and re-annotating section starts, stops, and names..."
            )
            break_found_ht = ht.filter(hl.is_defined(ht.breakpoint))
            logger.info("Writing out sections with single significant break...")
            break_found_ht.write(
                single_search_round_ht_path(
                    search_num=args.search_num,
                    is_break_found=True,
                    is_breakpoint_only=False,
                    freeze=args.freeze,
                ),
                overwrite=args.overwrite,
            )

            logger.info(
                "Filtering HT to sections without a significant break and writing..."
            )
            no_break_found_ht = ht.filter(hl.is_missing(ht.breakpoint))
            no_break_found_ht = no_break_found_ht.drop("breakpoint")
            no_break_found_ht.write(
                single_search_round_ht_path(
                    search_num=args.search_num,
                    is_break_found=False,
                    is_breakpoint_only=False,
                    freeze=args.freeze,
                ),
                overwrite=args.overwrite,
            )

        if args.merge_single_simul:
            hl.init(
                log=f"/round{args.search_num}_merge_single_simul.log",
                tmp_dir=TEMP_PATH_WITH_FAST_DEL,
                quiet=args.quiet,
            )
            # Get locus input to this simultaneous break search round from single search no-break HT
            single_no_break_ht = hl.read_table(
                single_search_round_ht_path(
                    search_num=args.search_num,
                    is_break_found=False,
                    is_breakpoint_only=False,
                    freeze=args.freeze,
                )
            )

            logger.info("Checking if simul breaks merged HT exists...")
            simul_results_path = simul_search_round_bucket_path(
                search_num=args.search_num,
                bucket_type="final_results",
                freeze=args.freeze,
            )
            simul_break_by_section_path = f"{simul_results_path}/merged.ht"
            if file_exists(simul_break_by_section_path):
                simul_exists = True
                logger.info(
                    "Converting merged simultaneous breaks HT from section-level to locus-level..."
                )
                simul_break_by_section_ht = hl.read_table(simul_break_by_section_path)

                # Filter locus-level table (from single search no-break results) to sections with simultaneous breaks
                single_no_break_ht = single_no_break_ht.annotate(
                    breakpoints=simul_break_by_section_ht[
                        single_no_break_ht.section
                    ].breakpoints
                )
                simul_break_ht = single_no_break_ht.filter(
                    hl.is_defined(single_no_break_ht.breakpoints)
                )
                logger.info(
                    "Annotating simul breaks with new sections and re-keying for next search..."
                )
                simul_break_ht = simul_break_ht.annotate(
                    section_1=hl.if_else(
                        simul_break_ht.locus.position > simul_break_ht.breakpoints[1],
                        hl.format(
                            "%s_%s_%s",
                            simul_break_ht.section.split("_")[0],
                            simul_break_ht.breakpoints[1] + 1,
                            simul_break_ht.section.split("_")[2],
                        ),
                        hl.if_else(
                            simul_break_ht.locus.position
                            > simul_break_ht.breakpoints[0],
                            hl.format(
                                "%s_%s_%s",
                                simul_break_ht.section.split("_")[0],
                                simul_break_ht.breakpoints[0] + 1,
                                simul_break_ht.breakpoints[1],
                            ),
                            hl.format(
                                "%s_%s_%s",
                                simul_break_ht.section.split("_")[0],
                                simul_break_ht.section.split("_")[1],
                                simul_break_ht.breakpoints[0],
                            ),
                        ),
                    )
                )
                simul_break_ht = simul_break_ht.key_by(
                    "locus", section=simul_break_ht.section_1
                ).drop("section_1", "breakpoints")
            else:
                simul_exists = False
                logger.info(
                    "No sections in round %i had breakpoints in simultaneous breaks search.",
                    args.search_num,
                )

            single_break_path = single_search_round_ht_path(
                search_num=args.search_num,
                is_break_found=True,
                is_breakpoint_only=False,
                freeze=args.freeze,
            )

            if file_exists(single_break_path):
                single_exists = True
                logger.info(
                    "Annotating single breaks with new sections and re-keying for next search..."
                )
                single_break_ht = hl.read_table(single_break_path)
                single_break_ht = single_break_ht.annotate(
                    section_1=hl.if_else(
                        single_break_ht.locus.position > single_break_ht.breakpoint,
                        hl.format(
                            "%s_%s_%s",
                            single_break_ht.section.split("_")[0],
                            single_break_ht.breakpoint + 1,
                            single_break_ht.section.split("_")[2],
                        ),
                        hl.format(
                            "%s_%s_%s",
                            single_break_ht.section.split("_")[0],
                            single_break_ht.section.split("_")[1],
                            single_break_ht.breakpoint,
                        ),
                    )
                )
                single_break_ht = single_break_ht.key_by(
                    "locus", section=single_break_ht.section_1
                ).drop("section_1", "breakpoint")
            else:
                single_exists = False
                logger.info(
                    "No sections in round %i had breakpoints in single search.",
                    args.search_num,
                )

            merged_path = merged_search_ht_path(
                search_num=args.search_num,
                freeze=args.freeze,
            )
            if single_exists and simul_exists:
                logger.info(
                    "Merging break results from single and simultaneous search and writing..."
                )
                merged_break_ht = single_break_ht.union(simul_break_ht)
                # TODO: Change break results bucket structure to have round first, then simul vs. single split
                merged_break_ht.write(merged_path, overwrite=args.overwrite)
            elif single_exists:
                single_break_ht.write(merged_path, overwrite=args.overwrite)
            elif simul_exists:
                simul_break_ht.write(merged_path, overwrite=args.overwrite)
            else:
                logger.info(
                    "No sections in round %i had breakpoints (neither in single nor in simultaneous search).",
                    args.search_num,
                )

            # DONE: 1. Annotate this newly found simul_ht with same annotations as on break_found_ht
            # DONE: 2. Merge simul_ht with break_found_ht and write
            # TODO: 3. Add validity checks that we haven't dropped any transcripts/sections - e.g. using sections with breaks expression
            # DONE: 4. Create and write final no break found ht for this round number

        if args.finalize:
            hl.init(
                log="/RMC_finalize.log",
                tmp_dir=TEMP_PATH_WITH_FAST_DEL,
                quiet=args.quiet,
            )

            logger.info("Checking round paths...")
            round_nums = check_break_search_round_nums(args.freeze)

            logger.info("Finalizing section-level RMC table...")
            rmc_ht = merge_rmc_hts(round_nums=round_nums, freeze=args.freeze)
            rmc_ht = rmc_ht.checkpoint(
                f"{TEMP_PATH_WITH_SLOW_DEL}/freeze{args.freeze}_rmc_results.ht",
                overwrite=args.overwrite,
                _read_if_exists=not args.overwrite,
            )

            logger.info("Removing outlier transcripts...")
            constraint_transcripts = get_constraint_transcripts(outlier=False)
            rmc_ht = rmc_ht.filter(constraint_transcripts.contains(rmc_ht.transcript))

            logger.info("Writing out RMC results...")
            rmc_ht.write(rmc_results.versions[args.freeze].path)

            logger.info("Getting transcripts without evidence of RMC...")
            create_no_breaks_he(freeze=args.freeze, overwrite=args.overwrite)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(LOGGING_PATH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script searches for regional missense constraint in gnomAD."
    )
    parser.add_argument(
        "--pre-process-data", help="Pre-process data.", action="store_true"
    )
    parser.add_argument(
        "--n-partitions",
        help="Desired number of partitions for output data.",
        type=int,
        default=40000,
    )
    parser.add_argument(
        "--high-cov-cutoff",
        help="Coverage threshold for a site to be considered high coverage",
        type=int,
        default=40,
    )
    parser.add_argument(
        "--prep-for-constraint",
        help="Prepare tables for constraint calculations.",
        action="store_true",
    )
    parser.add_argument(
        "--skip-calc-oe",
        help="Skip observed and expected variant calculations per transcript. Relevant only to gnomAD v2.1.1!",
        action="store_true",
    )
    parser.add_argument(
        "--search-for-single-break",
        help="Search for single significant break in transcripts or transcript subsections.",
        action="store_true",
    )
    parser.add_argument(
        "--search-num",
        help="Search iteration number (e.g., second round of searching for single break would be 2).",
        type=int,
    )
    parser.add_argument(
        "--save-chisq-ht",
        help="""
        Save temporary Table that contains chi square significance values
        for all possible loci. Note that chi square values will be missing for
        any loci that would divide a transcript into subsections with fewer than
        `MIN_EXP_MIS` expected missense variants.

        NOTE that this temporary Table should only get saved once.
        """,
        action="store_true",
    )
    parser.add_argument(
        "--merge-single-simul",
        help="""
        Merge those transcripts/transcript sections with significant breakpoints from
        single and simultaneous breaks searches into a single Table, and those without
        significant breakpoints into another Table.
        """,
        action="store_true",
    )
    parser.add_argument(
        "--finalize",
        help="Combine and reformat (finalize) RMC output.",
        action="store_true",
    )
    parser.add_argument(
        "--p-value",
        help="""
        p-value significance threshold for single break search.
        Used to determine chi square statistic thershold.

        If not specified, script will default to threshold set
        in `P_VALUE`.
        """,
        type=float,
    )
    parser.add_argument(
        "--freeze",
        help="RMC data freeze number",
        default=CURRENT_FREEZE,
    )
    parser.add_argument(
        "--overwrite", help="Overwrite existing data.", action="store_true"
    )
    parser.add_argument(
        "--slack-channel",
        help="Send message to Slack channel/user.",
    )
    parser.add_argument(
        "--quiet",
        help="Initialize Hail with `quiet=True` to print fewer log messages",
        action="store_true",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
