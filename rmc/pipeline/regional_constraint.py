import argparse
import logging

import hail as hl
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import (
    LOGGING_PATH,
    SINGLE_BREAK_TEMP_PATH,
    TEMP_PATH_WITH_FAST_DEL,
    TEMP_PATH_WITH_SLOW_DEL,
)
from rmc.resources.gnomad import (
    filtered_exomes,
    processed_exomes,
    prop_obs_coverage,
)
from rmc.resources.reference_data import filtered_context, gene_model
from rmc.resources.resource_utils import MISSENSE, NONSENSE, SYNONYMOUS
from rmc.resources.rmc import (
    CURRENT_FREEZE,
    P_VALUE,
    constraint_prep,
    merged_search_ht_path,
    rmc_results,
    simul_search_round_bucket_path,
    single_search_round_ht_path,
)
from rmc.slack_creds import slack_token
from rmc.utils.constraint import (
    add_obs_annotation,
    calculate_exp_from_mu,
    check_break_search_round_nums,
    create_context_with_oe,
    create_no_breaks_he,
    annotate_max_chisq_per_section,
    merge_rmc_hts,
    process_sections,
)
from rmc.utils.generic import (
    filter_context_using_gnomad,
    filter_to_region_type,
    generate_models,
    get_constraint_transcripts,
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
            logger.warning("Code currently only processes b37 data!")

            logger.info(
                "Preprocessing VEP context HT to filter to missense, nonsense, and"
                " synonymous variants in all canonical transcripts and add constraint"
                " annotations..."
            )
            # NOTE: Constraint outlier transcripts are not removed before computing RMC
            context_ht = process_context_ht(
                filter_csq=True, csq={MISSENSE, NONSENSE, SYNONYMOUS}
            )
            logger.info("Checkpointing context HT...")
            context_ht = context_ht.checkpoint(
                f"{TEMP_PATH_WITH_FAST_DEL}/processed_context.ht"
            )
            logger.info(
                "Filtering context HT to all covered sites not found or rare in gnomAD"
                " exomes..."
            )
            context_ht = filter_context_using_gnomad(context_ht, "exomes")
            logger.info("Writing out context HT...")
            context_ht.write(filtered_context.path, overwrite=args.overwrite)

            logger.info(
                "Filtering gnomAD exomes HT to missense, nonsense, and synonymous"
                " variants in all canonical transcripts..."
            )
            exome_ht = processed_exomes.ht()
            exome_ht = process_vep(
                exome_ht, filter_csq=True, csq={MISSENSE, NONSENSE, SYNONYMOUS}
            )
            logger.info("Checkpointing gnomAD exomes HT...")
            exome_ht = exome_ht.checkpoint(
                f"{TEMP_PATH_WITH_FAST_DEL}/processed_vep.ht"
            )
            logger.info(
                "Filtering gnomAD exomes HT to rare variants that pass filters and"
                " coverage criteria..."
            )
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
            logger.info("Writing out filtered gnomAD exomes HT...")
            exome_ht.write(filtered_exomes.path, overwrite=args.overwrite)

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

            logger.info("Calculating expected values per allele...")
            context_ht = (
                calculate_exp_from_mu(
                    filter_to_region_type(context_ht, "autosomes"),
                    locus_type="autosomes",
                )
                .union(
                    calculate_exp_from_mu(
                        filter_to_region_type(context_ht, "chrX"), locus_type="X"
                    )
                )
                .union(
                    calculate_exp_from_mu(
                        filter_to_region_type(context_ht, "chrY"), locus_type="Y"
                    )
                )
            )
            # TODO: Remove sites where expected is negative
            context_ht = context_ht.checkpoint(
                f"{TEMP_PATH_WITH_FAST_DEL}/context_exp.ht"
            )

            logger.info(
                "Annotating context HT with number of observed variants and writing"
                " out..."
            )
            context_ht = add_obs_annotation(context_ht)
            context_ht = context_ht.write(
                constraint_prep.path, overwrite=args.overwrite
            )
            # TODO: Repartition?

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
                "Searching for transcripts or transcript subsections with a single"
                " significant break..."
            )
            if args.search_num == 1:
                # TODO: Move existing HT to a different path so this can be remade
                all_loci_chisq_ht_path = f"{SINGLE_BREAK_TEMP_PATH}/all_loci_chisq.ht"
                if file_exists(all_loci_chisq_ht_path) and not args.save_chisq_ht:
                    logger.info("Reading in all loci chisq HT...")
                    ht = hl.read_table(all_loci_chisq_ht_path)
                    ht = annotate_max_chisq_per_section(ht, args.freeze)
                    ht = ht.annotate(
                        is_break=(
                            (ht.chisq == ht.section_max_chisq)
                            & (ht.chisq >= chisq_threshold)
                        )
                    )
                    run_single_search = False

                else:
                    # Read `constraint_prep` resource HT if this is the first search
                    logger.info(
                        "Reading in constraint prep HT, filtering to missenses, and"
                        " aggregating by locus..."
                    )
                    ht = constraint_prep.ht()
                    ht = ht.filter(ht.annotation == MISSENSE)
                    # Context HT is keyed by locus and allele, which means there is one row for every possible missense variant
                    # This means that any locus could be present up to three times (once for each possible missense)
                    ht = ht.group_by("locus", "transcript").aggregate(
                        mu_snp=hl.sum(ht.mu_snp),
                        observed=hl.sum(ht.observed),
                        expected=hl.sum(ht.expected),
                        # Locus should have the same coverage across the possible variants
                        coverage=ht.coverage[0],
                    )

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
                logger.info(
                    "Reading in merged single and simultaneous breaks results HT..."
                )
                ht = hl.read_table(
                    merged_search_ht_path(
                        search_num=args.search_num - 1,
                        freeze=args.freeze,
                    ),
                )

            if run_single_search:
                logger.info(
                    "Calculating nulls, alts, and chi square values and"
                    " checkpointing..."
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
                "Filtering to transcripts or transcript subsections with breaks and"
                " checkpointing..."
            )
            ht = ht.annotate(breakpoint=breakpoint_ht[ht.section].locus.position)
            # Possible checkpoint here if necessary
            logger.info(
                "Splitting at breakpoints and re-annotating section starts, stops, and"
                " names..."
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
                    "Converting merged simultaneous breaks HT from section-level to"
                    " locus-level..."
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
                    "Annotating simul breaks with new sections and re-keying for next"
                    " search..."
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
                    (
                        "No sections in round %i had breakpoints in simultaneous breaks"
                        " search."
                    ),
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
                    "Annotating single breaks with new sections and re-keying for next"
                    " search..."
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
                    "Merging break results from single and simultaneous search and"
                    " writing..."
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
                    (
                        "No sections in round %i had breakpoints (neither in single nor"
                        " in simultaneous search)."
                    ),
                    args.search_num,
                )

            # TODO: Add validity checks that we haven't dropped any transcripts/sections - e.g. using sections with breaks expression

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
            # Drop any unnecessary global fields
            rmc_ht = rmc_ht.select_globals()
            rmc_ht = rmc_ht.checkpoint(
                f"{TEMP_PATH_WITH_SLOW_DEL}/freeze{args.freeze}_rmc_results.ht",
                overwrite=args.overwrite,
                _read_if_exists=not args.overwrite,
            )

            # TODO: Remove this block to retain outlier transcripts initially
            logger.info("Removing outlier transcripts...")
            constraint_transcripts = get_constraint_transcripts(outlier=False)
            rmc_ht = rmc_ht.filter(constraint_transcripts.contains(rmc_ht.transcript))

            # Add p-value threshold to globals
            if not args.p_value:
                logger.warning(
                    "p-value threshold not specified! Defaulting to value stored in"
                    " `P_VALUE` constant..."
                )
                p_value = P_VALUE
            else:
                p_value = args.p_value
            rmc_ht = rmc_ht.annotate_globals(p_value=p_value)

            logger.info("Writing out RMC results...")
            rmc_ht.write(
                rmc_results.versions[args.freeze].path, overwrite=args.overwrite
            )

            logger.info("Getting transcripts without evidence of RMC...")
            create_no_breaks_he(freeze=args.freeze, overwrite=args.overwrite)

            logger.info("Creating OE-annotated context table...")
            create_context_with_oe(
                freeze=args.freeze, overwrite_temp=args.overwrite_temp
            )

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
        help=(
            "Skip observed and expected variant calculations per transcript. Relevant"
            " only to gnomAD v2.1.1!"
        ),
        action="store_true",
    )
    parser.add_argument(
        "--search-for-single-break",
        help=(
            "Search for single significant break in transcripts or transcript"
            " subsections."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--search-num",
        help=(
            "Search iteration number (e.g., second round of searching for single break"
            " would be 2)."
        ),
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
        Used to determine chi square threshold for likelihood ratio test.

        Also used to annotate globals of final RMC HT with p-value associated
        with RMC results.

        If not specified, script will default to threshold set
        in `P_VALUE`.
        """,
        type=float,
    )
    parser.add_argument(
        "--freeze",
        help="RMC data freeze number",
        type=int,
        default=CURRENT_FREEZE,
    )
    parser.add_argument(
        "--overwrite",
        help="""
        Overwrite existing output data.
        Applies to all outputs except OE-annotated context table (created in `--finalize`).
        """,
        action="store_true",
    )
    parser.add_argument(
        "--overwrite-temp",
        help="""
        Overwrite existing intermediate temporary data.
        Only applicable in creating OE-annotated context table.
        """,
        action="store_true",
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
