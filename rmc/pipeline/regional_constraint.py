import argparse
import logging

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import (
    LOGGING_PATH,
    TEMP_PATH,
    TEMP_PATH_WITH_DEL,
)
from rmc.resources.gnomad import (
    constraint_ht,
    filtered_exomes,
    processed_exomes,
    prop_obs_coverage,
)
from rmc.resources.reference_data import filtered_context, gene_model
from rmc.resources.rmc import (
    constraint_prep,
    merged_search_ht_path,
    multiple_breaks,
    simul_break,
    simul_search_round_bucket_path,
    single_search_round_ht_path,
)
from rmc.resources.resource_utils import MISSENSE
from rmc.slack_creds import slack_token
from rmc.utils.constraint import (
    add_obs_annotation,
    calculate_exp_per_transcript,
    calculate_observed,
    get_rescue_transcripts_and_create_no_breaks_ht,
    GROUPINGS,
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
            hl.init(log="/RMC_pre_process.log", tmp_dir=TEMP_PATH_WITH_DEL)
            # TODO: Add code to create annotations necessary for constraint_flag_expr and filter transcripts prior to running constraint
            logger.warning("Code currently only processes b37 data!")

            logger.info("Preprocessing reference fasta (context) HT...")
            context_ht = process_context_ht(args.trimers)

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
                pass_filters=exome_ht.pass_filters,
                exome_coverage=exome_ht.coverage.exomes.median,
                transcript_consequences=exome_ht.transcript_consequences,
            )
            exome_ht = exome_ht.filter(keep_criteria(exome_ht))
            exome_ht.write(filtered_exomes.path, overwrite=args.overwrite)

            logger.info("Done preprocessing files")

        if args.prep_for_constraint:
            hl.init(log="/RMC_prep_for_constraint.log", tmp_dir=TEMP_PATH_WITH_DEL)
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
                coverage_ht, coverage_x_ht, coverage_y_ht, trimers=args.trimers
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
            is_rescue = args.is_rescue
            hl.init(
                log=f"/round{args.search_num}_single_break_search.log",
                tmp_dir=TEMP_PATH_WITH_DEL,
            )

            logger.info(
                "Searching for transcripts or transcript subsections with a single significant break..."
            )
            if args.search_num == 1:
                assert not is_rescue, "No code to support rescue search round 1!"

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
                        is_rescue=is_rescue,
                        search_num=args.search_num - 1,
                        is_break_found=True,
                    ),
                )

            logger.info(
                "Calculating nulls, alts, and chi square values and checkpointing..."
            )
            ht = process_sections(
                ht=ht,
                search_num=args.search_num,
                chisq_threshold=args.chisq_threshold,
            )
            ht = ht.checkpoint(
                f"{TEMP_PATH_WITH_DEL}/round{args.search_num}_temp.ht", overwrite=True
            )

            logger.info(
                "Extracting breakpoints found in round %i...",
                args.search_num,
            )
            breakpoint_ht = ht.filter(ht.is_break)
            breakpoint_ht = breakpoint_ht.annotate_globals(
                chisq_threshold=args.chisq_threshold
            )
            breakpoint_ht = breakpoint_ht.key_by("section")
            breakpoint_ht = breakpoint_ht.checkpoint(
                single_search_round_ht_path(
                    is_rescue=is_rescue,
                    search_num=args.search_num,
                    is_break_found=True,
                    is_breakpoint_only=True,
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
                    is_rescue=is_rescue,
                    search_num=args.search_num,
                    is_break_found=True,
                    is_breakpoint_only=False,
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
                    is_rescue=is_rescue,
                    search_num=args.search_num,
                    is_break_found=False,
                    is_breakpoint_only=False,
                ),
                overwrite=args.overwrite,
            )

        if args.merge_single_simul:
            if args.is_rescue and args.search_num == 1:
                get_rescue_transcripts_and_create_no_breaks_ht(args.overwrite)
                return

            logger.info(
                "Converting merged simultaneous breaks HT from section-level to locus-level..."
            )
            simul_results_path = simul_search_round_bucket_path(
                is_rescue=args.is_rescue,
                search_num=args.search_num,
                bucket_type="final_results",
            )
            simul_by_section_ht = hl.read_table(f"{simul_results_path}/merged.ht")
            # Use no_break_found HT from single breaks results to get locus input to simultaneous break search
            single_no_break_ht = hl.read_table(
                single_search_round_ht_path(
                    is_rescue=args.is_rescue,
                    search_num=args.search_num,
                    is_break_found=False,
                    is_breakpoint_only=False,
                )
            )
            # Filter to sections with simultaneous breaks
            single_no_break_ht = single_no_break_ht.annotate(
                breakpoints=simul_by_section_ht[single_no_break_ht.section].breakpoints
            )
            simul_break_ht = single_no_break_ht.filter(
                hl.is_defined(single_no_break_ht.breakpoints)
            )

            logger.info("Annotating new sections and re-keying for next search...")
            single_break_ht = hl.read_table(
                single_search_round_ht_path(
                    is_rescue=args.is_rescue,
                    search_num=args.search_num,
                    is_break_found=True,
                    is_breakpoint_only=False,
                )
            )
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
                        simul_break_ht.locus.position > simul_break_ht.breakpoints[0],
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

            logger.info(
                "Merging break results from single and simultaneous search and writing..."
            )
            merged_break_ht = single_break_ht.union(simul_break_ht)
            merged_break_ht.write(
                # TODO: Change break results bucket structure to have round first, then simul vs. single split
                merged_search_ht_path(
                    is_rescue=args.is_rescue,
                    search_num=args.search_num,
                    is_break_found=True,
                ),
                overwrite=args.overwrite,
            )

            logger.info(
                "Merging no-break results from single and simultaneous search and writing..."
            )
            merged_no_break_ht = single_no_break_ht.filter(
                hl.is_missing(single_no_break_ht.breakpoints)
            ).drop("breakpoints")
            merged_no_break_ht.write(
                merged_search_ht_path(
                    is_rescue=args.is_rescue,
                    search_num=args.search_num,
                    is_break_found=False,
                ),
                overwrite=args.overwrite,
            )

            # DONE: 1. Annotate this newly found simul_ht with same annotations as on break_found_ht
            # DONE: 2. Merge simul_ht with break_found_ht and write
            # TODO: 3. Add validity checks that we haven't dropped any transcripts/sections - e.g. using sections with breaks expression
            # DONE: 4. Create and write final no break found ht for this round number

        if args.finalize:
            hl.init(log="/RMC_finalize.log", tmp_dir=TEMP_PATH_WITH_DEL)

            logger.info(
                "Getting start and end positions and total size for each transcript..."
            )
            if not file_exists(f"{TEMP_PATH}/transcript.ht"):
                raise DataException(
                    "Transcript HT doesn't exist. Please double check and recreate!"
                )

            if args.remove_outlier_transcripts:
                outlier_transcripts = get_constraint_transcripts(outlier=True)

            logger.info("Reading in context HT...")
            # Drop extra annotations from context HT
            context_ht = filtered_context.ht().drop(
                "obs_scan", "mu_scan", "forward_oe", "values"
            )
            context_ht = context_ht.filter(
                ~outlier_transcripts.contains(context_ht.transcript)
            )

            logger.info("Reading in results HTs...")
            multiple_breaks_ht = multiple_breaks.ht()
            multiple_breaks_ht = multiple_breaks_ht.filter(
                ~outlier_transcripts.contains(multiple_breaks_ht.transcript)
            )
            # TODO: Update simul breaks reformatting
            simul_breaks_ht = simul_break.ht()
            simul_breaks_ht = simul_breaks_ht.filter(
                ~outlier_transcripts.contains(simul_breaks_ht.transcript)
            )

            logger.info("Getting simultaneous breaks transcripts...")
            simul_break_transcripts = simul_breaks_ht.aggregate(
                hl.agg.collect_as_set(simul_breaks_ht.transcript), _localize=False
            )
            rmc_transcripts = simul_break_transcripts

            logger.info("Getting transcript information for breaks HT...")
            global_fields = multiple_breaks_ht.globals
            n_breaks = [
                int(x.split("_")[1]) for x in global_fields if "break" in x
            ].sort()
            rmc_transcripts = []
            for break_num in n_breaks:
                ht = hl.read_table(f"{TEMP_PATH}/break_{break_num}.ht")
                ht = ht.filter(ht.is_break)
                rmc_transcripts.append(
                    ht.aggregate(hl.agg.collect_as_set(ht.transcript), _localize=False)
                )

            logger.info("Removing overlapping transcript information...")
            # Extracting the unique set of transcripts for each break number
            # e.g., keeping only transcripts with a single break
            # and not transcripts also that had additional breaks in "break_1_transcripts"
            filtered_transcripts = {}

            logger.info(
                "Cycling through each set of transcripts with break information..."
            )
            for index, transcripts in enumerate(rmc_transcripts):
                # Cycle through the sets of transcripts for every number of breaks larger than current number of breaks
                for t in rmc_transcripts[index + 1 :]:
                    transcripts = transcripts.difference(t)

                # Separate transcripts with a single break from transcripts with multiple breaks
                # This is because transcripts with a single break are annotated into a separate struct
                # in the release HT globals
                if index == 0:
                    one_break_transcripts = transcripts
                else:
                    filtered_transcripts[f"break_{index + 1}_transcripts"] = transcripts

            # Getting total number of transcripts with evidence of rmc
            total_rmc_transcripts = simul_break_transcripts
            for break_num in filtered_transcripts:
                total_rmc_transcripts = total_rmc_transcripts.union(
                    filtered_transcripts[break_num]
                )
            logger.info(
                "Number of transcripts with evidence of RMC: %s",
                hl.eval(hl.len(total_rmc_transcripts)),
            )

            logger.info("Creating breaks HT...")
            breaks_ht = context_ht.filter(
                total_rmc_transcripts.contains(context_ht.transcript)
            )

            logger.info("Adding simultaneous breaks information to breaks HT...")
            breaks_ht = breaks_ht.annotate(
                simul_break_info=hl.struct(
                    is_break=simul_breaks_ht[breaks_ht.key].is_break,
                    window_end=simul_breaks_ht[breaks_ht.key].window_end,
                    post_window_pos=simul_breaks_ht[breaks_ht.key].post_window_pos,
                )
            )
            logger.info("Adding transcript information to breaks HT globals...")
            breaks_ht = breaks_ht.annotate_globals(
                single_break=hl.struct(transcripts=one_break_transcripts),
                multiple_breaks=hl.struct(**filtered_transcripts),
                simul_breaks=hl.struct(transcripts=simul_break_transcripts),
            )

            logger.info("Creating no breaks HT...")
            no_breaks_ht = context_ht.filter(
                ~total_rmc_transcripts.contains(context_ht.transcript)
            )

            # TODO: Add section chisq calculation here
            if (breaks_ht.count() + no_breaks_ht.count()) != context_ht.count():
                raise DataException(
                    "Row counts for breaks HT (one break, multiple breaks, simul breaks) and no breaks HT doesn't match context HT row count!"
                )

            logger.info("Checkpointing HTs...")
            breaks_ht = breaks_ht.checkpoint(f"{TEMP_PATH}/breaks.ht", overwrite=True)
            no_breaks_ht = no_breaks_ht.checkpoint(
                f"{TEMP_PATH}/no_breaks.ht", overwrite=True
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
        "--trimers", help="Use trimers instead of heptamers.", action="store_true"
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
        "--is-rescue",
        help="""
        Whether search is part of the 'rescue' pathway (pathway
        with lower chi square significance cutoff).
        """,
        action="store_true",
    )
    parser.add_argument(
        "--search-num",
        help="Search iteration number (e.g., second round of searching for single break would be 2).",
        type=int,
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
        "--chisq-threshold",
        help="Chi-square significance threshold. Value should be 6.6 (single break) or 9.2 (two breaks) (p = 0.01).",
        type=float,
        default=6.6,
    )
    parser.add_argument(
        "--overwrite", help="Overwrite existing data.", action="store_true"
    )
    parser.add_argument(
        "--slack-channel",
        help="Send message to Slack channel/user.",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
