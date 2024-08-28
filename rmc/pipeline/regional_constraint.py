import argparse
import logging

import hail as hl
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import (
    CONSTRAINT_PREFIX,
    LOGGING_PATH,
    TEMP_PATH_WITH_FAST_DEL,
    TEMP_PATH_WITH_SLOW_DEL,
)
from rmc.resources.resource_utils import (
    CURRENT_GNOMAD_VERSION,
    MISSENSE,
    NONSENSES,
    SYNONYMOUS,
)
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
    annotate_max_chisq_per_section,
    check_break_search_round_nums,
    create_constraint_prep_ht,
    create_context_with_oe,
    create_filtered_context_ht,
    create_no_breaks_he,
    create_rmc_release_downloads,
    format_rmc_browser_ht,
    merge_rmc_hts,
    process_sections,
    validate_rmc_release_downloads,
)
from rmc.utils.data_loading import create_transcript_ref
from rmc.utils.generic import get_constraint_transcripts

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_missense_constraint")
logger.setLevel(logging.INFO)


def main(args):
    """Call functions from `constraint.py` to calculate regional missense constraint."""
    try:
        if args.command == "create-transcript-refs":
            hl.init(
                log="/RMC_create_transcript_refs.log",
                tmp_dir=TEMP_PATH_WITH_FAST_DEL,
                quiet=args.quiet,
            )
            hl.default_reference("GRCh38")
            logger.info("Creating transcript reference resources...")
            create_transcript_ref(build="GRCh38", overwrite=args.overwrite)

        if args.command == "prep-filtered-context":
            hl.init(
                log="/RMC_pre_process.log",
                tmp_dir=TEMP_PATH_WITH_FAST_DEL,
                quiet=args.quiet,
            )
            hl.default_reference("GRCh38")
            logger.info("Creating filtered context HT...")
            n_partitions = 10000
            create_filtered_context_ht(
                canonical_only=args.filter_to_canonical,
                n_partitions=args.n_partitions if args.n_partitions else n_partitions,
                overwrite=args.overwrite_temp,
            )

        if args.command == "prep-constraint":
            hl.init(
                log="/RMC_prep_for_constraint.log",
                tmp_dir=TEMP_PATH_WITH_FAST_DEL,
                quiet=args.quiet,
            )
            hl.default_reference("GRCh38")
            logger.info("Creating constraint prep HT...")
            # Input to constraint prep step is filtered context HT created above
            # Constraint prep HT is filtered to missense variants by default
            # Use these args to run constraint prep on other variant consequences
            csq = {MISSENSE}
            if args.prep_nonsense:
                csq = NONSENSES
                logger.warning(
                    "Pipeline is currently set up to run on missenses. Please make sure"
                    " all missense-relevant files are deleted before running."
                )
            elif args.prep_synonymous:
                csq = SYNONYMOUS
                logger.warning(
                    "Pipeline is currently set up to run on missenses. Please make sure"
                    " all missense-relevant files are deleted before running."
                )
            n_partitions = 10000
            create_constraint_prep_ht(
                filter_csq=csq,
                n_partitions=args.n_partitions if args.n_partitions else n_partitions,
                overwrite=args.overwrite,
            )

        if args.command == "search-for-single-break":
            hl.init(
                log=f"/round{args.search_num}_single_break_search.log",
                tmp_dir=TEMP_PATH_WITH_FAST_DEL,
                quiet=args.quiet,
            )
            hl.default_reference("GRCh38")
            chisq_threshold = hl.eval(hl.qchisqtail(P_VALUE, 1))
            if args.p_value:
                chisq_threshold = hl.eval(hl.qchisqtail(args.p_value, 1))
            run_single_search = True

            logger.info(
                "Searching for transcripts or transcript subsections with a single"
                " significant break..."
            )
            if args.search_num == 1:
                all_loci_chisq_ht_path = f"{CONSTRAINT_PREFIX}/{CURRENT_GNOMAD_VERSION}/{args.freeze}/constraint_prep.ht"
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
                    logger.info("Reading in constraint prep HT...")
                    if args.n_partitions:
                        ht = hl.read_table(
                            constraint_prep.path,
                            _n_partitions=args.n_partitions,
                        )
                    else:
                        ht = constraint_prep.ht()

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
                _read_if_exists=not args.overwrite,
                overwrite=args.overwrite,
            )

            logger.info(
                "Filtering to transcripts or transcript subsections with breaks and"
                " checkpointing..."
            )
            ht = ht.annotate(breakpoint=breakpoint_ht[ht.section].locus.position)

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

        if args.command == "merge-single-simul":
            hl.init(
                log=f"/round{args.search_num}_merge_single_simul.log",
                tmp_dir=TEMP_PATH_WITH_FAST_DEL,
                quiet=args.quiet,
            )
            hl.default_reference("GRCh38")
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
            # NOTE: Not checking for row count of zero because am assuming
            # user didn't run merge step (`merge_hts.py`) if simul breaks search
            # didn't return any results
            simul_results_path = simul_search_round_bucket_path(
                search_num=args.search_num,
                bucket_type="final_results",
                freeze=args.freeze,
            )
            simul_break_by_section_path = f"{simul_results_path}/merged.ht"
            simul_exists = file_exists(simul_break_by_section_path)
            if simul_exists:
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
                logger.info(
                    "No sections in round %i had breakpoints in simultaneous breaks"
                    " search.",
                    args.search_num,
                )

            single_break_path = single_search_round_ht_path(
                search_num=args.search_num,
                is_break_found=True,
                is_breakpoint_only=False,
                freeze=args.freeze,
            )

            single_exists = file_exists(single_break_path)
            if single_exists:
                logger.info(
                    "Annotating single breaks with new sections and re-keying for next"
                    " search..."
                )
                single_break_ht = hl.read_table(single_break_path)
                if single_break_ht.count() > 0:
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
                    logger.info("Single break HT had zero rows!")
                    single_exists = False
            else:
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
                merged_break_ht.write(merged_path, overwrite=args.overwrite)
            elif single_exists:
                single_break_ht.write(merged_path, overwrite=args.overwrite)
            elif simul_exists:
                simul_break_ht.write(merged_path, overwrite=args.overwrite)
            else:
                logger.info(
                    "No sections in round %i had breakpoints (neither in single nor"
                    " in simultaneous search).",
                    args.search_num,
                )

            # TODO: Add validity checks that we haven't dropped any transcripts/sections - e.g. using sections with breaks expression

        if args.command == "finalize":
            hl.init(
                log="/RMC_finalize.log",
                tmp_dir=TEMP_PATH_WITH_FAST_DEL,
                quiet=args.quiet,
            )
            hl.default_reference("GRCh38")
            logger.info("Checking round paths...")
            round_nums = check_break_search_round_nums(args.freeze)

            logger.info("Finalizing section-level RMC table...")
            rmc_ht = merge_rmc_hts(
                round_nums=round_nums,
                freeze=args.freeze,
                overwrite_temp=args.overwrite_temp,
            )
            # Drop any unnecessary global fields
            rmc_ht = rmc_ht.select_globals()
            rmc_ht = rmc_ht.checkpoint(
                f"{TEMP_PATH_WITH_SLOW_DEL}/freeze{args.freeze}_rmc_results.ht",
                _read_if_exists=not args.overwrite,
                overwrite=args.overwrite,
            )

            if args.filter_to_canonical:
                logger.warning("Filtering to canonical transcripts only!")
                # NOTE: RMC search should be run on only canonical transcripts
                # rather than filtering to canonical transcripts at this step
                # for compute efficiency
                canonical_transcripts = get_constraint_transcripts(
                    all_transcripts=True,
                    filter_to_canonical=True,
                )
                rmc_ht = rmc_ht.filter(
                    canonical_transcripts.contains(rmc_ht.transcript)
                )

            if args.filter_outliers:
                logger.info("Removing outlier transcripts...")
                constraint_transcripts = get_constraint_transcripts(outlier=False)
                rmc_ht = rmc_ht.filter(
                    constraint_transcripts.contains(rmc_ht.transcript)
                )

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
            n_partitions = args.n_partitions if args.n_partitions else 1000
            rmc_ht = rmc_ht.naive_coalesce(n_partitions)
            rmc_ht.write(
                rmc_results.versions[args.freeze].path, overwrite=args.overwrite
            )

            logger.info("Getting transcripts without evidence of RMC...")
            create_no_breaks_he(freeze=args.freeze, overwrite=args.overwrite)
            # TODO: Create region-level table that combines `rmc_results` and `no_breaks_he`
            # (used in RMC assessment plots)

            logger.info("Reformatting RMC results for browser release...")
            # NOTE: `filter_to_canonical` is included here only to make
            # amino acid annotation with context HT more efficient if
            # RMC results were filtered to canonical transcripts
            format_rmc_browser_ht(
                args.freeze, args.overwrite_temp, args.filter_to_canonical
            )

        if args.command == "create-context-with-oe":
            hl.init(
                log="/RMC_create_context_with_oe.log",
                tmp_dir=TEMP_PATH_WITH_FAST_DEL,
                quiet=args.quiet,
            )
            hl.default_reference("GRCh38")
            logger.info("Creating OE-annotated context table...")
            create_context_with_oe(
                freeze=args.freeze,
                filter_to_canonical=args.filter_to_canonical,
                overwrite_temp=args.overwrite_temp,
                filter_outliers=args.filter_outliers,
            )

        if args.command == "create-rmc-release":
            logger.info(
                "Creating versions of files to be publicly released on gnomAD"
                " browser..."
            )
            create_rmc_release_downloads(args.freeze, args.overwrite)

        if args.command == "validate-rmc-release":
            logger.info("Running validity checks on public RMC files...")
            validate_rmc_release_downloads(args.freeze)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(LOGGING_PATH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script searches for regional missense constraint in gnomAD."
    )
    parser.add_argument(
        "--n-partitions",
        help="""
        Desired number of partitions for context Table or constraint prep Table depending on
        which is being prepared.
        """,
        type=int,
    )
    parser.add_argument(
        "--search-num",
        help="""
        Search iteration number (e.g., second round of searching for single break would be 2).
        """,
        type=int,
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
        Applies to all outputs except OE-annotated context table (created in `finalize`).
        """,
        action="store_true",
    )
    parser.add_argument(
        "--overwrite-temp",
        help="Overwrite existing intermediate temporary data.",
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
    parser.add_argument(
        "--filter-outliers",
        help="Remove constraint outlier transcripts from output.",
        action="store_true",
    )
    parser.add_argument(
        "--filter-to-canonical",
        help="Filter to canonical transcripts only.",
        action="store_true",
    )

    subparsers = parser.add_subparsers(title="command", dest="command", required=True)

    create_transcript_refs = subparsers.add_parser(
        "create-transcript-refs",
        help="Create transcript reference resources.",
    )

    prep_filtered_context = subparsers.add_parser(
        "prep-filtered-context",
        help="""
        Process VEP context HT to create allele-level filtered context HT with constraint annotations.
        """,
    )

    prep_constraint = subparsers.add_parser(
        "prep-constraint",
        help="""
        Prepare locus-level Table filtered to specific variant consequences for regional constraint
        break search.
        """,
    )
    prep_csq = parser.add_mutually_exclusive_group()
    prep_csq.add_argument(
        "--prep-nonsense",
        help="""
        "Filter to nonsense instead of missense variants in constraint prep Table.
        """,
        action="store_true",
    )
    prep_csq.add_argument(
        "--prep-synonymous",
        help="""
        "Filter to synonymous instead of missense variants in constraint prep Table.
        """,
        action="store_true",
    )

    single_break = subparsers.add_parser(
        "search-for-single-break",
        help="""
        Search for single significant break in transcripts or transcript subsections.
        """,
    )
    single_break.add_argument(
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

    merge_single_simul = subparsers.add_parser(
        "merge-single-simul",
        help="""
        Merge those transcripts/transcript sections with significant breakpoints from
        single and simultaneous breaks searches into a single Table, and those without
        significant breakpoints into another Table.
        """,
    )

    finalize = subparsers.add_parser(
        "finalize", help="Combine and reformat (finalize) RMC output."
    )

    create_oe_context = subparsers.add_parser(
        "create-context-with-oe",
        help="Create context Table with observed/expected values for RMC assessment.",
    )

    create_release = subparsers.add_parser(
        "create-rmc-release",
        help="Create RMC release files (to be publicly shared on gnomAD browser).",
    )
    validate_release = subparsers.add_parser(
        "validate-rmc-release",
        help="Check validity of RMC release files.",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
