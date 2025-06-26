"""
This script searches for two simultaneous breaks in groups of transcripts using Hail Query within Google Cloud Dataproc.

This script should be run only on transcripts that are greater than or equal
to the --transcript-len-threshold specified in `prepare_transcripts.py` if these transcripts are too slow or getting preempted in Hail Batch.
Transcripts smaller than --transcript-len-threshold  should be run using `run_batches.py` as they run quickly and inexpensively in Hail Batch.
"""
import argparse
import logging

import hail as hl
from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists

from rmc.resources.basics import LOGGING_PATH, TEMP_PATH_WITH_FAST_DEL
from rmc.resources.rmc import (
    CURRENT_FREEZE,
    P_VALUE,
    grouped_single_no_break_ht_path,
    simul_search_round_bucket_path,
    simul_sections_split_by_len_path,
)
from rmc.utils.simultaneous_breaks import process_section_group

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("run_batches_dataproc")
logger.setLevel(logging.INFO)


def main(args):
    """Search for two simultaneous breaks in transcripts without evidence of a single significant break."""
    try:
        logger.warning("This step should be run on an autoscaling cluster!")
        hl.init(
            log=(
                f"/round{args.search_num}search_for_two_breaks_run_batches_dataproc.log"
            ),
            tmp_dir=TEMP_PATH_WITH_FAST_DEL,
        )
        hl.default_reference("GRCh38")

        chisq_threshold = hl.eval(hl.qchisqtail(P_VALUE, 2))
        if args.p_value:
            chisq_threshold = hl.eval(hl.qchisqtail(args.p_value, 2))

        # Adding dummy values to variables to avoid pre-commit errors
        sections_to_run = list()
        over_threshold = None
        if args.run_all_sections:
            sections_to_run = list(
                hl.eval(
                    hl.experimental.read_expression(
                        simul_sections_split_by_len_path(
                            is_over_threshold=True,
                            search_num=args.search_num,
                            freeze=args.freeze,
                        )
                    )
                )
            ) + list(
                hl.eval(
                    hl.experimental.read_expression(
                        simul_sections_split_by_len_path(
                            is_over_threshold=False,
                            search_num=args.search_num,
                            freeze=args.freeze,
                        )
                    )
                )
            )

        if args.run_sections_over_threshold:
            over_threshold = True
            sections_to_run = list(
                hl.eval(
                    hl.experimental.read_expression(
                        simul_sections_split_by_len_path(
                            search_num=args.search_num,
                            is_over_threshold=True,
                            freeze=args.freeze,
                        )
                    )
                )
            )

        if args.run_sections_under_threshold:
            over_threshold = False
            sections_to_run = list(
                hl.eval(
                    hl.experimental.read_expression(
                        simul_sections_split_by_len_path(
                            search_num=args.search_num,
                            is_over_threshold=False,
                            freeze=args.freeze,
                        )
                    )
                )
            )

        logger.info(
            "Found %i transcripts or transcript sections to search...",
            len(sections_to_run),
        )

        if args.group_size:
            logger.info(
                "Splitting transcripts/transcript sections into groups of %i",
                args.group_size,
            )
            section_groups = [
                sections_to_run[x : x + args.group_size]
                for x in range(0, len(sections_to_run), args.group_size)
            ]
        else:
            logger.info("Running transcripts/transcript sections one at a time...")
            section_groups = [[section] for section in sections_to_run]

        raw_path = simul_search_round_bucket_path(
            search_num=args.search_num,
            bucket_type="raw_results",
            freeze=args.freeze,
        )
        for counter, group in enumerate(section_groups):
            if over_threshold is None:
                output_ht_path = f"{raw_path}/simul_break_dataproc_all_{counter}.ht"
            elif over_threshold:
                output_ht_path = f"{raw_path}/simul_break_dataproc_{counter}.ht"
            else:
                output_ht_path = f"{raw_path}/simul_break_dataproc_under_{counter}.ht"

            if file_exists(output_ht_path):
                raise DataException(
                    f"Output already exists at {output_ht_path}! Double check before"
                    " running script again."
                )

            process_section_group(
                ht_path=grouped_single_no_break_ht_path(
                    args.search_num,
                    args.freeze,
                ),
                section_group=group,
                count=counter,
                over_threshold=True,
                output_ht_path=output_ht_path,
                output_n_partitions=args.output_n_partitions,
                chisq_threshold=chisq_threshold,
                split_list_len=args.split_list_len,
                read_if_exists=args.read_if_exists,
                save_chisq_ht=args.save_chisq_ht,
                freeze=args.freeze,
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(LOGGING_PATH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        This regional missense constraint script searches for two simultaneous breaks in transcripts without evidence
        of a single significant break.
        """,
        # Add default values for args to help message
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--freeze",
        help="RMC data freeze number",
        type=int,
        default=CURRENT_FREEZE,
    )
    parser.add_argument(
        "--p-value",
        help="""
        p-value significance threshold for single break search.
        Used to determine chi square threshold for likelihood ratio test.

        If not specified, script will default to threshold set
        in `P_VALUE`.
        """,
        type=float,
    )
    parser.add_argument(
        "--output-n-partitions",
        help="Number of desired partitions for output Tables. Default is 10.",
        type=int,
        default=10,
    )
    parser.add_argument(
        "--search-num",
        help=(
            "Search iteration number (e.g., second round of searching for two"
            " simultaneous breaks would be 2)."
        ),
        type=int,
    )
    parser.add_argument(
        "--group-size",
        help="""
        Number of transcripts/transcript sections to include in each group to be run.
        """,
        type=int,
    )
    parser.add_argument(
        "--split-list-len",
        help=(
            "Max length to divide transcript/sections observed or expected missense and"
            " position lists into."
        ),
        type=int,
        default=500,
    )
    section_ids = parser.add_mutually_exclusive_group()
    section_ids.add_argument(
        "--run-all-sections",
        help="Search for simultaneous breaks in all sections (regardless of length).",
        action="store_true",
    )
    section_ids.add_argument(
        "--run-sections-over-threshold",
        help="Search for simultaneous breaks in sections that are over length cutoff.",
        action="store_true",
    )
    section_ids.add_argument(
        "--run-sections-under-threshold",
        help="Search for simultaneous breaks in sections that are under length cutoff.",
        action="store_true",
    )
    parser.add_argument(
        "--read-if-exists",
        help="Use temporary Tables if they already exist.",
        action="store_true",
    )
    parser.add_argument(
        "--save-chisq-ht",
        help="""
        Save temporary Table that contains chi square significance values
        for all possible loci. Note that chi square values will be missing for
        any loci that would divide a transcript into subsections with fewer than
        `MIN_EXP_MIS` expected missense variants.

        NOTE that this temporary Table should only get saved once per gnomAD version,
        (save when running RMC pipeline at strictest p-value threshold).
        """,
        action="store_true",
    )

    args = parser.parse_args()

    main(args)
