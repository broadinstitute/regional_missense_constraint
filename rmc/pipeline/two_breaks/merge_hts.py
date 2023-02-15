"""
This script merges all intermediate simultaneous breaks result Tables into a single Table.

This step should be run in Dataproc.
"""
import argparse
import logging

import hail as hl

from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import LOGGING_PATH, TEMP_PATH_WITH_FAST_DEL
from rmc.resources.rmc import simul_search_round_bucket_path
from rmc.utils.settings import SLACK_TOKEN
from rmc.utils.constraint import merge_simul_break_temp_hts


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("merge_hts")
logger.setLevel(logging.INFO)


def main(args):
    """Merge all simultaneous breaks intermediate results into single Table."""
    try:
        hl.init(
            log=f"/round{args.search_num}_search_for_two_breaks_merge_hts.log",
            tmp_dir=TEMP_PATH_WITH_FAST_DEL,
        )

        logger.info("Merging all temp HTs...")
        raw_path = simul_search_round_bucket_path(
            search_num=args.search_num,
            bucket_type="raw_results",
        )
        results_path = simul_search_round_bucket_path(
            search_num=args.search_num,
            bucket_type="final_results",
        )
        merge_simul_break_temp_hts(
            input_hts_path=raw_path,
            batch_phrase="under",
            query_phrase="dataproc",
            output_ht_path=f"{results_path}/merged.ht",
            overwrite=args.overwrite,
            google_project=args.google_project,
        )

        ht = hl.read_table(f"{results_path}/merged.ht")
        logger.info("Wrote temp simultaneous breaks HT with %i lines", ht.count())

        # Collect all transcripts and sections with two simultaneous breaks
        ht = ht.annotate(transcript=ht.section.split("_")[0])
        simul_break_transcripts = ht.aggregate(hl.agg.collect_as_set(ht.transcript))
        logger.info(
            "%i transcripts had two simultaneous breaks in this round",
            len(simul_break_transcripts),
        )
        simul_break_sections = ht.aggregate(hl.agg.collect_as_set(ht.section))
        hl.experimental.write_expression(
            simul_break_sections,
            f"{results_path}/sections.he",
            overwrite=args.overwrite,
        )
        logger.info(
            "%i transcript sections had two simultaneous breaks in this round",
            len(simul_break_sections),
        )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(LOGGING_PATH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This regional missense constraint script merges all intermediate simultaneous breaks results Tables into a single Table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--overwrite", help="Overwrite existing data.", action="store_true"
    )
    parser.add_argument(
        "--slack-channel",
        help="Send message to Slack channel/user.",
    )
    parser.add_argument(
        "--search-num",
        help="Search iteration number (e.g., second round of searching for two simultaneous breaks would be 2).",
        type=int,
    )
    parser.add_argument(
        "--google-project",
        help="""
            Google cloud project used to read from requester-pays buckets.
            """,
        default="broad-mpg-gnomad",
    )

    args = parser.parse_args()

    if args.slack_channel and SLACK_TOKEN:
        with slack_notifications(SLACK_TOKEN, args.slack_channel):
            main(args)
    else:
        main(args)
