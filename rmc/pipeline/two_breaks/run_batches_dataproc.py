"""
This script searches for two simultaneous breaks in groups of transcripts using Hail Query within Google Cloud Dataproc.

This script should be run only on transcripts that are greater than or equal
to the --transcript-len-threshold specified in `prepare_transcripts.py` if these transcripts are too slow or getting preempted in Hail Batch.
Transcripts smaller than --transcript-len-threshold  should be run using `run_batches.py` as they run quickly and inexpensively in Hail Batch.

If using this step to run TTN, use a large autoscaling cluster (highmem-8, scales to 100 preemptibles).
Otherwise, an autoscaling cluster of highmem-8s that scales to 50 preemptibles should suffice.
"""
import argparse
import logging

import hail as hl

from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import (
    LOGGING_PATH,
    not_one_break_grouped,
    simul_break_temp,
)
from rmc.slack_creds import slack_token
from rmc.utils.simultaneous_breaks import process_transcript_group

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
        hl.init(log="/search_for_two_breaks_run_batches_dataproc.log")
        transcripts_to_run = args.transcripts_to_run.split(",")
        if args.group_size:
            logger.info("Splitting transcripts into groups of %i", args.group_size)
            transcript_groups = [
                transcripts_to_run[x : x + args.group_size]
                for x in range(0, len(transcripts_to_run), args.group_size)
            ]
        else:
            logger.info("Running transcripts one at a time...")
            transcript_groups = [[transcript] for transcript in transcripts_to_run]

        for counter, group in enumerate(transcript_groups):
            process_transcript_group(
                ht_path=not_one_break_grouped.path,
                transcript_group=group,
                over_threshold=True,
                output_ht_path=f"{simul_break_temp}/hts/simul_break_dataproc_{counter}.ht",
                output_tsv_path=f"{simul_break_temp}/success_files",
                temp_ht_path=f"{simul_break_temp}",
                chisq_threshold=args.chisq_threshold,
                split_window_size=args.window_size,
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
        "--chisq-threshold",
        help="Chi-square significance threshold. Value should be 9.2 (value adjusted from ExAC code due to discussion with Mark).",
        type=float,
        default=9.2,
    )
    parser.add_argument(
        "--overwrite", help="Overwrite existing data.", action="store_true"
    )
    parser.add_argument(
        "--slack-channel",
        help="Send message to Slack channel/user.",
        default="@kc (she/her)",
    )

    parser.add_argument(
        "--group-size",
        help="""
        Number of transcripts to include in each group of transcripts to be run.
        Default is 50.
        """,
        type=int,
        default=50,
    )
    parser.add_argument(
        "--window-size",
        help="Size of windows to split transcripts. Default is 500.",
        type=int,
        default=500,
    )
    parser.add_argument(
        "--transcripts-to-run", help="Comma separated list of transcript IDs to run."
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
