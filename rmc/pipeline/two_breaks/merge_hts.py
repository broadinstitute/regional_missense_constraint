"""
This script merges all intermediate simultaneous breaks result Tables into a single Table.

This step should be run in Dataproc.
"""
import argparse
import logging
import subprocess

import hail as hl

from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import (
    LOGGING_PATH,
    no_breaks,
    not_one_break,
    simul_break,
    simul_break_temp,
)
from rmc.slack_creds import slack_token


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("merge_hts")
logger.setLevel(logging.INFO)


def main(args):
    """Merge all simultaneous breaks intermediate results into single Table."""
    try:
        hl.init(log="/search_for_two_breaks_merge_hts.log")

        logger.info("Collecting all HT paths...")
        intermediate_hts = []
        ht_bucket = f"{simul_break_temp}/hts/"
        temp_ht_paths = (
            subprocess.check_output(["gsutil", "ls", ht_bucket])
            .decode("utf8")
            .strip()
            .split("\n")
        )
        for ht_path in temp_ht_paths:
            ht_path = ht_path.strip("/")
            temp = hl.read_table(ht_path)
            intermediate_hts.append(temp)

        ht = intermediate_hts[0].union(*intermediate_hts[1:])
        ht = ht.checkpoint(simul_break.path, overwrite=args.overwrite)
        logger.info("Wrote simultaneous breaks HT with %i lines", ht.count())

        # Collect all transcripts with two simultaneous breaks
        simul_break_transcripts = ht.aggregate(hl.agg.collect_as_set(ht.transcript),)
        logger.info(
            "%i transcripts had two simultaneous breaks", len(simul_break_transcripts),
        )
        simul_break_transcripts = hl.literal(simul_break_transcripts)

        logger.info(
            "Getting transcripts with no evidence of regional missense constraint..."
        )
        context_ht = not_one_break.ht()
        context_ht = context_ht.filter(
            ~simul_break_transcripts.contains(context_ht.transcript)
        )
        context_ht.write(no_breaks.path, overwrite=args.overwrite)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(LOGGING_PATH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This regional missense constraint script merges all intermediate simultaneous breaks results Tables into a single Table."
    )
    parser.add_argument(
        "--overwrite", help="Overwrite existing data.", action="store_true"
    )
    parser.add_argument(
        "--slack-channel",
        help="Send message to Slack channel/user.",
        default="@kc (she/her)",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
