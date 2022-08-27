"""
This script merges all intermediate simultaneous breaks result Tables into a single Table.

This step should be run in Dataproc.
"""
import argparse
import logging
import subprocess

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import LOGGING_PATH, SIMUL_BREAK_TEMP_PATH
from rmc.resources.grch37.rmc import (
    no_breaks,
    not_one_break,
    simul_break,
)
from rmc.slack_creds import slack_token


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("merge_hts")
logger.setLevel(logging.INFO)


ANNOTATIONS = {"max_chisq", "start_pos", "end_pos"}
"""
Set of annotations added during two simultaneous breaks search.

`max_chisq`: Chi square value associated with two breaks.
`start_pos`: Start position of two break window.
`end_pos`: End position of two break window.
"""


def main(args):
    """Merge all simultaneous breaks intermediate results into single Table."""
    try:
        hl.init(log="/search_for_two_breaks_merge_hts.log")

        logger.info("Collecting all HT paths...")
        intermediate_hts = []
        ht_bucket = f"{SIMUL_BREAK_TEMP_PATH}/hts/"
        temp_ht_paths = (
            subprocess.check_output(["gsutil", "ls", ht_bucket])
            .decode("utf8")
            .strip()
            .split("\n")
        )
        ht_count = 0
        for ht_path in temp_ht_paths:
            ht_path = ht_path.strip("/")
            if ht_path.endswith("ht"):
                ht_count += 1
                logger.info("Working on %s", ht_path)
                temp = hl.read_table(ht_path)
                if temp.count() > 0:
                    # Tables containing transcripts that are over the transcript length threshold are keyed by transcript, i, j
                    # Tables containing transcripts that are under the length threshold are keyed only by transcript
                    # Rekey all tables here and select only the required fields to ensure the union on line 83 is able to work
                    # A normal `.key_by` should work here, since transcripts are already part of the key fields
                    # (see https://github.com/hail-is/hail/blob/master/hail/src/main/scala/is/hail/expr/ir/TableIR.scala#L812)
                    # However, using `.key_by_assert_sorted` to explicitly avoid shuffling on this rekey
                    temp = temp._key_by_assert_sorted("transcript")
                    row_fields = set(temp.row)
                    if len(ANNOTATIONS.intersection(row_fields)) < 3:
                        raise DataException(
                            f"The following fields are missing from the temp table: {ANNOTATIONS.difference(row_fields)}!"
                        )
                    temp = temp.select("max_chisq", "start_pos", "end_pos")
                    intermediate_hts.append(temp)
                else:
                    logger.warning("%s had 0 rows", ht_path)
        logger.info("Found %i HTs and appended %i", ht_count, len(intermediate_hts))

        if len(intermediate_hts) == 0:
            raise DataException(
                "All temp tables had 0 rows. Please double check the temp tables!"
            )
        ht = intermediate_hts[0].union(*intermediate_hts[1:])
        ht = ht.checkpoint(simul_break.path, overwrite=args.overwrite)
        logger.info("Wrote simultaneous breaks HT with %i lines", ht.count())

        # Collect all transcripts with two simultaneous breaks
        simul_break_transcripts = ht.aggregate(hl.agg.collect_as_set(ht.transcript))
        logger.info(
            "%i transcripts had two simultaneous breaks",
            len(simul_break_transcripts),
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
        description="This regional missense constraint script merges all intermediate simultaneous breaks results Tables into a single Table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
