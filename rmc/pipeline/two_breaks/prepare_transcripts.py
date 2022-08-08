"""
This script prepares the inputs to the two simultaneous breaks search.

This script has two possible steps:
- Create grouped version of not one break Table
- Split transcripts to search by number of possible missense variants per transcript.

Both steps should be run in Dataproc.
"""
import argparse
import logging

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import LOGGING_PATH
from rmc.resources.grch37.reference_data import gene_model
from rmc.resources.grch37.rmc import not_one_break, not_one_break_grouped
from rmc.slack_creds import slack_token
from rmc.utils.simultaneous_breaks import (
    group_not_one_break_ht,
    split_transcripts_by_len,
)


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("prepare_transcripts")
logger.setLevel(logging.INFO)


def main(args):
    """Prepare input Table and transcripts for two simultaneous breaks search."""
    try:
        if args.command == "create-grouped-ht":
            hl.init(log="/search_for_two_breaks_create_grouped_ht.log")

            if not file_exists(not_one_break_grouped.path) or args.create_grouped_ht:
                if args.min_num_obs == 0:
                    # Make sure user didn't specify a min obs of 0
                    raise DataException(
                        "Minimum number of observed variants must be greater than zero!"
                    )

                logger.info(
                    "Creating grouped HT with lists of cumulative observed and expected missense values..."
                )
                group_not_one_break_ht(
                    ht=not_one_break.ht(),
                    transcript_ht=gene_model.ht(),
                    get_min_window_size=args.get_min_window_size,
                    get_total_exome_bases=args.get_total_exome_bases,
                    get_total_gnomad_missense=args.get_total_gnomad_missense,
                    min_window_size=args.min_window_size,
                    min_num_obs=args.min_num_obs,
                )

        if args.command == "split-transcripts":
            logger.warning("This step should be run in Dataproc!")
            hl.init(log="/search_for_two_breaks_split_transcripts.log")
            split_transcripts_by_len(
                ht=not_one_break_grouped.ht(),
                transcript_len_threshold=args.transcript_len_threshold,
                ttn_id=args.ttn,
                overwrite=args.overwrite,
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(LOGGING_PATH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This regional missense constraint script prepares the input Table and transcripts for the two simultaneous breaks search.",
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

    # Create subparsers for each step
    # Need to specify `dest` to be able to check which subparser is being invoked
    # `dest`: https://docs.python.org/3/library/argparse.html#dest
    subparsers = parser.add_subparsers(title="command", dest="command", required=True)

    create_grouped_ht = subparsers.add_parser(
        "create-grouped-ht",
        help="""
        Create hail Table grouped by transcript with cumulative observed and expected missense values collected into lists.
        This step should be run in Dataproc.
        """,
    )
    min_window_group = create_grouped_ht.add_mutually_exclusive_group()
    min_window_group.add_argument(
        "--min-window-size",
        help="Smallest possible window size for simultaneous breaks. Determined by running --get-min-window-size.",
        type=int,
    )
    min_window_group.add_argument(
        "--get-min-window-size",
        help="Determine smallest possible window size for simultaneous breaks.",
        action="store_true",
    )
    create_grouped_ht.add_argument(
        "--get-total-exome-bases",
        help="Get total number of bases in the exome. If not set, will pull default value from TOTAL_EXOME_BASES.",
        action="store_true",
    )
    create_grouped_ht.add_argument(
        "--get-total-gnomad-missense",
        help="Get total number of missense variants in gnomAD. If not set, will pull default value from TOTAL_GNOMAD_MISSENSE.",
        action="store_true",
    )
    create_grouped_ht.add_argument(
        "--min-num-obs",
        help="Number of observed variants. Used when determining the smallest possible window size for simultaneous breaks.",
        type=int,
        default=10,
    )

    split_transcripts = subparsers.add_parser(
        "split-transcripts",
        help="""
        Split transcripts based on number of possible missense positions.
        This is used to create batches of transcripts to run through search for two breaks code.
        This step should be run in Dataproc.
        """,
    )
    split_transcripts.add_argument(
        "--transcript-len-threshold",
        help="Cutoff for number of possible missense positions in transcript. Used to create batches of transcripts.",
        type=int,
        default=5000,
    )
    split_transcripts.add_argument(
        "--ttn",
        help="TTN transcript ID. TTN is so large that it needs to be treated separately.",
        default="ENST00000589042",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
