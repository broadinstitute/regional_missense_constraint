"""
This script calculates missense badness.

Missense badness is a score that reflects the deleteriousness of each specific amino acid substitution.
"""
import argparse
import logging

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import LOGGING_PATH
from rmc.slack_creds import slack_token
from rmc.utils.missense_badness import prepare_amino_acid_ht


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("calculate_missense_badness")
logger.setLevel(logging.INFO)


def main(args):
    """Calculate missense badness."""
    if not args.command:
        raise DataException("Please specify command for this script!")

    try:

        if args.command == "prepare-amino-acid-ht":
            hl.init(log="/calc_misbad_prep_context_gamma_ht.log")
            prepare_amino_acid_ht()

        if args.command == "create-misbad":
            hl.init(log="/calc_misbad_create_score.log")

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(LOGGING_PATH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This regional missense constraint script two simultaneous breaks in transcripts without evidence of a single significant break."
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
    subparsers = parser.add_subparsers(title="command", dest="command")

    prepare_amino_acid_ht = subparsers.add_parser(
        "prepare-amino-acid-ht",
        help="""
        Prepare Table with all possible amino acid substitutions and their missense observed to expected (OE) ratio.
        """,
        action="store_true",
    )

    create_misbad = subparsers.add_parser(
        "create-misbad",
        help="""
        Calculate missense badness score for each possible amino acid substitution.
        """,
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
