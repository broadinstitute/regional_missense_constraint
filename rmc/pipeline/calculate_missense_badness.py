"""
This script calculates missense badness.

Missense badness is a score that reflects the deleteriousness of each specific amino acid substitution.
"""
import argparse
import logging

import hail as hl
from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import LOGGING_PATH, TEMP_PATH_WITH_FAST_DEL
from rmc.resources.rmc import CURRENT_FREEZE
from rmc.slack_creds import slack_token
from rmc.utils.missense_badness import calculate_misbad, prepare_amino_acid_ht

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("calculate_missense_badness")
logger.setLevel(logging.INFO)


def main(args):
    """Calculate missense badness."""
    temp_dir = f"{TEMP_PATH_WITH_FAST_DEL}/mb/"
    try:
        if args.command == "prepare-ht":
            hl.init(log="/calc_misbad_prep_context_gamma_ht.log", tmp_dir=temp_dir)
            prepare_amino_acid_ht(
                overwrite_temp=args.overwrite_temp,
                do_k_fold_training=args.do_k_fold_training,
                freeze=args.freeze,
            )

        if args.command == "create-misbad":
            hl.init(log="/calc_misbad_create_score.log", tmp_dir=temp_dir)
            calculate_misbad(
                use_exac_oe_cutoffs=args.use_exac_oe_cutoffs,
                overwrite_temp=args.overwrite_temp,
                do_k_fold_training=args.do_k_fold_training,
                freeze=args.freeze,
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(LOGGING_PATH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This regional missense constraint script calculates missense badness."
    )
    parser.add_argument(
        "--overwrite-temp",
        help="Overwrite existing intermediate temporary data.",
        action="store_true",
    )
    parser.add_argument(
        "--freeze",
        help="RMC data freeze number",
        type=int,
        default=CURRENT_FREEZE,
    )
    parser.add_argument(
        "--slack-channel",
        help="Send message to Slack channel/user.",
    )
    parser.add_argument(
        "--do-k-fold-training",
        help="""
        Generate k-fold missense badness models trained on training sets from respective folds.
        Otherwise, one model is generated, trained on all training transcripts.
        """,
        action="store_true",
    )
    parser.add_argument(
        "--use-exac-oe-cutoffs",
        help="""
        Use the same missense OE cutoffs as in ExAC missense badness calculation.
        This removes rows with 0.6 < missense OE <= 0.8.
        Only relevant for `create-misbad` step.
        """,
        action="store_true",
    )

    # Create subparsers for each step
    # Need to specify `dest` to be able to check which subparser is being invoked
    # `dest`: https://docs.python.org/3/library/argparse.html#dest
    subparsers = parser.add_subparsers(title="command", dest="command", required=True)

    prepare_ht = subparsers.add_parser(
        "prepare-ht",
        help="""
        Prepare Table with all possible amino acid substitutions and their missense observed to expected (OE) ratio.
        """,
    )

    create_misbad = subparsers.add_parser(
        "create-misbad",
        help="""
        Calculate missense badness score for each possible amino acid substitution.
        """,
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
