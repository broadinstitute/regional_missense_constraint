"""
This script prepares data or runs analyses to assess the MPC score.

MPC (missense badness, PolyPhen-2, and regional missense constraint) is a composite score
that predicts the deleteriousness of any given missense variant.
"""
import argparse
import logging

import hail as hl

from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import LOGGING_PATH, MPC_PREFIX
from rmc.resources.resource_utils import CURRENT_VERSION
from rmc.slack_creds import slack_token
from rmc.utils.mpc import (
    prep_mpc_comparison_ht,
    prep_mpc_histogram_tsv,
    prep_rate_ratio_tsv,
)


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("calculate_mpc")
logger.setLevel(logging.INFO)


def main(args):
    """Assess MPC (Missense badness, Polyphen-2, and Constraint) score."""
    try:
        if args.command == "prepare-histogram-tsv":
            hl.init(log="/prep_mpc_histogram.log")
            prep_mpc_histogram_tsv(
                case_ht=hl.read_table(args.case_path),
                control_ht=hl.read_table(args.control_path),
                output_tsv_path=args.output_tsv_path,
                keep_asd=args.asd,
                keep_dd=args.dd,
            )

        if args.command == "prepare-rate-ratio-tsv":
            hl.init(log="/prep_rate_ratio.log")
            prep_rate_ratio_tsv(
                output_tsv_path=args.output_tsv_path,
                keep_asd=args.asd,
                keep_dd=args.dd,
            )

        if args.command == "prepare-mpc-comparison-table":
            hl.init(log="/prepare_mpc_comparison_ht.log")
            prep_mpc_comparison_ht(
                case_ht=hl.read_table(args.case_path),
                control_ht=hl.read_table(args.control_path),
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(LOGGING_PATH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script assesses the MPC (missense badness, PolyPhen-2, and regional missense constraint) score."
    )
    parser.add_argument(
        "--overwrite", help="Overwrite existing data.", action="store_true"
    )
    parser.add_argument(
        "--slack-channel",
        help="Send message to Slack channel/user.",
    )
    parser.add_argument(
        "--asd",
        help="Specify this flag to keep variants from cases with Autism Spectrum Disorder (ASD).",
        action="store_true",
    )
    parser.add_argument(
        "--dd",
        help="Specify this flag to keep variants from cases with developmental disorders (DD).",
        action="store_true",
    )
    parser.add_argument(
        "--output-tsv-path",
        help="Output path for TSV (either for prepare-histogram-tsv or prepare-rate-ratio-tsv step).",
    )
    parser.add_argument(
        "--case-path",
        help="Path to Hail Table of de novo variants from NDD cases annotated with MPC.",
        default=f"{MPC_PREFIX}/{CURRENT_VERSION}/dd_case_mpc_annot.ht",
    )
    parser.add_argument(
        "--control-path",
        help="Path to Hail Table of de novo variants from NDD controls annotated with MPC.",
        default=f"{MPC_PREFIX}/{CURRENT_VERSION}/dd_control_mpc_annot.ht",
    )

    # Create subparsers for each step
    # Need to specify `dest` to be able to check which subparser is being invoked
    # `dest`: https://docs.python.org/3/library/argparse.html#dest
    subparsers = parser.add_subparsers(title="command", dest="command", required=True)

    prepare_hist_tsv = subparsers.add_parser(
        "prepare-histogram-tsv",
        help="""
        Prepare TSV using de novo variants (annotated with MPC) from neurodevelopmental disorders (NDD) cases or controls.

        TSV contains only two columns:
            - Case vs control status of variant
            - Variant's MPC score

        This step also creates case control Hail Table used downstream to calculate rate ratios of MPC in de novo variants
        from cases vs controls.
        """,
    )

    prepare_rate_tsv = subparsers.add_parser(
        "prepare-rate-ratio-tsv",
        help="""
        Prepare TSV of MPC bins and corresponding variant counts (
            total variants from cases, total variants from controls, total number of cases, total number of controls,
            rate per case, rate per control
        ).

        This TSV is used as input in a two-sided Poisson exact test to calculate MPC rate ratios in de novo variants from
        cases and controls.
        """,
    )

    prepare_mpc_comparison_ht = subparsers.add_parser(
        "prepare-mpc-comparison-table",
        help="""
        Prepare Hail Table of using de novo variants from neurodevelopmental disorders (NDD) cases or controls
        annotated with MPC, Polyphen-2, SIFT, CADD, and REVEL scores.

        Used to compare MPC performance against these scores.
        """,
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
