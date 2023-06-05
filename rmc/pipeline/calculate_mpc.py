"""
This script calculates the MPC score.

MPC (missense badness, PolyPhen-2, and regional missense constraint) is a composite score
that predicts the deleteriousness of any given missense variant.
"""
import argparse
import logging

import hail as hl
from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import LOGGING_PATH, MPC_PREFIX, TEMP_PATH_WITH_FAST_DEL
from rmc.resources.resource_utils import CURRENT_GNOMAD_VERSION
from rmc.resources.rmc import CURRENT_FREEZE
from rmc.slack_creds import slack_token
from rmc.utils.mpc import (
    annotate_mpc,
    create_mpc_release_ht,
    prepare_pop_path_ht,
    run_regressions,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("calculate_mpc")
logger.setLevel(logging.INFO)


def main(args):
    """Calculate MPC (Missense badness, Polyphen-2, and Constraint) score."""
    temp_dir = f"{TEMP_PATH_WITH_FAST_DEL}/mpc/"
    try:
        if args.command == "prepare-ht":
            hl.init(log="/write_pop_path_ht.log", tmp_dir=temp_dir)
            prepare_pop_path_ht(
                overwrite_temp=args.overwrite_temp,
                do_k_fold_training=args.do_k_fold_training,
                freeze=args.freeze,
            )

        if args.command == "run-glm":
            hl.init(log="/run_regressions_using_glm.log", tmp_dir=temp_dir)
            run_regressions(
                variables=args.variables.split(","),
                additional_variables=args.extra_variables.split(","),
                overwrite_temp=args.overwrite_temp,
                do_k_fold_training=args.do_k_fold_training,
                freeze=args.freeze,
            )

        if args.command == "create-mpc-release":
            hl.init(log="/create_mpc_release.log", tmp_dir=temp_dir)
            create_mpc_release_ht(
                overwrite_temp=args.overwrite_temp,
                freeze=args.freeze,
            )

        if args.command == "annotate-hts":
            hl.init(log="/annotate_hts.log", tmp_dir=temp_dir)
            mpc_bucket_path = f"{MPC_PREFIX}/{CURRENT_GNOMAD_VERSION}/{args.freeze}"
            if args.clinvar:
                from rmc.resources.reference_data import clinvar_plp_mis_haplo

                annotate_mpc(
                    ht=clinvar_plp_mis_haplo.ht(),
                    output_ht_path=f"{mpc_bucket_path}/clinvar_mpc_annot.ht",
                    overwrite_temp=args.overwrite_temp,
                    freeze=args.freeze,
                )

            if args.dd:
                from rmc.resources.reference_data import ndd_de_novo

                dd_ht = ndd_de_novo.ht()
                dd_ht = dd_ht.explode("sample_set")
                case_ht = dd_ht.filter(dd_ht.sample_set != "control")
                annotate_mpc(
                    ht=case_ht,
                    output_ht_path=f"{mpc_bucket_path}/dd_case_mpc_annot.ht",
                    overwrite_temp=args.overwrite_temp,
                    freeze=args.freeze,
                )
                control_ht = dd_ht.filter(dd_ht.sample_set == "control")
                shared_variants = case_ht.join(control_ht, how="inner")
                logger.info(
                    "%i variants overlapped between cases and controls",
                    shared_variants.count(),
                )
                annotate_mpc(
                    ht=control_ht,
                    output_ht_path=f"{mpc_bucket_path}/dd_control_mpc_annot.ht",
                    overwrite_temp=args.overwrite_temp,
                    freeze=args.freeze,
                )

            if args.gnomad_exomes:
                from gnomad.resources.grch37.gnomad import public_release

                ht = public_release("exomes").ht()
                annotate_mpc(
                    ht=ht,
                    output_ht_path=f"{mpc_bucket_path}/gnomAD_mpc_annot.ht",
                    overwrite_temp=args.overwrite_temp,
                    freeze=args.freeze,
                )

            if args.specify_ht:
                annotate_mpc(
                    ht=hl.read_table(args.ht_in_path),
                    output_ht_path=args.ht_out_path,
                    overwrite_temp=args.overwrite_temp,
                    freeze=args.freeze,
                )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(LOGGING_PATH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This regional missense constraint script calculates the MPC (missense badness, PolyPhen-2, and regional missense constraint) score."
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

    # Create subparsers for each step
    # Need to specify `dest` to be able to check which subparser is being invoked
    # `dest`: https://docs.python.org/3/library/argparse.html#dest
    subparsers = parser.add_subparsers(title="command", dest="command", required=True)

    prepare_ht = subparsers.add_parser(
        "prepare-ht",
        help="""
        Prepare Table with 'population' (common missense variants in gnomAD) and 'pathogenic'
        (ClinVar pathogenic/likely pathogenic missense variants in severe haploinsufficient genes) variants.

        This step joins gnomAD and ClinVar variants, annotates them with PolyPhen-2, missense badness,
        CADD (raw and phred), BLOSUM, Grantham, and missense observed/expected (OE) ratio, and removes
        any variants with undefined annotations.
        """,
    )

    run_glm = subparsers.add_parser(
        "run-glm",
        help="""
        Run logistic regressions on different models (single variable, joint).

        This step chooses a model based on the lowest AIC value and stores the
        model coefficients as a pickle.
        """,
    )
    run_glm.add_argument(
        "--variables",
        help="Comma separated string of variables to include in all logistic regression.",
        default="oe,misbad,polyphen",
    )
    run_glm.add_argument(
        "--extra-variables",
        help="Comma separated string of additional variables to include in single variable regressions.",
        default="blosum,grantham",
    )

    create_mpc_release = subparsers.add_parser(
        "create-mpc-release",
        help="""
        Create MPC release Table (VEP context Table filtered to missense variants in canonical, non-outlier transcripts).
        """,
    )

    annotate_hts = subparsers.add_parser(
        "annotate-hts", help="Annotate specified dataset with MPC."
    )
    annotate_hts.add_argument(
        "--clinvar", help="Calculate MPC for ClinVar variants", action="store_true"
    )
    annotate_hts.add_argument(
        "--dd",
        help="Calculate MPC for de novo variants from developmental disorder (DD) cases and controls",
        action="store_true",
    )
    annotate_hts.add_argument(
        "--gnomad-exomes",
        help="Calculate MPC for all gnomAD exomes variants",
        action="store_true",
    )
    annotate_hts.add_argument(
        "--specify-ht",
        help="Calculate MPC for variants in specified hail Table",
        action="store_true",
    )
    annotate_hts.add_argument(
        "--ht-in-path",
        help="Path to input hail Table for MPC calculations. Required if --specify-ht is set.",
    )
    annotate_hts.add_argument(
        "--ht-out-path",
        help="Output path for hail Table after adding MPC annotation. Required if --specify-ht is set.",
    )
    # TODO: Review if all these args are needed

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
