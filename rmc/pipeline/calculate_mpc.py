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
from rmc.resources.rmc import CURRENT_FREEZE, context_with_oe, mpc_release
from rmc.slack_creds import slack_token
from rmc.utils.mpc import annotate_mpc, prepare_pop_path_ht, run_regressions

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
                use_model_formula=args.use_model_formula,
                model_formula=args.model_formula,
                variables=args.variables.split(","),
                additional_variables=args.extra_variables.split(","),
                overwrite_temp=args.overwrite_temp,
                do_k_fold_training=args.do_k_fold_training,
                freeze=args.freeze,
            )

        if args.command == "create-mpc-release":
            hl.init(log="/create_mpc_release.log", tmp_dir=temp_dir)
            annotate_mpc(
                ht=context_with_oe.versions[args.freeze].ht().select(),
                output_ht_path=mpc_release.versions[args.freeze].path,
                temp_label="_release",
                use_release=False,
                overwrite_temp=args.overwrite_temp,
                freeze=args.freeze,
            )

        if args.command == "annotate-hts":
            hl.init(log="/annotate_hts.log", tmp_dir=temp_dir)
            mpc_bucket_path = f"{MPC_PREFIX}/{CURRENT_GNOMAD_VERSION}/{args.freeze}"
            # TODO: Add support for annotating with models from specific folds (all folds?)

            if args.dd:
                from rmc.resources.reference_data import ndd_de_novo

                dd_ht = ndd_de_novo.ht()
                annotate_mpc(
                    ht=dd_ht,
                    output_ht_path=f"{mpc_bucket_path}/dd_mpc_annot.ht",
                    temp_label="_dd",
                    model_train_fold=args.model_train_fold,
                    use_release=args.use_release,
                    overwrite_temp=args.overwrite_temp,
                    freeze=args.freeze,
                )

            if args.specify_ht:
                annotate_mpc(
                    ht=hl.read_table(args.ht_in_path),
                    output_ht_path=args.ht_out_path,
                    overwrite_temp=args.overwrite_temp,
                    model_train_fold=args.model_train_fold,
                    use_release=args.use_release,
                    temp_label=args.temp_label,
                    freeze=args.freeze,
                )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(LOGGING_PATH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This regional missense constraint script calculates the MPC (missense badness,"
        " PolyPhen-2, and regional missense constraint) score."
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
        Generate (or prepare generation of) k-fold MPC models trained on training sets from respective folds.

        Otherwise, one model is generated, trained on all training transcripts.
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
        "--use-model-formula",
        help="Use specified model formula in regression.",
        action="store_true",
    )
    run_glm.add_argument(
        "--model-formula",
        help=(
            "R-style model formula to use in regression. Required if"
            " --use-model-formula is set."
        ),
        default="pop_v_path ~ oe + misbad + oe:misbad + polyphen + oe:polyphen",
    )
    run_glm.add_argument(
        "--variables",
        help=(
            "Comma separated string of variables to include in all logistic regression."
        ),
        default="oe,misbad,polyphen",
    )
    run_glm.add_argument(
        "--extra-variables",
        help=(
            "Comma separated string of additional variables to include in single"
            " variable regressions."
        ),
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
        "--model-train-fold",
        help="Fold number in training set used to generate MPC model.",
    )
    annotate_hts.add_argument(
        "--use-release",
        help="""
        Use scores in MPC release hail Table to annotate.

        Otherwise, scores will be directly computed using specified MPC model.
        """,
        action="store_true",
    )
    annotate_hts.add_argument(
        "--temp-label",
        help="""
        Suffix to add to temporary data paths to avoid conflicting names for different models.
        """,
    )
    annotate_hts.add_argument(
        "--dd",
        help=(
            "Calculate MPC for de novo variants from developmental disorder (DD) cases"
            " and controls."
        ),
        action="store_true",
    )
    annotate_hts.add_argument(
        "--specify-ht",
        help="Calculate MPC for variants in specified hail Table.",
        action="store_true",
    )
    annotate_hts.add_argument(
        "--ht-in-path",
        help=(
            "Path to input hail Table for MPC calculations. Required if --specify-ht is"
            " set."
        ),
    )
    annotate_hts.add_argument(
        "--ht-out-path",
        help=(
            "Output path for hail Table after adding MPC annotation. Required if"
            " --specify-ht is set."
        ),
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
