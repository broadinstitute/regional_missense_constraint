import argparse
import logging

import hail as hl

from gnomad.utils.slack import slack_notifications
from gnomad_lof.constraint_utils.constraint_basics import prepare_ht
from rmc.resources.basics import LOGGING_PATH, temp_path
from rmc.resources.grch37.exac import filtered_exac
from rmc.resources.grch37.gnomad import (
    constraint_ht,
    filtered_exomes,
    processed_exomes,
    prop_obs_coverage,
)
from rmc.resources.grch37.reference_data import processed_context
from rmc.resources.resource_utils import MISSENSE
from rmc.slack_creds import slack_token
from rmc.utils.constraint import (
    calculate_exp_per_base,
    calculate_exp_per_transcript,
    calculate_observed,
    GROUPINGS,
)
from rmc.utils.generic import (
    filter_to_region_type,
    generate_models,
    get_coverage_correction_expr,
    keep_criteria,
    process_context_ht,
    process_vep,
)


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_missense_constraint")
logger.setLevel(logging.INFO)


def main(args):

    hl.init(log="/RMC.log")
    exac = args.exac

    # Add transcript to core grouping fields
    groupings = GROUPINGS.append("transcript")

    try:
        if args.pre_process_data:
            logger.warning("Code currently only processes b37 data!")
            logger.info(
                "Filtering gnomAD exomes HT to SNPs and annotating with variant type, methylation, and coverage..."
            )
            exome_ht = processed_exomes.ht()
            exome_ht = prepare_ht(exome_ht, args.trimers)

            logger.info(
                "Filtering gnomAD exomes HT to missense variants in canonical transcripts only..."
            )
            exome_ht = process_vep(exome_ht, filter_csq=True, csq=MISSENSE)
            exome_ht = exome_ht.filter(keep_criteria(exome_ht))
            exome_ht.write(filtered_exomes.path, overwrite=args.overwrite)

            logger.info("Preprocessing reference fasta (context) HT...")
            context_ht = process_context_ht("GRCh37", args.trimers)

            logger.info(
                "Filtering context HT to all sites not found in gnomAD exomes + all rare, covered sites in gnomAD"
            )
            exome_join = exome_ht[context_ht.key]
            context_ht = context_ht.filter(
                hl.is_missing(exome_join) | keep_criteria(exome_join)
            )

            # NOTE: should use ~30k-40k partitions here
            context_ht = context_ht.repartition(args.n_partitions)
            context_ht.write(processed_context.path, overwrite=args.overwrite)

            logger.info("Done preprocessing files")

        if args.prep_for_constraint:
            logger.info("Reading in exome HT...")
            if exac:
                exome_ht = filtered_exac.ht()

            else:
                exome_ht = filtered_exomes.ht()

            logger.info("Reading in context HT...")
            # TODO: context HT wrote out with only ~104 partitions? need to repartition
            context_ht = processed_context.ht()

            logger.info("Building plateau and coverage models...")
            coverage_ht = prop_obs_coverage.ht()
            coverage_x_ht = hl.read_table(
                prop_obs_coverage.path.replace(".ht", "_x.ht")
            )
            coverage_y_ht = hl.read_table(
                prop_obs_coverage.path.replace(".ht", "_y.ht")
            )

            (
                coverage_model,
                plateau_models,
                plateau_x_models,
                plateau_y_models,
            ) = generate_models(
                coverage_ht, coverage_x_ht, coverage_y_ht, trimers=args.trimers
            )

            context_ht = context_ht.annotate_globals(
                plateau_models=plateau_models,
                plateau_x_models=plateau_x_models,
                plateau_y_models=plateau_y_models,
                coverage_model=coverage_model,
            )
            context_ht = context_ht.annotate(
                coverage_correction=get_coverage_correction_expr(
                    context_ht.exome_coverage, coverage_model, args.high_cov_cutoff,
                )
            )

            if not args.skip_calc_oe:
                logger.info(
                    "Creating autosomes-only, chrX non-PAR-only, and chrY non-PAR-only HT versions..."
                )
                context_x_ht = filter_to_region_type(context_ht, "chrX")
                context_y_ht = filter_to_region_type(context_ht, "chrY")
                context_auto_ht = filter_to_region_type(context_ht, "autosomes")

                logger.info("Calculating expected values per transcript...")
                exp_ht = calculate_exp_per_transcript(
                    context_auto_ht, locus_type="autosomes", groupings=groupings
                )
                exp_x_ht = calculate_exp_per_transcript(
                    context_x_ht, locus_type="X", groupings=groupings
                )
                exp_y_ht = calculate_exp_per_transcript(
                    context_y_ht, locus_type="Y", groupings=groupings
                )
                exp_ht = exp_ht.union(exp_x_ht).union(exp_y_ht)

                logger.info(
                    "Aggregating total observed variant counts per transcript..."
                )
                obs_ht = calculate_observed(exome_ht, exac)

                logger.info(
                    "Annotating total observed and expected values and overall observed/expected value "
                    "(capped at 1) per transcript..."
                )
                context_ht = context_ht.annotate(
                    total_exp=exp_ht[context_ht.transcript].expected,
                    total_obs=obs_ht[context_ht.transcript].observed,
                )
                context_ht = context_ht.annotate(
                    overall_obs_exp=hl.min(
                        context_ht.total_obs / context_ht.total_exp, 1
                    )
                )

            else:
                logger.warning(
                    "Using observed and expected values calculated on gnomAD v2.1.1 exomes..."
                )
                gnomad_constraint_ht = (
                    constraint_ht.ht()
                    .key_by("transcript")
                    .select("obs_mis", "exp_mis", "oe_mis")
                )

                # Filter to canonical transcripts only and rename fields
                gnomad_constraint_ht = gnomad_constraint_ht.filter(
                    gnomad_constraint_ht.canonical
                )
                gnomad_constraint_ht = gnomad_constraint_ht.transmute(
                    total_obs=gnomad_constraint_ht.obs_mis,
                    total_exp=gnomad_constraint_ht.exp_mis,
                    overall_obs_exp=gnomad_constraint_ht.oe_mis,
                )
                context_ht = context_ht.annotate(
                    **gnomad_constraint_ht[context_ht.transcript]
                )

            logger.info(
                "Annotating context HT with number of observed and expected variants per site..."
            )
            context_ht = context_ht.annotate(_obs=exome_ht.index(context_ht.key))
            context_ht = context_ht.transmute(
                observed=hl.int(hl.is_defined(context_ht._obs))
            )

            context_ht = calculate_exp_per_base(context_ht, groupings)
            context_ht = context_ht.write(f"{temp_path}/context_obs_exp_annot.ht")

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(LOGGING_PATH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script searches for regional missense constraint in gnomAD"
    )
    parser.add_argument(
        "--trimers", help="Use trimers instead of heptamers", action="store_true"
    )
    parser.add_argument(
        "--exac", help="Use ExAC Table (not gnomAD Table)", action="store_true"
    )
    parser.add_argument(
        "n_partitions", help="Desired number of partitions for output data", type=int,
    )
    parser.add_argument(
        "--high_cov_cutoff",
        help="Coverage threshold for a site to be considered high coverage",
        type=int,
        default=40,
    )
    parser.add_argument(
        "--pre_process_data", help="Pre-process data", action="store_true"
    )
    parser.add_argument(
        "--prep_for_constraint",
        help="Prepare tables for constraint calculations",
        action="store_true",
    )
    parser.add_argument(
        "--skip_calc_oe",
        help="Skip observed and expected variant calculations per transcript. Relevant only to gnomAD v2.1.1!",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite", help="Overwrite existing data", action="store_true"
    )
    parser.add_argument(
        "--slack_channel", help="Send message to Slack channel/user", default="@kc"
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
