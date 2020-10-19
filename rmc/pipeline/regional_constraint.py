import argparse
import logging

import hail as hl

from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.slack import slack_notifications
from gnomad_lof.constraint_utils.constraint_basics import build_models, prepare_ht
from rmc.resources.basics import LOGGING_PATH
from rmc.resources.grch37.exac import filtered_exac
from rmc.resources.grch37.gnomad import (
    constraint_ht,
    filtered_exomes,
    processed_exomes,
    prop_obs_coverage,
)
from rmc.resources.grch37.reference_data import processed_context
from rmc.slack_creds import slack_token
from rmc.utils.constraint import (
    calculate_exp_per_base,
    calculate_observed,
)
from rmc.utils.generic import (
    filter_to_region_type,
    filter_to_missense,
    get_coverage_correction_expr,
    process_context_ht,
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

    try:
        if args.pre_process_data:
            logger.info("Preprocessing reference fasta and gencode files...")
            process_context_ht("GRCh37", args.trimers)

            logger.info(
                "Filtering gnomAD exomes HT to SNPs and annotating with variant type, methylation, and coverage..."
            )
            exome_ht = processed_exomes.ht()
            exome_ht = prepare_ht(exome_ht, args.trimers)

            logger.info(
                "Filtering gnomAD exomes HT to missense variants in canonical transcripts only..."
            )
            exome_ht = filter_to_missense(exome_ht)
            exome_ht.write(filtered_exomes.path, overwrite=args.overwrite)

            logger.info("Done preprocessing files")

        logger.info("Reading in exome HT...")
        if exac:
            exome_ht = filtered_exac.ht()

        else:
            exome_ht = filtered_exomes.ht()

        logger.info("Reading in context HT...")
        # TODO: context HT wrote out with only ~104 partitions? need to repartition
        context_ht = processed_context.ht()

        if args.test:
            logger.info("Inferring build of exome HT...")
            rg = get_reference_genome(exome_ht.locus)

            logger.info("Filtering to chr20 for testing...")
            contigs = rg.contigs[19]
            context_ht = hl.filter_intervals(
                context_ht, [hl.parse_locus_interval(contigs, reference_genome=rg)]
            )
            exome_ht = hl.filter_intervals(
                exome_ht, [hl.parse_locus_interval(contigs, reference_genome=rg)]
            )

        logger.info("Building plateau and coverage models...")
        coverage_ht = prop_obs_coverage.ht()
        coverage_x_ht = hl.read_table(prop_obs_coverage.path.replace(".ht", "_x.ht"))
        coverage_y_ht = hl.read_table(prop_obs_coverage.path.replace(".ht", "_y.ht"))
        coverage_model, plateau_models = build_models(
            coverage_ht, args.trimers, weighted=True
        )

        # TODO: make half_cutoff (for coverage cutoff) True for X/Y?
        # This would also mean saving a new coverage model
        _, plateau_x_models = build_models(coverage_x_ht, args.trimers, weighted=True)
        _, plateau_y_models = build_models(coverage_y_ht, args.trimers, weighted=True)

        context_ht = context_ht.annotate(
            coverage_correction=get_coverage_correction_expr(
                context_ht.exome_coverage, coverage_model, args.high_cov_cutoff,
            )
        )
        context_ht = context_ht.annotate_globals(
            plateau_models=plateau_models,
            plateau_x_models=plateau_x_models,
            plateau_y_models=plateau_y_models,
            coverage_model=coverage_model,
        )
        obs_ht = calculate_observed(exome_ht, exac)

        if not args.skip_calc_exp:
            logger.info(
                "Creating autosomes-only, chrX non-PAR-only, and chrY non-PAR-only HT versions..."
            )
            context_x_ht = filter_to_region_type(context_ht, "chrX")
            context_y_ht = filter_to_region_type(context_ht, "chrY")
            context_auto_ht = filter_to_region_type(context_ht, "autosomes")

            logger.info("Calculating expected values...")
            exp_ht = calculate_expected(context_auto_ht, coverage_model, plateau_models)
            exp_x_ht = calculate_expected(
                context_x_ht, coverage_model, plateau_x_models
            )
            exp_y_ht = calculate_expected(
                context_y_ht, coverage_model, plateau_y_models
            )
            exp_ht = exp_ht.union(exp_x_ht).union(exp_y_ht)

            logger.info(
                "Annotating total observed and expected values and overall observed/expected value "
                "(capped at 1) per transcript..."
            )
            context_ht = context_ht.annotate(
                total_exp=exp_ht[context_ht.transcript].expected,
                total_obs=obs_ht[context_ht.transcript].observed,
            )
            context_ht = context_ht.annotate(
                overall_obs_exp=hl.min(context_ht.total_obs / context_ht.total_exp, 1)
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
        context_ht = context_ht.annotate(_obs=obs_ht.index(context_ht.key))
        context_ht = context_ht.transmute(
            observed=hl.int(hl.is_defined(context_ht._obs))
        )

        context_ht = calculate_exp_per_base(context_ht, args.groupings.split(","))

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
        "--high_cov_cutoff",
        help="Coverage threshold for a site to be considered high coverage",
        type=int,
        default=40,
    )
    parser.add_argument(
        "--pre_process_data", help="Pre-process data", action="store_true"
    )
    parser.add_argument(
        "--test",
        help="Filter to chr22 (for code testing purposes)",
        action="store_true",
    )
    parser.add_argument(
        "--skip_calc_exp",
        help="Skip observed and expected variant calculations per transcript. Relevant only to gnomAD v2.1.1!",
        action="store_true",
    )
    parser.add_argument(
        "--groupings",
        help="Fields to group by when calculating expected variants per base",
        default="context,ref,alt,cpg,methylation_level,mu_snp,transcript,exome_coverage",
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
