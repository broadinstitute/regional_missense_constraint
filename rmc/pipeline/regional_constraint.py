import argparse
import logging

import hail as hl

from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.slack import slack_notifications
from gnomad_lof.constraint_utils.generic import prepare_ht
from rmc.resources.basics import logging_path
from rmc.resources.grch37.exac import filtered_exac
from rmc.resources.grch37.gnomad import filtered_exomes, processed_exomes
from rmc.resources.grch37.reference_data import processed_context
from rmc.slack_creds import slack_token
from rmc.utils.constraint import calculate_expected, calculate_observed
from rmc.utils.generic import (
    filter_to_region_type,
    filter_to_missense,
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
            exome_ht = prepare_ht(exome_ht, args.trimers)
            exome_ht = processed_exomes.ht()

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

        logger.info("Reading in context HT and gencode HT...")
        context_ht = processed_context.ht()

        if args.test:
            logger.info("Inferring build of exome HT...")
            rg = get_reference_genome(exome_ht.locus)

            logger.info("Filtering to chr22 for testing...")
            contigs = rg.contigs[21]
            context_ht = hl.filter_intervals(
                context_ht, [hl.parse_locus_interval(contigs, reference_genome=rg)]
            )
            exome_ht = hl.filter_intervals(
                exome_ht, [hl.parse_locus_interval(contigs, reference_genome=rg)]
            )

        else:
            logger.info(
                "Creating autosomes-only, chrX non-PAR-only, and chrY non-PAR-only HT versions..."
            )
            context_x_ht = filter_to_region_type(context_ht, "chrX")
            context_y_ht = filter_to_region_type(context_ht, "chrY")
            context_ht = filter_to_region_type(context_ht, "autosomes")

            exome_x_ht = filter_to_region_type(exome_ht, "chrX")
            exome_y_ht = filter_to_region_type(exome_ht, "chrY")
            exome_ht = filter_to_region_type(exome_ht, "chrX")

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script searches for regional missense constraint in gnomAD"
    )

    parser.add_argument(
        "--pre_process_data", help="Pre-process data", action="store_true"
    )
    parser.add_argument(
        "--calc_exp", help="Calculate expected variant counts", action="store_true"
    )
    parser.add_argument(
        "--calc_obs", help="Calculated observed variant counts", action="store_true"
    )

    parser.add_argument(
        "--trimers", help="Use trimers instead of heptamers", action="store_false"
    )
    parser.add_argument(
        "--exac", help="Use ExAC Table (not gnomAD Table)", action="store_true"
    )
    parser.add_argument(
        "--test",
        help="Filter to chr22 (for code testing purposes)",
        action="store_true",
    )
    parser.add_argument(
        "--pre_process_data", help="Pre-process data", action="store_true"
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
