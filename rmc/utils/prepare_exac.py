import argparse
import logging

import hail as hl

from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import CSQ_ORDER
from rmc.resources.grch37.exac import (
    coverage,
    exac,
    filtered_exac,
    filtered_exac_cov,
)
from rmc.resources.resource_utils import MISSENSE
from rmc.slack_creds import slack_token


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("prepare_exac")
logger.setLevel(logging.INFO)


def filter_to_missense(ht: hl.Table) -> hl.Table:
    """
    Filter ExAC ht to missense variants

    :param Table ht: ExAC ht
    :return: ExAC ht filtered to missense variants
    :rtype: Table
    """
    logger.info(f"HT count before filtration: {ht.count()}")

    csqs = hl.literal(CSQ_ORDER)
    ht = ht.explode(ht.info.CSQ)
    ht = ht.annotate(
        most_severe_consequence=csqs.find(lambda c: ht.info.CSQ.contains(c))
    )
    logger.info(
        f"Consequence counts:{ht.aggregate(hl.agg.counter(ht.most_severe_consequence))}"
    )

    logger.info("Filtering to missense variants...")
    ht = ht.filter(hl.literal(MISSENSE).contains(ht.most_severe_consequence))
    logger.info(f"HT count after filtration: {ht.count()}")
    logger.info("Deduplicating keys...")
    ht = ht.distinct()
    logger.info(f"HT count after dedup: {ht.count()}")
    return ht


def main(args):

    hl.init(log="/prepare_exac.log")

    if args.filter_ht:
        logger.info("Filtering ExAC ht to only missense variants")
        ht = exac.ht()
        ht = filter_to_missense(ht)
        rg = get_reference_genome(ht.locus, add_sequence=True)
        ht = ht.annotate(
            context=hl.get_sequence(
                ht.locus.contig,
                ht.locus.position,
                before=1,
                after=1,
                reference_genome=rg,
            )
        )
        ht.naive_coalesce(args.n_partitions).write(
            filtered_exac.path, overwrite=args.overwrite
        )

    if args.join_cov:
        ht = filtered_exac.ht()
        coverage_ht = coverage.ht()
        ht = ht.annotate(coverage=coverage_ht[ht.locus])
        ht = ht.naive_coalesce(args.n_partitions)
        ht.write(filtered_exac_cov.path, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script prepares the ExAC sites vcf for RMC testing"
    )

    parser.add_argument(
        "--filter_ht", help="Filter ExAC ht to missense variants", action="store_true"
    )
    parser.add_argument(
        "--join_cov", help="Annotate ExAC ht with coverage", action="store_true"
    )
    parser.add_argument(
        "--n_partitions",
        help="Desired number of partitions for output HTs",
        default=500,
        type=int,
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
