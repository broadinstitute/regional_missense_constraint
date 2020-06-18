import argparse
import logging

import hail as hl

from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import CSQ_ORDER
from rmc.resources.basics import (
    exac_cov_path,
    exac_ht,
    exac_tsv_path,
    exac_vcf,
    filt_exac_cov_ht,
    filtered_exac_ht,
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

    if args.import_vcf:
        logger.info("Importing ExAC VCF")
        ht = hl.import_vcf(
            exac_vcf, force_bgz=True, min_partitions=args.min_partitions
        ).rows()
        ht = ht.naive_coalesce(args.n_partitions).write(
            exac_ht, overwrite=args.overwrite
        )

    if args.filter_ht:
        logger.info("Filtering ExAC ht to only missense variants")
        ht = hl.read_table(exac_ht)
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
            filtered_exac_ht, overwrite=args.overwrite
        )

    # chroms = ['X', 'Y']
    # for i in range(1, 23):
    #    chroms.append[i]
    chroms = ["22"]

    # NOTE: only calculated for chr22
    if args.import_cov:
        for chrom in chroms:
            tsv = f"{exac_tsv_path}/Panel.chr{chrom}.coverage.txt.gz"
            out = f"{exac_cov_path}/{chrom}_coverage.ht"
            ht = hl.import_table(
                tsv, min_partitions=args.min_partitions, impute=True, force_bgz=True
            )
            ht = ht.transmute(
                locus=hl.parse_locus(hl.format("%s:%s", ht["#chrom"], ht.pos))
            )
            ht = ht.rename(
                {
                    "1": "over_1",
                    "5": "over_5",
                    "10": "over_10",
                    "15": "over_15",
                    "20": "over_20",
                    "25": "over_25",
                    "30": "over_30",
                    "50": "over_50",
                    "100": "over_100",
                }
            )
            ht = ht.key_by("locus")
            ht = ht.naive_coalesce(args.n_partitions)
            ht.write(out, overwrite=args.overwrite)

    if args.join_cov:
        ht = hl.read_table(filtered_exac_ht)
        for chrom in chroms:
            cov_path = f"{exac_cov_path}/{chrom}_coverage.ht"
            cov_ht = hl.read_table(cov_path)
            ht = ht.annotate(coverage=cov_ht[ht.locus])
        ht = ht.naive_coalesce(args.n_partitions)
        ht.write(filt_exac_cov_ht, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script prepares the ExAC sites vcf for RMC testing"
    )

    parser.add_argument(
        "--import_vcf", help="Import ExAC VCF and write to ht", action="store_true"
    )
    parser.add_argument(
        "--filter_ht", help="Filter ExAC ht to missense variants", action="store_true"
    )
    parser.add_argument(
        "--import_cov", help="Import coverage files", action="store_true"
    )
    parser.add_argument(
        "--join_cov", help="Annotate ExAC ht with coverage", action="store_true"
    )
    parser.add_argument(
        "--min_partitions",
        help="Minimum number of partitions for imported data",
        default=100,
        type=int,
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
