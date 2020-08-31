import argparse
import logging

import hail as hl

from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.slack import slack_notifications
from rmc.resources.grch37.exac import (
    coverage,
    exac,
    filtered_exac,
)
from rmc.slack_creds import slack_token


CSQ = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info|context|ancestral"
"""
VEP CSQ format taken from ExAC release 1.0 VCF header.
"""

KEEP_FIELDS = [
    "context",
    "ac",
    "vqslod",
    "transcript",
    "exon",
    "coverage",
    "a_index",
    "was_split",
]
"""
Fields to select from the ExAC HT.
"""


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("prepare_exac")
logger.setLevel(logging.INFO)


def filter_to_missense(ht: hl.Table, csq: str = CSQ) -> hl.Table:
    """
    Filter ExAC Table to missense variants in canonical transcripts.

    :param hl.Table ht: ExAC Table.
    :param str csq: String with format of VEP CSQ field. Default is CSQ.
    :return: Filtered ExAC Table. 
    :rtype: hl.Table
    """
    logger.info("Splitting multiallelic variants and filtering to SNPs...")
    ht = hl.split_multi(ht)
    ht = ht.filter(hl.is_snp(ht.alleles[0], ht.alleles[1]))

    logger.info("Adding canonical, transcript, exon, and consequence information...")
    csq = csq.split("|")
    ht = ht.annotate(
        canonical=ht.info.CSQ[ht.a_index - 1].split("\|")[csq.index("CANONICAL")]
        == "YES",
        transcript=ht.info.CSQ[ht.a_index - 1].split("\|")[csq.index("Feature")],
        exon=ht.info.CSQ[ht.a_index - 1].split("\|")[csq.index("EXON")],
        consequence=ht.info.CSQ[ht.a_index - 1].split("\|")[csq.index("Consequence")],
    )

    # NOTE: Using most severe consequence didn't work properly when testing in a notebook
    logger.info("Filtering to missense variants in canonical transcripts only...")
    return ht.filter((ht.canonical) & (ht.consequence.contains("missense_variant")))


def main(args):

    hl.init(log="/prepare_exac.log")

    logger.info("Filtering ExAC HT to missense variants on chr22...")
    ht = exac.ht()
    ht = hl.filter_intervals(ht, [hl.parse_locus_interval("22")])
    ht = filter_to_missense(ht)

    # Move necessary annotations out of info struct and into top level annotations
    # Also add coverage annotation
    coverage_ht = coverage.ht()
    ht = ht.transmute(
        ac=ht.info.AC_Adj[ht.a_index - 1],
        vqslod=ht.info.VQSLOD,
        coverage=coverage_ht[ht.locus],
    )

    # Add context bases
    rg = get_reference_genome(ht.locus, add_sequence=True)
    ht = ht.annotate(
        context=hl.get_sequence(
            ht.locus.contig, ht.locus.position, before=1, after=1, reference_genome=rg,
        )
    ).select(*KEEP_FIELDS)
    ht.naive_coalesce(args.n_partitions).write(
        filtered_exac.path, overwrite=args.overwrite
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script prepares the ExAC sites vcf for RMC testing"
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
