import argparse
import logging

import hail as hl

from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.slack import slack_notifications
from gnomad_lof.constraint_utils.generic import prepare_ht
from rmc.resources.basics import logging_path
from rmc.resources.grch37.exac import coverage, filtered_exac_cov
from rmc.resources.grch37.gnomad import filtered_exomes, processed_exomes
from rmc.resources.grch37.reference_data import processed_context, processed_gencode
from rmc.slack_creds import slack_token
from rmc.utils.constraint import calculate_expected, calculate_observed
from rmc.utils.generic import (
    filter_alt_decoy,
    filter_to_missense,
    process_context_ht,
    process_gencode_ht,
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

            logger.info("Preprocessing GENCODE GTF information...")
            process_gencode_ht("GRCh37")

            logger.info(
                "Filtering gnomAD exomes HT to missense variants in canonical transcripts only..."
            )
            exome_ht = processed_exomes.ht()
            exome_ht = filter_to_missense(exome_ht)
            exome_ht.write(filtered_exomes.path, overwrite=args.overwrite)

            logger.info("Done preprocessing files")

        logger.info("Reading in exome HT...")
        if exac:
            exome_ht = filtered_exac_cov.ht().select("context", "coverage")
            coverage_ht = coverage.ht()
            csq = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info|context|ancestral"
            csq = csq.split("|")

            # Move necessary annotations out of info struct and into top level annotations
            exome_ht = exome_ht.transmute(
                ac=exome_ht.info.AC_Adj[exome_ht.a_index - 1],
                VQSLOD=exome_ht.info.VQSLOD,
                coverage=exome_ht.coverage.median,
                canonical=exome_ht.info.CSQ[exome_ht.a_index - 1].split("\|")[
                    csq.index("CANONICAL")
                ]
                == "YES",
                transcript=exome_ht.info.CSQ[exome_ht.a_index - 1].split("\|")[
                    csq.index("Feature")
                ],
                exon=exome_ht.info.CSQ[exome_ht.a_index - 1].split("\|")[csq.index("EXON")],
            )

            logger.info("Filtering to canonical transcripts only...")
            # NOTE: filter_vep_to_canonical_transcripts doesn't work because of vep format in ExAC HT (in csq format from vcf export)
            # Searching for 'YES' based on this line of code from vep.py: "canonical": hl.cond(element.canonical == 1, "YES", "")
            exome_ht = exome_ht.filter(exome_ht.canonical)

        else:
            exome_ht = prepare_ht(filtered_exomes.ht(), args.trimers)
            coverage_ht = None

        logger.info("Reading in context HT and gencode HT...")
        context_ht = processed_context.ht()
        #gencode_ht = processed_gencode.ht()

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
                "Removing alt/decoy contig calls and annotating with region type..."
            )
            context_ht = filter_alt_decoy(context_ht)
            exome_ht = filter_alt_decoy(exome_ht)

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
