import argparse
import logging

import hail as hl

from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.slack import slack_notifications
from gnomad_lof.constraint_utils.generic import count_variants, prepare_ht
from rmc.resources.basics import logging_path
from rmc.resources.grch37.exac import coverage, filtered_exac_cov
from rmc.resources.grch37.gnomad import filtered_exomes, processed_exomes
from rmc.resources.grch37.reference_data import processed_context, processed_gencode
from rmc.slack_creds import slack_token
from rmc.utils.utils import (
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


def _keep_criteria(ht: hl.Table, exac: bool) -> hl.expr.BooleanExpression:
    """
    Returns Boolean expression to filter variants in input Table.

    :param hl.Table ht: Input Table.
    :param bool exac: Whether input Table is ExAC data.
    :return: Keep criteria Boolean expression.
    :rtype: hl.expr.BooleanExpression
    """
    # ExAC keep criteria: adjusted AC <= 123 and VQSLOD >= -2.632
    # Also remove variants with median depth < 1
    if exac:
        keep_criteria = (
            (ht.ac <= 123) & (ht.ac > 0) & (ht.VQSLOD >= -2.632) & (ht.coverage > 1)
        )
    else:
        # TODO: check about impose_high_af_cutoff upfront
        keep_criteria = (ht.ac > 0) & (ht.pass_filters)

    return keep_criteria


def calculate_expected(context_ht: hl.Table, coverage_ht: hl.Table,) -> hl.Table:
    """
    Annotates context Table with expected variants count (adjusted by mutation rate, divergence score, and region type).

    .. note::
        This function is currently ExAC-specific.
        TODO: look into whether can reuse models from gnomad_lof for gnomAD v2 exomes

    :param Table context_ht: Context Table.
    :return: Context Table with expected variant annotations.
    :rtype: hl.Table
    """
    logger.info("Annotating context HT with median ExAC coverage...")
    context_ht = context_ht.annotate(exac_cov=coverage_ht[context_ht.key].median)

    logger.info("Adjusting mutation rate with divergence scores...")
    # p_mut2 = p_mut1*(1 + (0.31898*div_score))
    context_ht = context_ht.annotate(
        mu=context_ht.mu_snp * (1 + (0.31898 * context_ht.div)) * 0.9822219
    )

    logger.info("Processing context HT to calculate number of expected variants...")
    # Annotate prediction flag based on context HT locus type
    # Region type annotation added when filtering alt and decoy contigs
    context_ht = context_ht.annotate(
        raw_expected=hl.case()
        .when(
            context_ht.region_type == "x_nonpar", (0.3715167 + 7796945 * context_ht.mu)
        )
        .when(context_ht.region_type == "y", (0.05330181 + 2457366 * context_ht.mu))
        .default(0.4190964 + 11330208 * context_ht.mu)
    )

    # Adjust expected counts based on depth
    context_ht = context_ht.transmute(
        expected=hl.case()
        .when(context_ht.exac_cov < 1, 0.089 * context_ht.raw_expected)
        .when(
            (context_ht.exac_cov >= 1) & (context_ht.exac_cov < 50),
            (0.089 + 0.217 * hl.log(context_ht.raw_expected)),
        )
        .default(context_ht.raw_expected)
    )
    return context_ht


def calculate_observed(
    ht: hl.Table,
    gencode_ht: hl.Table,
    scan: bool,
    exac: bool,
    omit_methylation: bool,
    n_partitions: int,
) -> hl.Table:
    """
    Groups input Table by transcript and exon and calculates observed variants count.

    If `scan` is true, annotates input Table with an array containing cumulative variant counts.
    Otherwise, annotates input Table with observed variants counts per exon.

    :param Table ht: Input Table.
    :param Table gencode_ht: Gencode Table grouped by transcript and exon number.
    :param bool scan: Whether to perform a scan on the Table.
    :param bool exac: Whether the input Table is ExAC data.
    :param bool omit_methylation: Whether to omit grouping the Table by methylation. Must be true if `exac` is True.
    :param int n_partitions: Number of partitions used during the `group_by` operation.
    :return: Table annotated with observed variant counts.
    :rtype: hl.Table
    """
    ht = ht.filter(_keep_criteria(ht, exac))

    # Create dictionary of additional values to group variants by
    # Need these additional groupings to do scans correctly downstream
    groupings = {
        "transcript": ht.transcript,
        "exon": ht.exon,
        "chrom": ht.locus.contig,
        "pos": ht.locus.position,
        "coverage": ht.coverage.median,
    }
    ht = ht.annotate(**groupings)

    # Reformat exon number annotation on HT
    ht = ht.transmute(exon_number=ht.exon.split("\/")[0])

    # NOTE: count variants from gnomAD lof repo
    # https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/generic.py#L68
    # Need force_grouping to be true to return a Table and not a struct
    if exac:
        omit_methylation = True
    ht = count_variants(
        ht,
        omit_methylation=omit_methylation,
        additional_grouping=list(groupings),
        partition_hint=n_partitions,
        force_grouping=True,
    )
    if scan:
        # Group HT by transcript and exon and collect variants + variant counts into an array of structs
        ht = ht.group_by(ht.transcript, ht.exon,).aggregate(
            variant_info=hl.agg.collect(
                hl.struct(chrom=ht.chrom, pos=ht.chrom, variant_count=ht.variant_count,)
            )
        )
        # Sort the array of structs by chrom, then by position 
        ht = ht.transmute(variant_info=hl.sorted(ht.variant_info, key=lambda x: (x[0], x[1])))
        ht = ht.annotate(
            scan_sum=hl.array_scan(lambda i, j: i + j, 0, ht.variant_info.variant_count)
        )

    else:
        # Group HT by transcript and exon and get overall counts per exon
        ht = ht.group_by(ht.transcript, ht.exon,).aggregate(
            observed=hl.agg.sum(ht.variant_count)
        )

    # Left join the observed counts onto the gencode HT and return
    # TODO: exons with no observed counts have missing values after the join. should I set to 0?
    return gencode_ht.join(ht.key_by("transcript", "exon_number"), how="left")


def main(args):

    hl.init(log="/RMC.log")
    exac = args.exac

    try:
        if args.pre_process_data:
            logger.info("Preprocessing reference fasta and gencode files...")
            process_context_ht("GRCh37", args.trimers)

            logger.info("Preprocessing GENCODE GTF information...")
            process_gencode_ht("GRCh37")

            logger.info("Filtering gnomAD exomes HT to missense variants only...")
            exome_ht = processed_exomes.ht()
            # TODO: make filter to missense filter to canonical transcripts only
            exome_ht = filter_to_missense(exome_ht)
            exome_ht.write(filtered_exomes.path, overwrite=args.overwrite)

            logger.info("Done preprocessing files")

        logger.info("Reading in exome HT...")
        if exac:
            exome_ht = filtered_exac_cov.ht().select("context", "coverage")
            coverage_ht = coverage.ht()
            csq = "Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|SYMBOL|SYMBOL_SOURCE|HGNC_ID|BIOTYPE|CANONICAL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SIFT|PolyPhen|EXON|INTRON|DOMAINS|HGVSc|HGVSp|GMAF|AFR_MAF|AMR_MAF|ASN_MAF|EUR_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF_info|LoF_flags|LoF_filter|LoF".split(
                "|"
            )

            # Move necessary annotations out of info struct and into top level annotations
            exome_ht = exome_ht.transmute(
                ac=exome_ht.info.AC_Adj[0],
                VQSLOD=exome_ht.info.VQSLOD,
                coverage=exome_ht.coverage.median,
                canonical=exome_ht.info.CSQ.split("\|")[csq.index("CANONICAL")]
                == "YES",
                transcript=exome_ht.info.CSQ.split("\|")[csq.index("Feature")],
                exon=exome_ht.info.CSQ.split("\|")[csq.index("EXON")],
                gene=exome_ht.info.CSQ.split("\|")[csq.index("Gene")],
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
        gencode_ht = processed_gencode.ht()

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
