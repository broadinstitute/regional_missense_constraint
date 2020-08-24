import hail as hl

from gnomad_lof.constraint_utils.generic import count_variants
from rmc.utils.generic import keep_criteria


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)


def calculate_expected(context_ht: hl.Table, coverage_ht: hl.Table) -> hl.Table:
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
    exac: bool,
    omit_methylation: bool,
    n_partitions: int,
) -> hl.Table:
    """
    Groups input Table by transcript and exon and calculates observed variants count.

    Annotates input Table with an array containing cumulative variant counts.
    Also adds annotation with total number of observed variants per exon.

    :param Table ht: Input Table.
    :param Table gencode_ht: Gencode Table grouped by transcript and exon number.
    :param bool exac: Whether the input Table is ExAC data.
    :param bool omit_methylation: Whether to omit grouping the Table by methylation. Must be true if `exac` is True.
    :param int n_partitions: Number of partitions used during the `group_by` operation.
    :return: Table annotated with observed variant counts.
    :rtype: hl.Table
    """
    ht = ht.filter(keep_criteria(ht, exac))

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

    # Group HT by transcript and exon and collect variants + variant counts into an array of structs
    ht = ht.group_by(ht.transcript, ht.exon,).aggregate(
        variant_info=hl.agg.collect(
            hl.struct(chrom=ht.chrom, pos=ht.chrom, variant_count=ht.variant_count,)
        )
    )
    # Sort the array of structs by chrom, then by position
    ht = ht.transmute(
        variant_info=hl.sorted(ht.variant_info, key=lambda x: (x[0], x[1]))
    )
    ht = ht.annotate(
        scan_sum=hl.array_scan(lambda i, j: i + j, 0, ht.variant_info.variant_count)
    )

    # Add overall count per exon to HT
    ht = ht.annotate(observed=ht.scan_sum[-1])

    # Left join the observed counts onto the gencode HT and return
    # TODO: exons with no observed counts have missing values after the join. should I set to 0?
    return gencode_ht.join(ht.key_by("transcript", "exon_number"), how="left")
