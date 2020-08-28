import logging

import hail as hl

# from gnomad_lof.constraint_utils.generic import count_variants
from rmc.utils.generic import keep_criteria


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)


def calculate_expected(context_ht: hl.Table, coverage_ht: hl.Table) -> hl.Table:
    """
    Returns table of transcripts and the total number of expected variants per transcript.

    Expected variants count is adjusted by mutation rate, divergence score, and region type.

    .. note::
        This function is currently ExAC-specific.
        TODO: look into whether can reuse models from gnomad_lof for gnomAD v2 exomes

    :param Table context_ht: Context Table.
    :return: Table grouped by transcript with expected variant counts per transcript.
    :rtype: hl.Table
    """
    logger.info("Annotating context HT with median ExAC coverage...")
    context_ht = context_ht.annotate(exac_cov=coverage_ht[context_ht.key].median)

    logger.info("Adjusting mutation rate with divergence scores...")
    # p_mut2 = p_mut1*(1 + (0.31898*div_score))
    context_ht = context_ht.annotate(
        mu=context_ht.mu_snp * (1 + (0.31898 * context_ht.div)) * 0.9822219
    )

    logger.info("Grouping by transcript...")
    context_ht = context_ht.group_by(
        context_ht.transcript, context_ht.region_type, context_ht.coverage
    ).aggregate(raw_expected=hl.agg.sum(context_ht.mu))

    logger.info("Processing context HT to calculate number of expected variants...")
    # Annotate prediction flag based on context HT locus type
    # Region type annotation added when filtering alt and decoy contigs
    context_ht = context_ht.transmute(
        expected=hl.case()
        .when(
            context_ht.region_type == "x_nonpar",
            (0.3715167 + 7796945 * context_ht.raw_expected),
        )
        .when(
            context_ht.region_type == "y",
            (0.05330181 + 2457366 * context_ht.raw_expected),
        )
        .default(0.4190964 + 11330208 * context_ht.raw_expected)
    )

    # Adjust expected counts based on depth
    context_ht = context_ht.transmute(
        expected=hl.case()
        .when(context_ht.coverage < 1, context_ht.expected * 0.089)
        .when(
            ((context_ht.coverage >= 1) & (context_ht.coverage < 50)),
            (0.089 + 0.217 * hl.log(context_ht.coverage)) * context_ht.expected,
        )
        .default(context_ht.expected)
    )

    context_ht = context_ht.group_by("transcript").aggregate(
        expected=hl.agg.sum(context_ht.expected)
    )
    return context_ht


def calculate_observed(ht: hl.Table, exac: bool) -> hl.Table:
    """
    Groups input Table by transcript and aggregates observed variants count per transcript.

    :param Table ht: Input Table.
    :param bool exac: Whether the input Table is ExAC data.
    :return: Table annotated with observed variant counts.
    :rtype: hl.Table
    """
    ht = ht.filter(keep_criteria(ht, exac))
    return ht.group_by(ht.transcript).aggregate(observed=hl.agg.count())
