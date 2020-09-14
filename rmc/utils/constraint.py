import logging
from typing import Tuple, Union

import hail as hl

# from gnomad_lof.constraint_utils.generic import count_variants
from rmc.utils.generic import keep_criteria


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)


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

    return context_ht.group_by("transcript").aggregate(
        expected=hl.agg.sum(context_ht.expected)
    )


def annotate_observed_expected(
    ht: hl.Table, obs_ht: hl.Table, exp_ht: hl.Table, group_by_transcript: bool = True,
) -> hl.Table:
    """
    Annotates input Table with total number of observed and expected variants.

    Groups Table by transcript or section of transcript. 
    Also annotates Table with overall observed/expected value.

    .. note::
        - obs_ht will be grouped by section of transcript only after finding first breakpoint in transcript.
        - Expects that keys of input ht and obs_ht match.

    :param hl.Table ht: Input Table. 
    :param hl.Table obs_ht: Table grouped by transcript with observed variant counts per transcript.
        Expects observed counts field to be named `observed`.
    :param hl.Table exp_ht: Table grouped by transcript with expected variant counts per transcript.
        Expects expected counts field to be named `expected`.
    :param bool group_by_transcript: Whether the observed Table is grouped by transcript. Default is True.
        If False, expects that observed Table is grouped by section of transcript.
    :return: Table annotated with observed, expected, and observed/expected values.
    :rtype: hl.Table.
    """
    logger.info("Annotating HT with observed counts...")
    ht = ht.annotate(_obs=obs_ht.index(ht.key, all_matches=True))
    ht = ht.transmute(observed=hl.int(hl.is_defined(ht._obs)))

    logger.info("Annotating HT with total expected/observed counts...")
    if group_by_transcript:
        group_expr = ht.transcript
    else:
        group_expr = ht.section
    ht = ht.annotate(
        total_exp=exp_ht[group_expr].expected, total_obs=obs_ht[group_expr].observed,
    )

    logger.info(
        "Annotating overall observed/expected value (capped at 1) and returning HT..."
    )
    return ht.annotate(overall_obs_exp=hl.min(ht.total_obs / ht.total_exp, 1))


def get_cumulative_scan_expr(
    transcript_expr: hl.expr.StringExpression,
    mu_expr: hl.expr.Float64Expression,
    observed_expr: hl.expr.Int64Expression,
    prediction_flag: Tuple(float, int),
) -> hl.expr.StructExpression:
    """
    Creates struct with cumulative number of observed and expected variants.

    :param hl.expr.StringExpression transcript_expr: Transcript expression.
    :param hl.expr.Float64Expression mu_expr: Mutation rate expression.
    :param hl.expr.Int64Expression observed_expr: Observed variants expression.
    :return: Struct containing the cumulative number of observed and expected variants.
    :param Tuple(float, int) prediction_flag: Adjustments to mutation rate based on chromosomal location
        (autosomes/PAR, X non-PAR, Y non-PAR). 
        E.g., prediction flag for autosomes/PAR is (0.4190964, 11330208)
    :rtype: hl.expr.StructExpression
    """
    # TODO: Create functions for nulls/alt + reverse null/alt
    return hl.struct(
        cumulative_expected=hl.scan.group_by(
            transcript_expr,
            prediction_flag[0] + prediction_flag[1] * hl.scan.sum(mu_expr),
        ),
        cumulative_observed=hl.scan.group_by(
            transcript_expr, hl.scan.sum(observed_expr)
        ),
    )


def search_for_break(
    context_ht: hl.Table,
    obs_ht: hl.Table,
    exp_ht: hl.Table,
    transcript: str,
    prediction_flag: Tuple(float, int),
    chisq_threshold: float,
) -> Union[hl.Table, None]:
    """
    Searches for breakpoints in a transcript. 

    Currently designed for one transcript at a time.

    Expects context HT to contain the following fields:
        - locus
        - alleles
        - transcript
        - coverage (median)
        - mu
    Also expects:
        - multiallelic variants in context HT have been split
        - context HT is autosomes/PAR only, X non-PAR only, or Y non-PAR only.

    Returns HT filtered to lines with maximum chisq if chisq >= max_value, otherwise returns None.

    :param hl.Table context_ht: Context Table. 
    :param hl.Table obs_ht: Table grouped by transcript with observed variant counts per transcript.
        Expects observed counts field to be named `observed`.
    :param hl.Table exp_ht: Table grouped by transcript with expected variant counts per transcript.
        Expects expected counts field to be named `expected`.
    :param str transcript: Transcript of interest.
    :param Tuple(float, int) prediction_flag: Adjustments to mutation rate based on chromosomal location
        (autosomes/PAR, X non-PAR, Y non-PAR). 
        E.g., prediction flag for autosomes/PAR is (0.4190964, 11330208)
    :param float chisq_threshold: Chisq threshold for significance. 
        Value should be 10.8 (single break) and 13.8 (two breaks) (values from ExAC RMC code).
    :return: Table filtered to rows with maximum chisq value IF max chisq is larger than chisq_threshold.
        Otherwise, returns None.
    :rtype: Union[hl.Table, None]
    """
    ht = annotate_observed_expected(context_ht, obs_ht, exp_ht)

    logger.info(
        "Annotating HT with cumulative expected/observed counts per transcript..."
    )
    ht = ht.annotate(
        cumulative_expected=hl.scan.group_by(
            ht.transcript, prediction_flag[0] + prediction_flag[1] * hl.scan.sum(ht.mu)
        ),
        cumulative_observed=hl.scan.group_by(ht.transcript, hl.scan.sum(ht.observed)),
    )

    logger.info("Annotating HT with forward scan section observed/expected value...")
    # NOTE: Capping observed/expected values at 1
    ht = ht.annotate(
        obs_exp=hl.or_missing(
            hl.len(ht.cumulative_observed) != 0,
            hl.min(
                ht.cumulative_observed[transcript] / ht.cumulative_expected[transcript],
                1,
            ),
        ),
    )

    logger.info("Adding forward scan section nulls and alts...")
    # Add forwards sections (going through positions from smaller to larger)
    # section_null = stats.dpois(section_obs, section_exp*overall_obs_exp)[0]
    # section_alt = stats.dpois(section_obs, section_exp*section_obs_exp)[0]
    ht = ht.annotate(
        null=hl.or_missing(
            hl.len(ht.cumulative_observed) != 0,
            hl.dpois(
                ht.cumulative_observed[transcript],
                ht.cumulative_expected[transcript] * ht.overall_obs_exp,
            ),
        ),
        alt=hl.or_missing(
            hl.len(ht.cumulative_observed) != 0,
            hl.dpois(
                ht.cumulative_observed[transcript],
                ht.cumulative_expected[transcript] * ht.obs_exp,
            ),
        ),
    )

    logger.info("Adding reverse section observeds and expecteds...")
    # reverse value = total value - cumulative value
    ht = ht.annotate(
        reverse_obs=hl.or_missing(
            hl.len(ht.cumulative_observed) != 0,
            ht.total_obs - ht.cumulative_observed[transcript],
        ),
        reverse_exp=hl.or_missing(
            hl.len(ht.cumulative_expected) != 0,
            ht.total_exp - ht.cumulative_expected[transcript],
            # 0
        ),
    )

    # Set reverse o/e to missing if reverse expected value is 0 (to avoid NaNs)
    # Also cap reverse observed/expected at 1
    ht = ht.annotate(
        reverse_obs_exp=hl.or_missing(
            ht.reverse_exp != 0, hl.min(ht.reverse_obs / ht.reverse_exp, 1),
        )
    )

    logger.info("Adding reverse section nulls and alts...")
    ht = ht.annotate(
        reverse_null=hl.or_missing(
            hl.is_defined(ht.reverse_obs),
            hl.dpois(ht.reverse_obs, ht.reverse_exp * ht.overall_obs_exp),
        ),
        reverse_alt=hl.or_missing(
            hl.is_defined(ht.reverse_obs),
            hl.dpois(ht.reverse_obs, ht.reverse_exp * ht.reverse_obs_exp),
        ),
    )

    logger.info("Multiplying all section nulls and all section alts...")
    # Kaitlin stores all nulls/alts in section_null and section_alt and then multiplies
    # e.g., p1 = prod(section_null_ps)
    ht = ht.annotate(
        section_null=ht.null * ht.reverse_null, section_alt=ht.alt * ht.reverse_alt,
    )

    logger.info("Adding chisq value and getting max chisq...")
    ht = ht.annotate(chisq=(2 * (hl.log(ht.section_alt) - hl.log(ht.section_null))))

    # "The default chi-squared value for one break to be considered significant is
    # 10.8 (p ~ 10e-3) and is 13.8 (p ~ 10e-4) for two breaks. These currently cannot
    # be adjusted."
    max_chisq = ht.aggregate(hl.agg.max(ht.chisq))
    if max_chisq >= chisq_threshold:
        return ht.filter(ht.chisq == max_chisq)

    return None
