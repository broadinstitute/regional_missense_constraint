import logging
from typing import Dict, Tuple, Union

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


def get_cumulative_scan_expr(
    search_expr: hl.expr.StringExpression,
    observed_expr: hl.expr.Int64Expression,
    mu_expr: hl.expr.Float64Expression,
    prediction_flag: Tuple(float, int),
) -> hl.expr.StructExpression:
    """
    Creates struct with cumulative number of observed and expected variants.

    .. note::
        This function can produce the scan when searching for the first break or when searching for a second additional break.
            - When searching for the first break, this function should group by the transcript name (e.g., 'ENST00000255882').
            - When searching for an additional break, this function should group by the section of the transcript 
                (e.g., 'first' for before the first breakpoint or 'second' for after the first breakpoint).

    :param hl.expr.StringExpression search_expr: Expression containing transcript if searching for first break.
        Otherwise, expression containing transcript section if searching for second additional break.
    :param hl.expr.Float64Expression mu_expr: Mutation rate expression.
    :param hl.expr.Int64Expression observed_expr: Observed variants expression.
    :return: Struct containing the cumulative number of observed and expected variants.
    :param Tuple(float, int) prediction_flag: Adjustments to mutation rate based on chromosomal location
        (autosomes/PAR, X non-PAR, Y non-PAR). 
        E.g., prediction flag for autosomes/PAR is (0.4190964, 11330208)
    :return: Struct containing scan expressions for cumulative observed and expected variant counts.
    :rtype: hl.expr.StructExpression
    """
    return hl.struct(
        cumulative_observed=hl.scan.group_by(search_expr, hl.scan.sum(observed_expr)),
        cumulative_expected=hl.scan.group_by(
            search_expr, prediction_flag[0] + prediction_flag[1] * hl.scan.sum(mu_expr),
        ),
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


def get_reverse_obs_exp_expr(
    cond_expr: hl.expr.BooleanExpression,
    total_obs_expr: hl.expr.Int64Expression,
    total_exp_expr: hl.expr.Float64Expression,
    scan_obs_expr: Dict[hl.expr.StringExpression, hl.expr.Int64Expression],
    scan_exp_expr: Dict[hl.expr.StringExpression, hl.expr.Float64Expression],
) -> hl.expr.StructExpression:
    """
    Returns the "reverse" section observed and expected variant counts.

    The reverse counts are the counts moving from larger to smaller positions 
    (backwards from the end of the transcript back to the beginning of the transcript).
    reverse value = total value - cumulative value

    .. note::
        This function is designed to run on one transcript at a time.

    :param hl.expr.BooleanExpression cond_expr: Conditional expression to check before calculating reverse observed or expected value.
        Should be that the cumulative scan expression length isn't 0 when searching for the first break, or
        that the length of the cumulative scan expression length is 2 when searching for an additional break.
    :param hl.expr.Int64Expression total_obs_expr: Expression containing total number of observed variants for transcript.
    :param hl.expr.Float64Expression total_exp_expr: Expression containing total number of expected variants for transcript.
    :param Dict[hl.expr.StringExpression, hl.expr.Int64Expression] scan_obs_expr: Expression containing cumulative number of observed variants for transcript.
    :param Dict[hl.expr.StringExpression, hl.expr.Float64Expression] scan_expr_expr: Expression containing cumulative number of expected variants for transcript.
    :return: Struct with reverse observed and expected variant counts.
    :rtype: hl.expr.StructExpression
    """
    return hl.struct(
        obs=hl.or_missing(cond_expr, total_obs_expr - scan_obs_expr),
        exp=hl.or_missing(cond_expr, total_exp_expr - scan_exp_expr),
    )


def get_null_alt_expr(
    cond_expr: hl.expr.BooleanExpression,
    overall_oe_expr: hl.expr.Float64Expression,
    section_oe_expr: hl.expr.Float64Expression,
    obs_expr: Union[
        Dict[hl.expr.StringExpression, hl.expr.Int64Expression], hl.expr.Int64Expression
    ],
    exp_expr: Union[
        Dict[hl.expr.StringExpression, hl.expr.Float64Expression],
        hl.expr.Float64Expression,
    ],
) -> hl.expr.StructExpression:
    """
    Calculates null and alt values in preparation for chi-squared test to find significant breaks.

    All parameter values depend on the direction of calculation (forward/reverse) and 
    number of breaks (searching for first break or searching for additional break).

    For forward null/alts, values for obs_expr and and exp_expr should be:
        - Expression containing cumulative numbers for entire transcript.
        - Expression containing cumulative numbers for section of transcript 
            between the beginning or end of the transcript and the first breakpoint.
    For reverse null/alts, values for obs_expr and and exp_expr should be:
        - Reverse counts for entire transcript.
        - Reverse counts for section of transcript.
    
    For forward null/alts, values for overall_oe_expr and section_oe_expr should be:
        - Expression containing observed/expected value for entire transcript and
            expression containing observed/expected value calculated on cumulative observed and expected
            variants at each position.
        - Expression containing observed/expected value for section of transcript.
    For reverse null/alts, values for overall_oe_expr and section_oe_expr should be:
        - Expression containing observed/expected value for entire transcript and
            expression containing observed/expected value calculated on reverse observed variants value
            (total observed - cumulative observed count).
        - Expression containing observed/expected value for section of transcript and 
            expression containing reverse observed/expected value for section of transcript.

    For forward null/alts, cond_expr should check:
        - That the length of the obs_expr isn't 0 when searching for the first break.
        - That the length of the obs_expr is 2 when searching for a second additional break.
    For reverse null/alts, cond_expr should check:
        - That the reverse observed value for the entire transcript is defined when searching for the first break.
        - That the reverse observed value for the section between the first breakpoint and the end of the transcript
             is defined when searching for a second additional break.

    :param hl.expr.BooleanExpression cond_expr: Conditional expression to check before calculating null and alt values.
    :param hl.expr.Float64Expression overall_oe_expr: Expression of overall observed/expected value.
    :param hl.expr.Float64Expression section_oe_expr: Expression of section observed/expected value.
    :param Union[Dict[hl.expr.StringExpression, hl.expr.Int64Expression], hl.expr.Int64Expression] obs_expr: Expression containing observed variants count.
    :param Union[Dict[hl.expr.StringExpression, hl.expr.Float64Expression], hl.expr.Float64Expression] exp_expr: Expression containing expected variants count.
    :return: Struct containing forward or reverse null and alt values (either when searching for first or second break).
    :rtype: hl.expr.StructExpression
    """
    return hl.struct(
        null=hl.or_missing(cond_expr, hl.dpois(obs_expr, exp_expr * overall_oe_expr)),
        alt=hl.or_missing(cond_expr, hl.dpois(obs_expr, exp_expr * section_oe_expr)),
    )


def search_for_break(
    context_ht: hl.Table,
    obs_ht: hl.Table,
    exp_ht: hl.Table,
    search_str: str,
    prediction_flag: Tuple(float, int),
    chisq_threshold: float,
    group_by_transcript: bool,
) -> Union[hl.Table, None]:
    """
    Searches for breakpoints in a transcript. 

    Currently designed to run one transcript at a time.

    Expects context HT to contain the following fields:
        - locus
        - alleles
        - transcript or section
        - coverage (median)
        - mu
    Also expects:
        - multiallelic variants in context HT have been split.
        - context HT is autosomes/PAR only, X non-PAR only, or Y non-PAR only.

    Returns HT filtered to lines with maximum chisq if chisq >= max_value, otherwise returns None.

    :param hl.Table context_ht: Context Table. 
    :param hl.Table obs_ht: Table grouped by transcript with observed variant counts per transcript.
        Expects observed counts field to be named `observed`.
    :param hl.Table exp_ht: Table grouped by transcript with expected variant counts per transcript.
        Expects expected counts field to be named `expected`.
    :param str search_str: String that describes search space. Either transcript when searching for first break,
        or break section name (e.g., 'first', 'second') if searching for second additional break.
    :param Tuple(float, int) prediction_flag: Adjustments to mutation rate based on chromosomal location
        (autosomes/PAR, X non-PAR, Y non-PAR). 
        E.g., prediction flag for autosomes/PAR is (0.4190964, 11330208)
    :param float chisq_threshold: Chisq threshold for significance. 
        Value should be 10.8 (single break) and 13.8 (two breaks) (values from ExAC RMC code).
    :param bool group_by_transcript: Whether function should group by transcript. Should be True if searching 
        for first break and False to search for an additional break. 
    :return: Table filtered to rows with maximum chisq value IF max chisq is larger than chisq_threshold.
        Otherwise, returns None.
    :rtype: Union[hl.Table, None]
    """
    ht = annotate_observed_expected(context_ht, obs_ht, exp_ht, group_by_transcript)

    logger.info(
        "Annotating HT with cumulative expected/observed counts per transcript..."
    )
    ht = ht.annotate(
        scan_counts=get_cumulative_scan_expr(
            search_expr=ht[search_str],
            observed_expr=ht.observed,
            mu_expr=ht.mu,
            prediction_flag=prediction_flag,
        )
    )

    logger.info("Annotating HT with forward scan section observed/expected value...")
    # NOTE: Capping observed/expected values at 1
    ht = ht.annotate(
        obs_exp=hl.or_missing(
            hl.len(ht.scan_counts.cumulative_observed) != 0,
            hl.min(
                ht.scan_counts.cumulative_observed[search_str]
                / ht.scan_counts.cumulative_expected[search_str],
                1,
            ),
        ),
    )

    logger.info("Adding forward scan section nulls and alts...")
    # Add forwards sections (going through positions from smaller to larger)
    # section_null = stats.dpois(section_obs, section_exp*overall_obs_exp)[0]
    # section_alt = stats.dpois(section_obs, section_exp*section_obs_exp)[0]
    ht = ht.annotate(
        forward=get_null_alt_expr(
            cond_expr=hl.len(ht.cumulative_observed) != 0,
            overall_oe_expr=ht.overall_obs_exp,
            section_oe_expr=ht.obs_exp,
            obs_expr=ht.cumulative_observed[search_str],
            exp_expr=ht.cumulative_expected[search_str],
        )
    )

    logger.info("Adding reverse section observeds and expecteds...")
    # reverse value = total value - cumulative value
    ht = ht.annotate(
        reverse_counts=get_reverse_obs_exp_expr(
            cond_expr=hl.len(ht.cumulative_observed) != 0,
            total_obs_expr=ht.total_obs,
            total_exp_expr=ht.total_exp,
            scan_obs_expr=ht.cumulative_observed[search_str],
            scan_exp_expr=ht.cumulative_expected[search_str],
        )
    )

    # Set reverse o/e to missing if reverse expected value is 0 (to avoid NaNs)
    # Also cap reverse observed/expected at 1
    ht = ht.annotate(
        reverse_obs_exp=hl.or_missing(
            ht.reverse_counts.exp != 0,
            hl.min(ht.reverse_counts.obs / ht.reverse_counts.exp, 1),
        )
    )

    logger.info("Adding reverse section nulls and alts...")
    ht = ht.annotate(
        reverse=get_null_alt_expr(
            cond_expr=hl.len(ht.scan_counts.cumulative_observed) != 0,
            overall_oe_expr=ht.overall_obs_exp,
            section_oe_expr=ht.obs_exp,
            obs_expr=ht.scan_counts.cumulative_observed[search_str],
            exp_expr=ht.scan_counts.cumulative_expected[search_str],
        )
    )

    logger.info("Multiplying all section nulls and all section alts...")
    # Kaitlin stores all nulls/alts in section_null and section_alt and then multiplies
    # e.g., p1 = prod(section_null_ps)
    ht = ht.annotate(
        section_null=ht.forward.null * ht.reverse.null,
        section_alt=ht.forward.alt * ht.reverse.alt,
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
