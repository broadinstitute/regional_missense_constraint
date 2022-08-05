import logging
from typing import Dict, List, Set, Tuple, Union

import hail as hl


from rmc.resources.basics import (
    multiple_breaks,
    rmc_browser,
    rmc_results,
    temp_path,
)


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)


GROUPINGS = [
    "context",
    "ref",
    "alt",
    "methylation_level",
    "annotation",
    "modifier",
    "transcript",
    "gene",
    "exome_coverage",
]
"""
Core fields to group by when calculating expected variants per variant type.

Fields taken from gnomAD LoF repo.
"""

CONSTRAINT_ANNOTATIONS = [
    "mu_snp",
    "total_exp",
    "_mu_scan",
    "total_mu",
    "cumulative_obs",
    "observed",
    "cumulative_exp",
    "total_obs",
    "reverse",
    "forward_oe",
    "overall_oe",
    "start_pos",
    "end_pos",
    "transcript_size",
]
"""
List of annotations required to calculate constraint.

Used to drop unnecessary fields when searching for simultaneous breaks.
"""

FINAL_ANNOTATIONS = [
    "mu_snp",
    "observed",
    "total_exp",
    "total_mu",
    "total_obs",
    "cumulative_obs",
    "_mu_scan",
    "cumulative_exp",
    "break_list",
    "start_pos",
    "end_pos",
]
"""
List of annotations to keep when finalizing release HT.
"""


def calculate_observed(ht: hl.Table) -> hl.Table:
    """
    Group input Table by transcript, filter based on `keep_criteria`, and aggregate observed variants count per transcript.

    .. note::
        Assumes input HT has been filtered using `keep_criteria`.

    :param hl.Table ht: Input Table.
    :return: Table annotated with observed variant counts.
    :rtype: hl.Table
    """
    return ht.group_by(ht.transcript_consequences.transcript_id).aggregate(
        observed=hl.agg.count()
    )


def get_cumulative_mu_expr(
    transcript_expr: hl.expr.StringExpression,
    mu_expr: hl.expr.Float64Expression,
) -> hl.expr.DictExpression:
    """
    Return annotation with the cumulative mutation rate probability, shifted by one.

    Value is shifted by one due to the nature of `hl.scan` and needs to be corrected later.

    This function can produce the scan when searching for the first break or when searching for additional break(s).

    :param hl.expr.StringExpression transcript_expr: StringExpression containing transcript information.
    :param hl.expr.Float64Expression mu_expr: FloatExpression containing mutation rate probability per variant.
    :return: DictExpression containing scan expressions for cumulative mutation rate probability.
    :rtype: hl.expr.DictExpression
    """
    return hl.scan.group_by(transcript_expr, hl.scan.sum(mu_expr))


def adjust_mu_expr(
    cumulative_mu_expr: hl.expr.DictExpression,
    mu_expr: hl.expr.Int32Expression,
    transcript_expr: hl.expr.StringExpression,
) -> hl.expr.DictExpression:
    """
    Adjust the scan with the cumulative number mutation rate probability.

    This adjustment is necessary because scans are always one line behind, and we want the values to match per line.

    This function can correct the scan created when searching for the first break or when searching for additional break(s).

    :param hl.expr.DictExpression cumulative_mu_expr: DictExpression containing scan expressions for cumulative mutation rate probability.
    :param hl.expr.Float64Expression mu_expr: FloatExpression containing mutation rate probability per variant.
    :param hl.expr.StringExpression transcript_expr: StringExpression containing transcript information.
    :return: Adjusted cumulative mutation rate expression.
    :rtype: hl.expr.DictExpression
    """
    return hl.if_else(
        # Check if the current transcript/section exists in the _mu_scan dictionary
        # If it doesn't exist, that means this is the first line in the HT for that particular transcript/section
        # The first line of a scan is always missing, but we want it to exist
        # Thus, set the cumulative_mu equal to the current mu_snp value
        hl.is_missing(cumulative_mu_expr.get(transcript_expr)),
        {transcript_expr: mu_expr},
        # Otherwise, add the current mu_snp to the scan to make sure the cumulative value isn't one line behind
        {transcript_expr: cumulative_mu_expr[transcript_expr] + mu_expr},
    )


def translate_mu_to_exp_expr(
    cumulative_mu_expr: hl.expr.DictExpression,
    transcript_expr: hl.expr.StringExpression,
    total_mu_expr: hl.expr.Float64Expression,
    total_exp_expr: hl.expr.Float64Expression,
) -> hl.expr.DictExpression:
    """
    Translate cumulative mutation rate probability into cumulative expected count per base.

    Expected variants counts are produced per base by first calculating the fraction of probability of mutation per base,
    then multiplying that fraction by the total expected variants count for a transcript or transcript sub-section.

    The expected variants count for section of interest is mutation rate per SNP adjusted by location in the genome/CpG status
    (plateau model) and coverage (coverage model).

    :param hl.expr.DictExpression cumulative_mu_expr: DictExpression containing scan expressions for cumulative mutation rate probability.
    :param hl.expr.StringExpression transcript_expr: StringExpression containing transcript information.
    :param hl.expr.Float64Expression total_mu_expr: FloatExpression describing total sum of mutation rate probabilities per transcript.
    :param hl.expr.Float64Expression total_exp_expr: FloatExpression describing total expected variant counts per transcript.
    :return: Cumulative expected variants count expression.
    :rtype: hl.expr.DictExpression
    """
    return (cumulative_mu_expr[transcript_expr] / total_mu_expr) * total_exp_expr


def get_obs_exp_expr(
    cond_expr: hl.expr.BooleanExpression,
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
) -> hl.expr.Float64Expression:
    """
    Return observed/expected annotation based on inputs.

    Cap observed/expected (OE) value at 1. This is to avoid pulling out regions that are enriched for missense variation.
    Code in this pipeline is looking for missense constraint, so regions with an OE >= 1.0 can be grouped together.

    Function can generate observed/expected values across the entire transcript or section of a transcript depending on inputs.
    Function can also generate 'forward' (moving from smaller to larger positions") or 'reverse' (moving from larger to smaller positions)
    section obs/exp values.

    .. note::
        `cond_expr` should vary depending on size/direction of section being annotated.

    :param hl.expr.BooleanExpression cond_expr: Condition to check prior to adding obs/exp expression.
    :param hl.expr.Int64Expression obs_expr: Expression containing number of observed variants.
    :param hl.expr.Float64Expression exp_expr: Expression containing number of expected variants.
    :return: Observed/expected expression.
    :rtype: hl.expr.Float64Expression
    """
    return hl.or_missing(cond_expr, hl.min(obs_expr / exp_expr, 1))


def get_cumulative_obs_expr(
    transcript_expr: hl.expr.StringExpression,
    observed_expr: hl.expr.Int64Expression,
) -> hl.expr.DictExpression:
    """
    Return annotation with the cumulative number of observed variants, shifted by one.

    Value is shifted by one due to the nature of `hl.scan` and needs to be corrected later.
    This function can produce the scan when searching for the first break or when searching for additional break(s).
    :param hl.expr.StringExpression transcript_expr: StringExpression containing transcript information.
    :param hl.expr.Int64Expression observed_expr: Observed variants expression.
    :return: Struct containing the cumulative number of observed and expected variants.
    :return: DictExpression containing scan expressions for cumulative observed variant counts for `search_expr`.
    :rtype: hl.expr.DictExpression
    """
    return hl.scan.group_by(transcript_expr, hl.scan.sum(observed_expr))


def adjust_obs_expr(
    cumulative_obs_expr: hl.expr.DictExpression,
    obs_expr: hl.expr.Int64Expression,
    transcript_expr: hl.expr.StringExpression,
) -> hl.expr.DictExpression:
    """
    Adjust the scan with the cumulative number of observed variants.

    This adjustment is necessary because scans are always one line behind, and we want the values to match per line.

    This function can correct the scan created when searching for the first break or when searching for additional break(s).

    .. note::
        This function expects:
            - cumulative_obs_expr is a DictExpression keyed by transcript.

    :param hl.expr.DictExpression cumulative_obs_expr: DictExpression containing scan expression with cumulative observed counts per base.
    :param hl.expr.Int32Expression obs_expr: IntExpression with value of either 0 (no observed variant at site) or 1 (variant found in gnomAD).
    :param hl.expr.StringExpression transcript_expr: StringExpression containing transcript information.
    :return: Adjusted cumulative observed counts expression.
    :rtype: hl.expr.DictExpression
    """
    return hl.if_else(
        # Check if the current transcript/section exists in the _obs_scan dictionary
        # If it doesn't exist, that means this is the first line in the HT for that particular transcript
        # The first line of a scan is always missing, but we want it to exist
        # Thus, set the cumulative_obs equal to the current observed value
        hl.is_missing(cumulative_obs_expr.get(transcript_expr)),
        {transcript_expr: obs_expr},
        # Otherwise, add the current obs to the scan to make sure the cumulative value isn't one line behind
        {transcript_expr: cumulative_obs_expr[transcript_expr] + obs_expr},
    )


def get_reverse_obs_exp_expr(
    total_obs_expr: hl.expr.Int64Expression,
    total_exp_expr: hl.expr.Float64Expression,
    scan_obs_expr: hl.expr.DictExpression,
    cumulative_exp_expr: hl.expr.Float64Expression,
) -> hl.expr.StructExpression:
    """
    Return the "reverse" section observed and expected variant counts.

    The reverse counts are the counts moving from larger to smaller positions
    (backwards from the end of the transcript back to the beginning of the transcript).
    reverse value = total value - cumulative value

    :param hl.expr.Int64Expression total_obs_expr: Expression containing total number of observed variants for transcript.
    :param hl.expr.Float64Expression total_exp_expr: Expression containing total number of expected variants for transcript.
    :param hl.expr.DictExpression scan_obs_expr: Expression containing cumulative number of observed variants for transcript.
    :param hl.expr.Float64Expression cumulative_exp_expr: Expression containing cumulative number of expected variants for transcript.
    :return: Struct with reverse observed and expected variant counts.
    :rtype: hl.expr.StructExpression
    """
    return hl.struct(
        # NOTE: Adding hl.max to exp expression to make sure reverse exp is never negative
        # Without this, ran into errors where reverse exp was -5e-14
        # Picked 1e-09 here as tiny number that is not 0
        # ExAC code also did not allow reverse exp to be zero, as this breaks the likelihood ratio tests
        obs=total_obs_expr - scan_obs_expr,
        exp=hl.max(total_exp_expr - cumulative_exp_expr, 1e-09),
    )


def get_fwd_exprs(
    ht: hl.Table,
    transcript_str: str,
    obs_str: str,
    mu_str: str,
    total_mu_str: str,
    total_exp_str: str,
) -> hl.Table:
    """
    Annotate input Table with the forward section cumulative observed, expected, and observed/expected values.

    .. note::
        'Forward' refers to moving through the transcript from smaller to larger chromosomal positions.

    Expects:
        - Input HT is annotated with transcript, observed, mutation rate, total mutation rate (per section),
        and total expected counts (per section).

    :param hl.Table ht: Input Table.
    :param str transcript_str: Name of field containing transcript information.
    :param str obs_str: Name of field containing observed variants counts.
    :param str mu_str: Name of field containing mutation rate probability per variant.
    :param str total_mu_str: Name of field containing total mutation rate per section of interest (transcript or sub-section of transcript).
    :param str total_exp_str: Name of field containing total expected variants count per section of interest (transcript or sub-section of transcript).
    :param str search_field: Name of field to group by prior to running scan. Should be 'transcript' if searching for the first break.
        Otherwise, should be transcript section if searching for additional breaks.
    :return: Table with forward values (cumulative obs, exp, and forward o/e) annotated.
    :rtype: hl.Table
    """
    logger.info("Getting cumulative observed variant counts...")
    ht = ht.annotate(
        _obs_scan=get_cumulative_obs_expr(
            transcript_expr=ht[transcript_str],
            observed_expr=ht[obs_str],
        )
    )
    ht = ht.annotate(
        cumulative_obs=adjust_obs_expr(
            cumulative_obs_expr=ht._obs_scan,
            obs_expr=ht[obs_str],
            transcript_expr=ht[transcript_str],
        )
    )

    logger.info("Getting cumulative expected variant counts...")
    # Get scan of mu_snp
    ht = ht.annotate(_mu_scan=get_cumulative_mu_expr(ht[transcript_str], ht[mu_str]))
    # Adjust scan of mu_snp
    ht = ht.annotate(
        _mu_scan=adjust_mu_expr(ht._mu_scan, ht[mu_str], ht[transcript_str])
    )
    ht = ht.annotate(
        cumulative_exp=translate_mu_to_exp_expr(
            ht._mu_scan, ht[transcript_str], ht[total_mu_str], ht[total_exp_str]
        )
    )

    logger.info("Getting forward observed/expected count and returning...")
    # NOTE: adding cond_expr here because get_obs_exp_expr expects it
    # cond_expr is necessary for reverse obs/exp, which is why the function has it
    ht = ht.annotate(cond_expr=True)
    ht = ht.annotate(
        forward_oe=get_obs_exp_expr(
            ht.cond_expr,
            ht.cumulative_obs[ht[transcript_str]],
            ht.cumulative_exp,
        )
    )
    return ht.drop("cond_expr")


def get_reverse_exprs(
    ht: hl.Table,
    total_obs_expr: hl.expr.Int64Expression,
    total_exp_expr: hl.expr.Float64Expression,
    scan_obs_expr: Dict[hl.expr.StringExpression, hl.expr.Int64Expression],
    scan_exp_expr: Dict[hl.expr.StringExpression, hl.expr.Float64Expression],
) -> hl.Table:
    """
    Call `get_reverse_obs_exp_expr` and `get_obs_exp_expr` to add the reverse section cumulative observed, expected, and observed/expected values.

    .. note::
        'Reverse' refers to moving through the transcript from larger to smaller chromosomal positions.

    :param hl.Table ht: Input Table.
    :param hl.expr.Int64Expression total_obs_expr: Expression containing total number of observed variants per transcript (if searching for first break)
        or per section (if searching for additional breaks).
    :param hl.expr.Float64Expression total_exp_expr: Expression containing total number of expected variants per transcript (if searching for first break)
        or per section (if searching for additional breaks).
    :param Dict[hl.expr.StringExpression, hl.expr.Int64Expression] scan_obs_expr: Expression containing cumulative number of observed variants per transcript
        (if searching for first break) or per section (if searching for additional breaks).
    :param Dict[hl.expr.StringExpression, hl.expr.Float64Expression] scan_exp_expr: Expression containing cumulative number of expected variants per transcript
        (if searching for first break) or per section (if searching for additional breaks).
    :return: Table with reverse values annotated
    :rtype: hl.Table
    """
    # reverse value = total value - cumulative value
    ht = ht.annotate(
        reverse=get_reverse_obs_exp_expr(
            total_obs_expr=total_obs_expr,
            total_exp_expr=total_exp_expr,
            scan_obs_expr=scan_obs_expr,
            cumulative_exp_expr=scan_exp_expr,
        )
    )

    # Set reverse o/e to missing if reverse expected value is 0 (to avoid NaNs)
    return ht.annotate(
        reverse_obs_exp=get_obs_exp_expr(
            (ht.reverse.exp != 0), ht.reverse.obs, ht.reverse.exp
        )
    )


def get_dpois_expr(
    cond_expr: hl.expr.BooleanExpression,
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
    Calculate null and alt values in preparation for chi-squared test to find significant breaks.

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
    :param hl.expr.Float64Expression section_oe_expr: Expression of section observed/expected value.
    :param Union[Dict[hl.expr.StringExpression, hl.expr.Int64Expression], hl.expr.Int64Expression] obs_expr: Expression containing observed variants count.
    :param Union[Dict[hl.expr.StringExpression, hl.expr.Float64Expression], hl.expr.Float64Expression] exp_expr: Expression containing expected variants count.
    :return: Struct containing forward or reverse null and alt values (either when searching for first or second break).
    :rtype: hl.expr.StructExpression
    """
    return hl.or_missing(cond_expr, hl.dpois(obs_expr, exp_expr * section_oe_expr))


def get_section_expr(
    dpois_expr: hl.expr.ArrayExpression,
) -> hl.expr.Float64Expression:
    """
    Build null or alt model by multiplying all section null or alt distributions.

    For example, when checking for the first break in a transcript, the transcript is broken into two sections:
    pre-breakpoint and post-breakpoint. Each section's null and alt distributions must be multiplied
    to later test the significance of the break.

    :param hl.expr.ArrayExpression dpois_expr: ArrayExpression that contains all section nulls or alts.
        Needs to be reduced to single float by multiplying each number in the array.
    :return: Overall section distribution.
    :rtype: hl.expr.Float64Expression
    """
    return hl.fold(lambda i, j: i * j, 1, dpois_expr)


def search_for_break(
    ht: hl.Table,
    search_field: hl.str,
    break_num: int,
    chisq_threshold: float = 10.8,
) -> hl.Table:
    """
    Search for breakpoints in a transcript or within a transcript subsection.

    Expects input context HT to contain the following fields:
        - locus
        - alleles
        - transcript or section
        - coverage (median)
        - mu_snp
        - scan_counts struct
        - overall_oe
        - forward_oe
        - reverse struct
        - reverse_obs_exp
        - total_obs
        - total_exp
    If searching for simultaneous breaks, expects HT to have the following fields:
        - pre_obs
        - pre_exp
        - pre_oe
        - window_obs
        - window_exp
        - window_oe
        - post_obs
        - post_exp
        - post_oe
    Also expects:
        - multiallelic variants in input HT have been split.
        - Input HT is autosomes/PAR only, X non-PAR only, or Y non-PAR only.

    Return HT filtered to lines with maximum chisq if chisq >= max_value, otherwise returns None.

    :param hl.Table ht: Input context Table.
    :param hl.expr.StringExpression search_field: Field of table to search. Value should be either 'transcript' or 'section'.
    :param float chisq_threshold: Chi-square significance threshold.
        Value should be 10.8 (single break) and 13.8 (two breaks) (values from ExAC RMC code).
        Default is 10.8.
    :param int break_num: Break search round number (Temporary addition for additional breaks check)
    :return: Table annotated with whether position is a breakpoint.
    :rtype: hl.Table
    """
    logger.info(
        "Creating section null (no regional variability in missense depletion)\
        and alt (evidence of domains of missense constraint) expressions..."
    )

    # Split transcript or transcript subsection into two sections
    # Split transcript when searching for first break
    # Split transcript subsection when searching for additional breaks
    ht = ht.annotate(
        section_nulls=[
            # Add forwards section null (going through positions from smaller to larger)
            # section_null = stats.dpois(section_obs, section_exp*overall_obs_exp)[0]
            get_dpois_expr(
                cond_expr=hl.len(ht.cumulative_obs) != 0,
                section_oe_expr=ht.overall_oe,
                obs_expr=ht.cumulative_obs[ht[search_field]],
                exp_expr=ht.cumulative_exp,
            ),
            # Add reverse section null (going through positions larger to smaller)
            get_dpois_expr(
                cond_expr=hl.is_defined(ht.reverse.obs),
                section_oe_expr=ht.overall_oe,
                obs_expr=ht.reverse.obs,
                exp_expr=ht.reverse.exp,
            ),
        ],
        section_alts=[
            # Add forward section alt
            # section_alt = stats.dpois(section_obs, section_exp*section_obs_exp)[0]
            get_dpois_expr(
                cond_expr=hl.len(ht.cumulative_obs) != 0,
                section_oe_expr=ht.forward_oe,
                obs_expr=ht.cumulative_obs[ht[search_field]],
                exp_expr=ht.cumulative_exp,
            ),
            # Add reverse section alt
            get_dpois_expr(
                cond_expr=hl.is_defined(ht.reverse.obs),
                section_oe_expr=ht.reverse_obs_exp,
                obs_expr=ht.reverse.obs,
                exp_expr=ht.reverse.exp,
            ),
        ],
    )

    logger.info("Multiplying all section nulls and all section alts...")
    # Kaitlin stores all nulls/alts in section_null and section_alt and then multiplies
    # e.g., p1 = prod(section_null_ps)
    ht = ht.annotate(
        total_null=get_section_expr(ht.section_nulls),
        total_alt=get_section_expr(ht.section_alts),
    )

    logger.info("Adding chisq value and getting max chisq...")
    ht = ht.annotate(chisq=(2 * (hl.log10(ht.total_alt) - hl.log10(ht.total_null))))

    # "The default chi-squared value for one break to be considered significant is
    # 10.8 (p ~ 10e-3) and is 13.8 (p ~ 10e-4) for two breaks. These currently cannot
    # be adjusted."
    group_ht = ht.group_by(search_field).aggregate(max_chisq=hl.agg.max(ht.chisq))
    group_ht = group_ht.checkpoint(
        f"{temp_path}/break_{break_num}_max_chisq.ht", overwrite=True
    )
    ht = ht.annotate(max_chisq=group_ht[ht.transcript].max_chisq)
    return ht.annotate(
        is_break=((ht.chisq == ht.max_chisq) & (ht.chisq >= chisq_threshold))
    )


def process_transcripts(ht: hl.Table, chisq_threshold: float):
    """
    Annotate each position in Table with whether that position is a significant breakpoint.

    Also annotate input Table with cumulative observed, expected, and observed/expected values
    for both forward (moving from smaller to larger positions) and reverse (moving from larger to
    smaller positions) directions.

    Expects input Table to have the following fields:
        - locus
        - alleles
        - mu_snp
        - transcript
        - observed
        - expected
        - coverage_correction
        - cpg
        - cumulative_obs
        - cumulative_exp
        - forward_oe
    Also expects:
        - multiallelic variants in context HT have been split.
        - Global annotations contains plateau models.

    :param hl.Table ht: Input Table. annotated with observed and expected variants counts per transcript.
    :param float chisq_threshold: Chi-square significance threshold.
        Value should be 10.8 (single break) and 13.8 (two breaks) (values from ExAC RMC code).
    :return: Table with cumulative observed, expected, and observed/expected values annotated for forward and reverse directions.
        Table also annotated with boolean for whether each position is a breakpoint.
    :rtype: hl.Table
    """
    logger.info(
        "Annotating HT with reverse observed and expected counts for each transcript...\n"
        "(transcript-level reverse (moving from larger to smaller positions) values)"
    )
    ht = get_reverse_exprs(
        ht=ht,
        total_obs_expr=ht.total_obs,
        total_exp_expr=ht.total_exp,
        scan_obs_expr=ht.cumulative_obs[ht.transcript],
        scan_exp_expr=ht.cumulative_exp,
    )

    return search_for_break(ht, "transcript", chisq_threshold)


def get_subsection_exprs(
    ht: hl.Table,
    section_str: str = "section",
    obs_str: str = "observed",
    mu_str: str = "mu_snp",
    total_mu_str: str = "total_mu",
    total_exp_str: str = "total_exp",
) -> hl.Table:
    """
    Annotate total observed, expected, and observed/expected (OE) counts for each section of a transcript.

    .. note::
        Assumes input Table is annotated with:
            - section
            - observed variants count per site
            - mutation rate probability per site
            - total mutation rate probability per transcript
            - total expected variant counts per transcript
        Names of annotations must match section_str, obs_str, mu_str, total_mu_str, and total_exp_str.

    :param hl.Table ht: Input Table.
    :param str section_str: Name of section annotation.
    :param str obs_str: Name of observed variant counts annotation.
    :param str mu_str: Name of mutation rate probability per site annotation.
    :param str total_mu_str: Name of annotation containing sum of mutation rate probabilities per transcript.
    :param str total_exp_str: Name of annotation containing total expected variant counts per transcript.
    :return: Table annotated with section observed, expected, and OE counts.
    :return: hl.Table
    """
    logger.info(
        "Getting total observed and expected counts for each transcript section..."
    )
    # Get total obs and mu per section
    section_counts = ht.group_by(ht[section_str]).aggregate(
        obs=hl.agg.sum(ht[obs_str]),
        mu=hl.agg.sum(ht[mu_str]),
    )

    # Translate total mu to total expected per section
    ht = ht.annotate(
        section_exp=(section_counts[ht[section_str]].mu / ht[total_mu_str])
        * ht[total_exp_str],
        section_obs=section_counts[ht[section_str]].obs,
    )

    logger.info("Getting observed/expected value for each transcript section...")
    return ht.annotate(
        break_oe=get_obs_exp_expr(
            cond_expr=hl.is_defined(ht[section_str]),
            obs_expr=ht.section_obs,
            exp_expr=ht.section_exp,
        )
    )


def process_sections(ht: hl.Table, chisq_threshold: float, break_num: int):
    """
    Search for breaks within given sections of a transcript.

    Expects that input Table has the following annotations:
        - context
        - ref
        - alt
        - cpg
        - observed
        - mu_snp
        - coverage_correction
        - methylation_level
        - section

    Also assumes that Table's globals contain plateau and coverage models.

    :param hl.Table ht: Input Table.
    :param float chisq_threshold: Chi-square significance threshold.
        Value should be 10.8 (single break) and 13.8 (two breaks) (values from ExAC RMC code).
    :param int break_num: Break search round number (Temporary addition for additional breaks check)
    :return: Table annotated with whether position is a breakpoint.
    :rtype: hl.Table
    """
    ht = get_subsection_exprs(ht)

    # Rename break_oe (each section's observed/expected value) to be overall_oe
    # This is because the overall OE ratio used in searching for additional breaks should be
    # transcript section OE and not transcript overall OE
    ht = ht.transmute(overall_oe=ht.break_oe)

    logger.info("Splitting HT into pre and post breakpoint sections...")
    pre_ht = ht.filter(ht.section.contains("pre"))
    post_ht = ht.filter(ht.section.contains("post"))

    logger.info(
        "Annotating post breakpoint HT with cumulative observed and expected counts..."
    )
    post_ht = get_fwd_exprs(
        ht=post_ht,
        transcript_str="transcript",
        obs_str="observed",
        mu_str="mu_snp",
        total_mu_str="total_mu",
        total_exp_str="total_exp",
    )

    logger.info(
        "Annotating pre and post breakpoint HTs with reverse observed and expected counts..."
    )
    pre_ht = get_reverse_exprs(
        ht=pre_ht,
        total_obs_expr=pre_ht.section_obs,
        total_exp_expr=pre_ht.section_exp,
        scan_obs_expr=pre_ht.cumulative_obs[pre_ht.transcript],
        scan_exp_expr=pre_ht.cumulative_exp,
    )
    post_ht = get_reverse_exprs(
        ht=post_ht,
        total_obs_expr=post_ht.section_obs,
        total_exp_expr=post_ht.section_exp,
        scan_obs_expr=post_ht.cumulative_obs[post_ht.transcript],
        scan_exp_expr=post_ht.cumulative_exp,
    )

    logger.info("Searching for a break in each section and returning...")
    pre_ht = search_for_break(
        pre_ht,
        search_field="transcript",
        break_num=break_num,
        chisq_threshold=chisq_threshold,
    )
    post_ht = search_for_break(
        post_ht,
        search_field="transcript",
        break_num=break_num,
        chisq_threshold=chisq_threshold,
    )
    # Adjust is_break annotation in both HTs
    # to prevent this function from continually finding previous significant breaks
    pre_ht = pre_ht.annotate(
        is_break=hl.if_else(pre_ht.break_list.any(lambda x: x), False, pre_ht.is_break)
    )
    post_ht = post_ht.annotate(
        is_break=hl.if_else(
            post_ht.break_list.any(lambda x: x), False, post_ht.is_break
        )
    )
    return pre_ht.union(post_ht)


def process_additional_breaks(
    ht: hl.Table, break_num: int, chisq_threshold: float
) -> hl.Table:
    """
    Search for additional breaks in a transcript after finding one significant break.

    Expects that input Table is filtered to only transcripts already containing one significant break.
    Expects that input Table has the following annotations:
        - cpg
        - observed
        - mu_snp
        - coverage_correction
        - transcript

    Also assumes that Table's globals contain plateau models.

    :param hl.Table ht: Input Table.
    :param int break_num: Number of additional break: 2 for second break, 3 for third, etc.
    :param float chisq_threshold: Chi-square significance threshold.
        Value should be 10.8 (single break) and 13.8 (two breaks) (values from ExAC RMC code).
    :return: Table annotated with whether position is a breakpoint.
    :rtype: hl.Table
    """
    logger.info(
        "Generating table keyed by transcripts (used to get breakpoint position later)..."
    )
    break_ht = ht.filter(ht.is_break).key_by("transcript")

    logger.info(
        "Renaming is_break field to prepare to search for an additional break..."
    )
    # Rename because this will be overwritten when searching for additional break
    annot_expr = {"is_break": f"is_break_{break_num - 1}"}
    ht = ht.rename(annot_expr)

    logger.info(
        "Splitting each transcript into two sections: pre-first breakpoint and post..."
    )
    ht = ht.annotate(
        section=hl.if_else(
            ht.locus.position > break_ht[ht.transcript].locus.position,
            hl.format("%s_%s", ht.transcript, "post"),
            hl.format("%s_%s", ht.transcript, "pre"),
        ),
    )
    return process_sections(ht, chisq_threshold, break_num)


def calculate_section_chisq(
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
) -> hl.expr.Float64Expression:
    """
    Create expression checking if transcript section is significantly different than the null model (no evidence of regional missense constraint).

    Formula is: (section obs - section exp)^2 / section exp. Taken from ExAC RMC code.

    :param hl.expr.Int64Expression obs_expr:
    :param hl.expr.Float64Expression exp_expr:
    :return: Transcript section chi-squared value.
    :rtype: hl.expr.Float64Expression
    """
    return ((obs_expr - exp_expr) ** 2) / exp_expr


def get_all_breakpoint_pos(ht: hl.Table) -> hl.GroupedTable:
    """
    Get all breakpoint positions per transcript.

    :param hl.Table ht: Input Table.
            Example schema (truncated at two breaks for space reasons):
            ----------------------------------------
            Row fields:
                'locus': locus<GRCh37>
                'transcript': str
                'mu_snp': float64
                'observed': int32
                'total_exp': float64
                'total_mu': float64
                'total_obs': int64
                'max_chisq': float64
                'break_list': array<bool>
                'break_1_null': float64
                'break_1_alt': float64
                'break_1_chisq': float64
                'break_2_chisq': float64
                'break_2_max_chisq': float64
                'break_2_null': float64
                'break_2_alt': float64
                'break_pos': array<int32>
            ----------------------------------------
            Key: ['locus', 'transcript']
            ----------------------------------------
    :return: Table grouped by transcript, with all breakpoint positions annotated as a list.
    :rtype: hl.GroupedTable
    """
    ht = ht.filter(ht.break_list.any(lambda x: x))
    return ht.group_by("transcript").aggregate(
        break_pos=hl.sorted(hl.agg.collect(ht.locus.position))
    )


def get_section_info(
    ht: hl.Table,
    section_num: int,
    section_type: str,
    indices: Tuple[int],
) -> hl.Table:
    """
    Get the number of observed variants, number of expected variants, and chi square value for transcript section.

    .. note::
        Assumes that the input Table is annotated with a list of breakpoint positions (`break_pos`) and
        with each transcript's start and end positions (`start_pos`, `end_pos`).

    :param hl.Table ht: Input Table.
    :param int section_num: Transcript section number (e.g., 1 for first section, 2 for second, 3 for third, etc.).
    :param str section_type: Transcript section type. Must be one of 'first', 'middle', or 'end'.
    :param Tuple[int] indices: List of indices pointing to breakpoints.
    :return: Table containing transcript and new section obs, exp, and chi square annotations.
    """
    assert section_type in {
        "first",
        "middle",
        "last",
    }, "section_type must be one of 'first', 'middle', 'last'!"

    logger.info("Getting info for section of transcript between two breakpoints...")
    if section_type == "first":
        logger.info(
            "Getting info for first section of transcript (up to and including smallest breakpoint pos)..."
        )
        ht = ht.filter(ht.locus.position <= ht.break_pos[0])
        ht = ht.annotate(
            # Start position is transcript start if this is the first section
            section_start_pos=ht.start_pos,
            section_end_pos=ht.break_pos[0],
        )
    elif section_type == "middle":
        ht = ht.filter(
            (ht.locus.position > ht.break_pos[indices[0]])
            & (ht.locus.position <= ht.break_pos[indices[1]])
        )

        # Add section start/end position annotations
        ht = ht.annotate(
            # Add 1 to break_pos[indices[0]] since it isn't actually included in the region
            section_start_pos=ht.break_pos[indices[0]] + 1,
            section_end_pos=ht.break_pos[indices[1]],
        )

    else:
        logger.info(
            "Getting info for last section of transcript (after largest breakpoint pos)..."
        )
        ht = ht.filter(ht.locus.position > ht.break_pos[-1])
        ht = ht.annotate(
            # Get last position from break_pos list and add 1 since it isn't included in the region
            # Use the transcript end position as the end position since this is the last section
            section_start_pos=ht.break_pos[-1] + 1,
            section_end_pos=ht.end_pos,
        )

    ht = ht.annotate(section=hl.format("%s_%s", ht.transcript, str(section_num)))
    ht = get_subsection_exprs(
        ht, "transcript", "observed", "mu_snp", "total_mu", "total_exp"
    )
    return ht.annotate(
        section_chisq=calculate_section_chisq(ht.section_obs, ht.section_exp)
    )


def annotate_transcript_sections(
    ht: hl.Table,
    max_n_breaks: int,
) -> hl.Table:
    """
    Annotate each transcript section with observed, expected, OE, and section chi square values.

    .. note::
        Needs to be run for each break number. For example, this function needs to be run for transcripts with three breaks,
        and it needs to be run again for transcripts with four breaks.

    :param hl.Table ht: Input Table.
    :param int max_n_breaks: Largest number of breaks.
    :return: Table with section and section values annotated.
    :rtype: hl.Table
    """
    logger.info("Get section information for first section of each transcript...")
    count = 1
    section_ht = get_section_info(
        ht,
        section_num=count,
        section_type="first",
        indices=None,
    )

    # Check sections between breakpoint positions
    while count <= max_n_breaks:
        if count == 1:
            # One break transcripts only get divided into two sections
            # Thus, increment counter and continue
            count += 1
            continue
        else:
            temp_ht = get_section_info(
                ht,
                section_num=count,
                section_type="middle",
                indices=(count - 2, count - 1),
            )
            section_ht = section_ht.union(temp_ht)
            count += 1
    end_ht = get_section_info(ht, section_num=count, section_type="last", indices=None)
    return section_ht.union(end_ht)


def get_unique_transcripts_per_break(
    ht: hl.Table,
    max_n_breaks: int,
) -> Dict[int, Union[Set[str], hl.expr.SetExpression]]:
    """
    Return the set of transcripts unique to each break number.

    If the set is empty, return an empty SetExpression.

    .. note::
        - This function will only get unique transcripts for transcripts with one or one + additional breaks.
        - This will not work for transcripts with two simultaneous breaks.
        - Assumes input Table is annotated with list containing booleans for whether that locus is a breakpoint
        (`break_list`).

    :param hl.Table: Input Table.
    :param int max_n_breaks: Largest number of breaks.
    :return: Dictionary with break number (key) and set of transcripts unique to that break number or empty SetExpression (value).
    :rtype: Dict[int, Union[Set[str], hl.expr.SetExpression]]
    """
    transcripts_per_break = {}
    # Use `break_list` annotation (list of booleans for whether a row is a breakpoint)
    # to filter one/additional breaks ht to rows that are significant breakpoints ONLY
    ht = ht.filter(ht.break_list.any(lambda x: x))

    # Group HT (filtered to breakpoint positions only) by transcript and
    # count the number of breakpoints associated with each transcript
    group_ht = ht.group_by("transcript").aggregate(n_breaks=hl.agg.count())

    # Checkpoint to force hail to finish this group by computation
    group_ht = group_ht.checkpoint(
        f"{temp_path}/breaks_per_transcript.ht", overwrite=True
    )

    for i in range(1, max_n_breaks + 1):
        temp_ht = group_ht.filter(group_ht.n_breaks == i)
        transcripts = temp_ht.aggregate(hl.agg.collect_as_set(temp_ht.transcript))
        if len(transcripts > 0):
            transcripts_per_break[i] = transcripts
        else:
            transcripts_per_break[i] = hl.missing(hl.tset(hl.tstr))
    return transcripts_per_break


def reformat_annotations_for_release(ht: hl.Table) -> hl.Table:
    """
    Reformat annotations in input HT for release.

    This reformatting is necessary to load data into the gnomAD browser.

    Desired schema:
    ---------------------------------------
    Row fields:
        'transcript_id': str
        'regions': array<struct {
            'start': int32
            'stop': int32
            'observed_missense': int32
            'expected_missense': float64
            'chisq': float64
        }>

    ----------------------------------------
    Key: ['transcript_id']
    ----------------------------------------

    :param hl.Table ht: Input Table.
    :return: Table with schema described above.
    :rtype: hl.Table
    """
    ht = ht.annotate(
        regions=hl.struct(
            start=ht.section_start,
            stop=ht.section_end,
            observed_missense=ht.section_obs,
            expected_missense=ht.section_exp,
            chisq=ht.section_chisq,
        )
    )
    # Group Table by transcript
    return ht.group_by(transcript_id=ht.transcript).aggregate(
        regions=hl.agg.collect(ht.regions)
    )


def finalize_multiple_breaks(
    ht: hl.Table,
    max_n_breaks: int,
    annotations: List[str] = FINAL_ANNOTATIONS,
) -> hl.Table:
    """
    Organize table of transcripts with multiple breaks.

    Get number of transcripts unique to each break number and drop any extra annotations.
    Also calculate each section's observed, expected, OE, and chi-square values.

    Assumes:
        - Table is annotated with set of transcripts per break (e.g., `break_1_transcripts`)

    :param hl.Table ht: Input Table. Example schema (truncated at two breaks for space reasons):
        ---------------------------------------
        Row fields:
            'locus': locus<grch37>
            'transcript': str
            'mu_snp': float64
            'observed': int32
            'total_exp': float64
            'total_mu': float64
            'total_obs': int64
            'max_chisq': float64
            'break_list': array<bool>
            'break_1_null': float64
            'break_1_alt': float64
            'break_1_chisq': float64
            'break_2_chisq': float64
            'break_2_max_chisq': float64
            'break_2_null': float64
            'break_2_alt': float64
        ----------------------------------------
        Key: ['locus', 'transcript']
        ----------------------------------------
    :param int max_n_breaks: Largest number of breakpoints in any transcript.
    :param List[str] annotations: List of annotations to keep from input Table.
        Default is FINAL_ANNOTATIONS.
    :return: Table annotated with transcript subsection values and breakpoint positions.
    :rtype: hl.Table
    """
    logger.info("Getting transcripts associated with each break number...")
    # Get number of transcripts UNIQUE to each break number
    # Transcript sets in globals currently are not unique
    # i.e., `break_1_transcripts` could contain transcripts that are also present in `break_2_transcripts`
    ht = ht.select_globals()
    transcripts_per_break = get_unique_transcripts_per_break(ht, max_n_breaks)

    # Print number of transcripts per break to output
    # This is used to create TSV input to `n_transcripts_per_break.R`
    for break_num in transcripts_per_break:
        logger.info(
            "Break number %i has %i transcripts",
            break_num,
            len(transcripts_per_break[break_num]),
        )

    logger.info("Selecting only relevant annotations from HT and checkpointing...")
    ht = ht.select(*annotations)
    ht = ht.checkpoint(f"{temp_path}/multiple_breaks.ht", overwrite=True)

    logger.info("Getting all breakpoint positions...")
    break_ht = get_all_breakpoint_pos(ht)
    break_ht = break_ht.checkpoint(
        f"{temp_path}/multiple_breaks_breakpoints.ht", overwrite=True
    )
    ht = ht.annotate(break_pos=break_ht[ht.transcript].break_pos)

    logger.info("Get transcript section annotations (obs, exp, OE, chisq)...")
    hts = []
    for i in range(1, max_n_breaks + 1):
        transcripts = hl.literal(transcripts_per_break[i])
        if hl.is_missing(transcripts):
            logger.info("Break number %i has no transcripts. Continuing...")

        # Filter HT to transcripts associated with this break only
        # and annotate section information
        temp_ht = ht.filter(transcripts.contains(ht.transcript))
        temp_ht = annotate_transcript_sections(temp_ht, i)

        # Add section oe
        # Do not cap section oe value here (this is for browser display)
        temp_ht = temp_ht.annotate(section_oe=temp_ht.section_obs / temp_ht.section_exp)
        temp_ht = temp_ht.checkpoint(
            f"gs://regional_missense_constraint/temp/break_{i}_sections.ht",
            overwrite=True,
        )
        hts.append(temp_ht)

    # Merge all HTs together
    ht = hts[0].union(*hts[1:])
    ht = ht.annotate_globals(**transcripts_per_break)
    ht = ht.checkpoint(multiple_breaks.path, overwrite=True)
    return ht


def finalize_simul_breaks(ht: hl.Table) -> hl.Table:
    """
    Process Table containing transcripts with two simultaneous breaks.

    Calculate each transcript subsection observed, expected, OE, and chi-square values.

    :param hl.Table ht: Input Table with simultaneous breaks results.
    :return: Table annotated with transcript subsection values.
    :rtype: hl.Table
    """
    logger.info("Get transcript section annotations (obs, exp, OE, chisq)...")
    ht = annotate_transcript_sections(ht, max_n_breaks=2)
    ht = ht.checkpoint(f"{temp_path}/simul_break_sections.ht", overwrite=True)
    return ht


def finalize_all_breaks_results(
    breaks_ht: hl.Table,
    simul_breaks_ht: hl.Table,
    max_n_breaks: int,
    annotations: List[str] = FINAL_ANNOTATIONS,
) -> None:
    """
    Finalize all breaks results.

    Finalize results for multiple breaks (first/additional breaks) and simultaneous breaks,
    then merge results into single Table.

    Also reformat Table to match desired release HT schema.

    :param hl.Table breaks_ht: Input Table with multiple break results.
    :param hl.Table simul_breaks_ht: Input Table with simultaneous breaks results.
    :param int max_n_breaks: Largest number of breakpoints in any transcript. Used only for multiple breaks results.
    :param List[str] annotations: List of annotations to keep from input Table.
        Default is FINAL_ANNOTATIONS.
    :return: None; writes Table to resource path.
    """
    logger.info("Finalizing multiple breaks results...")
    breaks_ht = finalize_multiple_breaks(breaks_ht, max_n_breaks, annotations)

    logger.info("Finalizing simultaneous breaks results...")
    simul_breaks_ht = finalize_simul_breaks(simul_breaks_ht)
    ht = breaks_ht.union(simul_breaks_ht, unify=False)
    ht = ht.checkpoint(f"{temp_path}/breaks.ht", overwrite=True)

    logger.info("Reformatting for browser release...")
    ht = reformat_annotations_for_release(ht)
    ht.write(rmc_browser.path, overwrite=True)


def check_loci_existence(ht1: hl.Table, ht2: hl.Table, annot_str: str) -> hl.Table:
    """
    Check if loci from `ht1` are present in `ht2`.

    Annotate `ht1` with int showing whether locus is present in `ht2`.
    Annotation will be "0" if locus isn't present in `ht2` and "1" if locus is present.

    :param hl.Table ht1: Table to be annotated.
    :param hl.Table ht2: Table to check for loci from `ht1`.
    :param str annot_str: Name of annotation to be added to `ht1` designating whether locus is present in `ht2`.
    :return: Annotated version of `ht1`.
    :rtype: hl.Table
    """
    return ht1.annotate(**{f"{annot_str}": hl.int(hl.is_defined(ht2[ht1.locus]))})


def constraint_flag_expr(
    obs_syn_expr: hl.expr.Int64Expression,
    obs_mis_expr: hl.expr.Int64Expression,
    obs_lof_expr: hl.expr.Int64Expression,
    exp_syn_expr: hl.expr.Float64Expression,
    exp_mis_expr: hl.expr.Float64Expression,
    exp_lof_expr: hl.expr.Float64Expression,
    raw_syn_z_expr: hl.expr.Float64Expression,
    raw_mis_z_expr: hl.expr.Float64Expression,
    raw_lof_z_expr: hl.expr.Float64Expression,
) -> hl.expr.StructExpression:
    """
    Return struct with constraint flags.

    Flags are designed to mark outlier transcripts, and explanation of constraint flags is in the gnomAD browser FAQ:
    https://gnomad.broadinstitute.org/faq#why-are-constraint-metrics-missing-for-this-gene-or-annotated-with-a-note

    :param hl.expr.Int64Expression obs_syn_expr: Expression containing number of observed synonymous variants in gnomAD.
    :param hl.expr.Int64Expression obs_mis_expr: Expression containing number of observed missense variants in gnomAD.
    :param hl.expr.Int64Expression obs_lof_expr: Expression containing number of observed loss-of-function (LoF) variants in gnomAD.
    :param hl.expr.Float64Expression exp_syn_expr: Expression containing number of expected synonymous variants.
    :param hl.expr.Float64Expression exp_mis_expr: Expression containing number of expected missense variants.
    :param hl.expr.Float64Expression exp_lof_expr: Expression containing number of expected LoF variants.
    :param hl.expr.Float64Expression raw_syn_z_expr: Expression containing number of Z score for synonymous variants.
    :param hl.expr.Float64Expression raw_mis_z_expr: Expression containing number of Z score for missense variants.
    :param hl.expr.Float64Expression raw_lof_z_expr: Expression containing number of Z score for LoF variants.
    :return: StructExpression containing constraint flags.
    :rtype: hl.expr.StructExpression
    """
    return hl.struct(
        no_variants=(
            hl.or_else(obs_syn_expr, 0)
            + hl.or_else(obs_mis_expr, 0)
            + hl.or_else(obs_lof_expr, 0)
        )
        == 0,
        no_exp_syn=exp_syn_expr == 0,
        no_exp_mis=exp_mis_expr == 0,
        no_exp_lof=exp_lof_expr == 0,
        syn_outlier=hl.abs(raw_syn_z_expr) > 5,
        mis_too_many=raw_mis_z_expr < -5,
        lof_too_many=raw_lof_z_expr < -5,
    )


def group_rmc_ht_by_section(overwrite: bool = False) -> hl.Table:
    """
    Group RMC results Table by transcript subsection and return interval and section missense o/e.

    .. note::
        - Function reads RMC results Table from resource path.
        - Assumes RMC HT is annotated with `locus`, `transcript`, `section`, `section_start_pos`,
        `section_end_pos`, and `section_oe`.
        - Assumes `transcript` is one of RMC HT's key fields.

    :param bool overwrite: Whether to overwrite temporary checkpointed Table if it exists.
        Default is False.
    :return: RMC results Table keyed by interval and annotated with transcript and section o/e.
    """
    rmc_ht = (
        rmc_results.ht()
        .select_globals()
        .select("section", "section_start_pos", "section_end_pos", "section_oe")
    )
    rmc_ht = rmc_ht.group_by("section").aggregate(
        transcript=hl.agg.take(rmc_ht.transcript, 1)[0],
        start_pos=hl.agg.take(rmc_ht.section_start_pos, 1)[0],
        end_pos=hl.agg.take(rmc_ht.section_end_pos, 1)[0],
        contig=hl.agg.take(rmc_ht.locus.contig, 1)[0],
        section_oe=hl.agg.take(rmc_ht.section_oe, 1)[0],
    )
    rmc_ht = rmc_ht.checkpoint(
        f"{temp_path}/rmc_group_by_section.ht",
        _read_if_exists=not overwrite,
        overwrite=overwrite,
    )
    rmc_ht = rmc_ht.transmute(
        interval=hl.parse_locus_interval(
            hl.format(
                "[%s:%s-%s]",
                rmc_ht.contig,
                rmc_ht.start_pos,
                rmc_ht.end_pos,
            )
        ),
    )
    return rmc_ht.key_by("interval").drop("section")
