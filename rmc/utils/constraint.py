import logging
from typing import Dict, List, Optional, Union

import hail as hl

from gnomad.resources.resource_utils import DataException

from gnomad.utils.reference_genome import get_reference_genome

from gnomad_lof.constraint_utils.generic import annotate_variant_types

from rmc.resources.basics import temp_path
from rmc.utils.generic import (
    get_coverage_correction_expr,
    get_exome_bases,
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
    transcript_expr: hl.expr.StringExpression, mu_expr: hl.expr.Float64Expression,
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


def calculate_exp_per_transcript(
    context_ht: hl.Table, locus_type: str, groupings: List[str] = GROUPINGS,
) -> hl.Table:
    """
    Return the total number of expected variants and aggregate mutation rate per transcript.

    .. note::
        - Assumes that context_ht is annotated with all of the fields in `groupings` and that the names match exactly.
        - Assumes that input table is filtered to autosomes/PAR only, X nonPAR only, or Y nonPAR only.
        - Expects that input table contains coverage and plateau models in its global annotations (`coverage_model`, `plateau_models`).
        - Expects that input table has multiple fields for mutation rate probabilities:
            `mu_snp` for mutation rate probability adjusted by coverage, and
            `raw_mu_snp` for raw mutation rate probability.

    :param hl.Table context_ht: Context Table.
    :param str locus_type: Locus type of input table. One of "X", "Y", or "autosomes".
        NOTE: will treat any input other than "X" or "Y" as autosomes.
    :param List[str] groupings: List of Table fields used to group Table to adjust mutation rate.
        Table must be annotated with these fields. Default is GROUPINGS.
    :return: Table grouped by transcript with expected counts per search field.
    :rtype: hl.Table
    """
    logger.info(f"Grouping by {groupings}...")
    group_ht = context_ht.group_by(*groupings).aggregate(
        mu_agg=hl.agg.sum(context_ht.raw_mu_snp)
    )

    logger.info("Adding CpG annotations...")
    group_ht = annotate_variant_types(group_ht)

    logger.info("Adjusting aggregated mutation rate with plateau model...")
    if locus_type == "X":
        model = group_ht.plateau_x_models["total"][group_ht.cpg]
    elif locus_type == "Y":
        model = group_ht.plateau_y_models["total"][group_ht.cpg]
    else:
        model = group_ht.plateau_models["total"][group_ht.cpg]

    group_ht = group_ht.annotate(mu_adj=group_ht.mu_agg * model[1] + model[0])

    logger.info(
        "Adjusting aggregated mutation rate with coverage correction to get expected counts..."
    )
    group_ht = group_ht.annotate(
        coverage_correction=get_coverage_correction_expr(
            group_ht.exome_coverage, group_ht.coverage_model
        ),
    )
    group_ht = group_ht.annotate(
        _exp=group_ht.mu_adj * group_ht.coverage_correction,
        mu=group_ht.mu_agg * group_ht.coverage_correction,
    )

    logger.info("Getting expected counts per transcript and returning...")
    return group_ht.group_by("transcript").aggregate(
        expected=hl.agg.sum(group_ht._exp), mu_agg=hl.agg.sum(group_ht.mu)
    )


def get_obs_exp_expr(
    cond_expr: hl.expr.BooleanExpression,
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
) -> hl.expr.Float64Expression:
    """
    Return observed/expected annotation based on inputs.

    Cap observed/expected value at 1.

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
    transcript_expr: hl.expr.StringExpression, observed_expr: hl.expr.Int64Expression,
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

    .. note::
        This function is designed to run on one transcript at a time.

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
            transcript_expr=ht[transcript_str], observed_expr=ht[obs_str],
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
            ht.cond_expr, ht.cumulative_obs[ht[transcript_str]], ht.cumulative_exp,
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


def get_section_expr(dpois_expr: hl.expr.ArrayExpression,) -> hl.expr.Float64Expression:
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
    simul_break: bool = False,
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

    Returns HT filtered to lines with maximum chisq if chisq >= max_value, otherwise returns None.

    :param hl.Table ht: Input context Table.
    :param hl.expr.StringExpression search_field: Field of table to search. Value should be either 'transcript' or 'section'.
    :param bool simul_break: Whether this function is searching for simultaneous breaks. Default is False.
    :param float chisq_threshold: Chi-square significance threshold.
        Value should be 10.8 (single break) and 13.8 (two breaks) (values from ExAC RMC code).
        Default is 10.8.
    :return: Table annotated with whether position is a breakpoint.
    :rtype: hl.Table
    """
    logger.info(
        "Creating section null (no regional variability in missense depletion)\
        and alt (evidence of domains of missense constraint) expressions..."
    )

    # Split transcript into three sections if searching for simultaneous breaks
    if simul_break:
        ht = ht.annotate(
            section_nulls=[
                # Get null expression for section of transcript pre-window
                get_dpois_expr(
                    cond_expr=True,
                    section_oe_expr=ht.overall_oe,
                    obs_expr=ht.pre_obs,
                    exp_expr=ht.pre_exp,
                ),
                # Get null expression for window of constraint
                get_dpois_expr(
                    cond_expr=True,
                    section_oe_expr=ht.overall_oe,
                    obs_expr=ht.window_obs,
                    exp_expr=ht.window_exp,
                ),
                # Get null expression for section of transcript post-window
                get_dpois_expr(
                    cond_expr=True,
                    section_oe_expr=ht.overall_oe,
                    obs_expr=ht.post_obs,
                    exp_expr=ht.post_exp,
                ),
            ]
        )
        ht = ht.annotate(
            section_alts=[
                # Get alt expression for section of transcript pre-window
                get_dpois_expr(
                    cond_expr=True,
                    section_oe_expr=ht.pre_oe,
                    obs_expr=ht.pre_obs,
                    exp_expr=ht.pre_exp,
                ),
                # Get alt expression for window of constraint
                get_dpois_expr(
                    cond_expr=True,
                    section_oe_expr=ht.window_oe,
                    obs_expr=ht.window_obs,
                    exp_expr=ht.window_exp,
                ),
                # Get alt expression for section of transcript post-window
                get_dpois_expr(
                    cond_expr=True,
                    section_oe_expr=ht.post_oe,
                    obs_expr=ht.post_obs,
                    exp_expr=ht.next_values.reverse_exp,
                ),
            ]
        )

    # Otherwise, split transcript only into two sections (when searching for first/additional breaks)
    else:
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

    return search_for_break(ht, "transcript", False, chisq_threshold)


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
        obs=hl.agg.sum(ht[obs_str]), mu=hl.agg.sum(ht[mu_str]),
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


def process_sections(ht: hl.Table, chisq_threshold: float):
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
    :return: Table annotated with whether position is a breakpoint.
    :rtype: hl.Table
    """
    ht = get_subsection_exprs(ht)

    # TODO: Move this calculation to release HT generation step
    logger.info("Getting section chi-squared values...")
    ht = ht.annotate(
        section_chisq=calculate_section_chisq(
            obs_expr=ht.section_obs, exp_expr=ht.section_exp,
        )
    )

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
        "Annotating post breakpoint HT with reverse observed and expected counts..."
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
        simul_break=False,
        chisq_threshold=chisq_threshold,
    )
    # Adjust is_break annotation in pre_ht
    # to prevent this function from continually finding previous significant breaks
    pre_ht = pre_ht.annotate(
        is_break=hl.if_else(pre_ht.break_list.any(lambda x: x), False, pre_ht.is_break)
    )
    post_ht = search_for_break(
        post_ht,
        search_field="transcript",
        simul_break=False,
        chisq_threshold=chisq_threshold,
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
    return process_sections(ht, chisq_threshold)


def get_avg_bases_between_mis(ht: hl.Table) -> int:
    """
    Return average number of bases between observed missense variation.

    For example, if the total number of bases is 30, and the total number of missense variants is 10,
    this function will return 3.

    This function is used to determine the minimum size window to check for significant missense depletion
    when searching for two simultaneous breaks.

    .. note::
        Assumes input Table has been filtered to missense variants in canonical protein-coding transcripts only.

    :param hl.Table ht: Input gnomAD exomes Table.
    :return: Average number of bases between observed missense variants, rounded to the nearest integer,
    :rtype: int
    """
    logger.info("Getting total number of bases in the exome from full context HT...")
    total_bases = get_exome_bases(build=get_reference_genome(ht.locus).name)

    logger.info("Getting total number of missense variants in gnomAD...")
    total_variants = ht.count()
    logger.info(f"Total number of bases in the exome: {total_bases}")
    logger.info(f"Total number of missense variants in gnomAD exomes: {total_variants}")

    logger.info("Getting average bases between missense variants and returning...")
    return round(total_bases / total_variants)


def search_for_two_breaks(
    ht: hl.Table,
    break_size: int,
    max_break_size: Optional[int] = None,
    chisq_threshold: float = 13.8,
) -> hl.Table:
    """
    Search for evidence of constraint within a set window size/number of base pairs.

    Function is designed to search in transcripts that didn't have one single significant break.

    Assumes that:
        - Input Table has a field named 'transcript'.
        - Input Table is annotated with a field named 'max_window_size' if `max_break_size` argument is not None.

    :param hl.Table ht: Input Table.
    :param int break_size: Number of bases to search for constraint (window size for simultaneous breaks).
    :param int max_break_size: Maximum window size to search for constraint. Only used when expanding
        boundaries of simultaneous breaks windows. Default is None.
    :param float chisq_threshold: Chi-square significance threshold.
        Value should be 10.8 (single break) and 13.8 (two breaks) (values from ExAC RMC code).
    :return: Table annotated with is_break at the *end* position of a simultaneous break window.
    :rtype: hl.Table
    """
    if max_break_size:
        # Double check input HT has max_window_size annotation if max_break_size is set
        if "max_window_size" not in ht.row:
            raise DataException(
                "HT must be annotated with max_window_size if max_break_size arg is set!"
            )

    logger.info(
        f"Annotating each position with end position for window \
        if position + {break_size - 1} bases is less than or equal to the transcript stop pos..."
    )
    ht = ht.annotate(
        window_end=hl.or_missing(
            # Window start is defined only if position + (break_size - 1) is less than or
            # equal to the transcript end position
            # e.g., if position is 8, break size is 3, and end position is 12, then window_end is 8 + (3-1) = 10
            (ht.locus.position + (break_size - 1) <= ht.end_pos),
            ht.locus.position + (break_size - 1),
        )
    )

    def _get_pos_per_transcript(
        ht: hl.Table, transcripts: hl.expr.SetExpression
    ) -> hl.Table:
        """
        Filter input HT to input transcripts and gather all positions per transcript.

        Input HT is HT containing all transcripts without a single significant break and stripped of all annotations
        aside from locus and transcript.

        :param hl.Table ht: Input HT.
        :param hl.expr.SetExpression transcripts: SetExpression of transcripts to extract.
        :return: Table grouped by transcript with list of positions per transcript.
        :rtype: hl.Table
        """
        logger.info("Gathering all positions in each transcript...")
        pos_ht = ht.filter(transcripts.contains(ht.transcript))
        pos_ht = pos_ht.group_by(transcript=pos_ht.transcript).aggregate(
            positions=hl.sorted(hl.agg.collect(pos_ht.locus.position))
        )
        pos_ht = pos_ht.checkpoint(f"{temp_path}/pos_per_transcript.ht", overwrite=True)
        return pos_ht

    def _get_post_window_pos(
        ht: hl.Table, pos_ht: hl.Table, has_end: bool = False
    ) -> hl.Table:
        """
        Get first position of transcript outside of window of constraint.

        Run `hl.binary_search` to find index of first post-window position.

        :param hl.Table ht: Input Table.
        :param hl.Table pos_ht: Input GroupedTable grouped by transcript with list of all positions per transcript.
        :param bool has_end: Whether the input HT has window ends defined in the HT. Default is False.
        :return: Table annotated with post-window position.
        :rtype: hl.Table
        """
        logger.info("Annotating the positions per transcript back onto the input HT...")
        # This is to run `hl.binary_search` using the window end position as input
        ht = ht.key_by("transcript").join(pos_ht)
        ht = ht.transmute(pos_per_transcript=ht.positions)

        logger.info(
            "Running hl.binary_search to find the index of the first position post-window..."
        )
        # hl.binary search will return the index of the position that is closest to the search element (window_end)
        # If window_end exists in the HT, hl.binary_search will return the index of the window_end
        # If window_end does not exist in the HT, hl.binary_search will return the index of the next largest position
        # e.g., pos_per_transcript = [1, 4, 8]
        # hl.binary_search(pos_per_transcript, 4) will return 1
        # hl.binary_search(pos_per_transcript, 5) will return 2
        ht = ht.annotate(
            post_window_index=hl.binary_search(
                ht.pos_per_transcript, ht.window_end.position
            ),
            n_pos_per_transcript=hl.len(ht.pos_per_transcript),
        )
        if has_end:
            return ht.annotate(
                post_window_pos=hl.case()
                .when(
                    ht.post_window_index + 1 < ht.n_pos_per_transcript - 1,
                    ht.pos_per_transcript[ht.post_window_index + 1],
                )
                # hl.binary_search will only return an index larger than the length of a list
                # if that position is larger than the largest element in the list
                # e.g., pos_per_transcript = [1, 4, 8]
                # hl.binary_search(pos_per_transcript, 10) will return 3
                # post-window positions should never be larger than the largest position of a transcript
                # (so this should never happen)
                .when(
                    ht.post_window_index + 1 == ht.n_pos_per_transcript - 1, ht.end_pos
                )
                .or_missing()
            )
        return ht.annotate(
            post_window_pos=hl.case()
            .when(
                ht.post_window_index < ht.n_pos_per_transcript - 1,
                ht.pos_per_transcript[ht.post_window_index],
            )
            .when(ht.post_window_index == ht.n_pos_per_transcript - 1, ht.end_pos)
            .or_missing(),
        )

    logger.info(
        "Checking how many window end positions are actually defined in the input HT..."
    )
    # The input HT contains every possible missense variant
    # and will be missing positions if no missense variant was possible at that position
    # Duplicate HT to check window ends
    window_ht = ht.select().select_globals()
    window_ht = window_ht.key_by("locus")

    # Keep version of HT with all relevant annotations and strip HT of all annotations
    annotation_ht = ht.select(
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
    )
    ht = ht.select(
        "window_end", "start_pos", "end_pos", "transcript_size"
    ).select_globals()
    ht = ht.annotate(window_end=hl.locus(ht.locus.contig, ht.window_end))

    logger.info(
        "Checkpointing HT with sites that have their window ends defined in the HT..."
    )
    end_ht = ht.key_by("window_end").join(window_ht, how="inner")
    end_ht = end_ht.checkpoint(
        f"{temp_path}/simul_break_prep_end_def.ht", overwrite=True
    )
    end_transcripts = end_ht.aggregate(
        hl.agg.collect_as_set(end_ht.transcript), _localize=False
    )
    logger.info(f"Sites with defined window ends: {end_ht.count()}")

    logger.info(
        "Adding post-window position for sites that have their window ends defined..."
    )
    # NOTE: The annotation post_window_pos must actually be outside the window of constraint
    # Otherwise, `search_for_break` will undercount the obs/exp variant counts within the window of constraint
    pos_ht = _get_pos_per_transcript(window_ht, end_transcripts)
    end_ht = _get_post_window_pos(end_ht, pos_ht, has_end=True)
    end_ht = end_ht.checkpoint(
        f"{temp_path}/simul_break_prep_end_def_ready.ht", overwrite=True
    )

    logger.info(
        "Working on sites that don't have their window ends defined in the HT..."
    )
    no_end_ht = ht.key_by("window_end").anti_join(end_ht.key_by("window_end"))
    no_end_ht.write(f"{temp_path}/simul_break_prep_no_end.ht", overwrite=True)
    no_end_ht = hl.read_table(
        f"{temp_path}/simul_break_prep_no_end.ht", _n_partitions=1000
    )
    no_end_transcripts = no_end_ht.aggregate(
        hl.agg.collect_as_set(no_end_ht.transcript), _localize=False
    )
    logger.info(f"Sites without defined window ends: {no_end_ht.count()}")

    logger.info(
        "Adding post-window position for sites that don't have defined window ends..."
    )
    pos_ht = _get_pos_per_transcript(window_ht, no_end_transcripts)
    no_end_ht = _get_post_window_pos(no_end_ht, pos_ht)
    no_end_ht.write(f"{temp_path}/simul_break_prep_no_end_ready.ht", overwrite=True)
    no_end_ht = hl.read_table(
        f"{temp_path}/simul_break_prep_no_end_ready.ht", _n_partitions=5000
    )

    logger.info("Joining no end HT with end HT...")
    end_ht = end_ht.drop(
        "pos_per_transcript", "post_window_index", "n_pos_per_transcript"
    )
    no_end_ht = no_end_ht.drop(
        "pos_per_transcript", "post_window_index", "n_pos_per_transcript"
    )
    ht = end_ht.key_by("locus", "transcript").join(
        no_end_ht.key_by("locus", "transcript"), how="outer"
    )
    ht = ht.annotate(
        post_window_pos=hl.coalesce(ht.post_window_pos, ht.post_window_pos_1),
        **annotation_ht[ht.key],
    )
    ht = ht.checkpoint(f"{temp_path}/simul_break_prep.ht", overwrite=True)
    logger.info(f"HT count: {ht.count()}")

    # Check if post window pos is ever smaller than window_end
    # Convert window end from locus into position
    ht = ht.annotate(window_end=ht.window_end.position)
    check_end = ht.aggregate(hl.agg.count_where(ht.window_end > ht.post_window_pos))
    if check_end > 0:
        ht = ht.filter(ht.window_end > ht.post_window_pos)
        ht.show()
        raise DataException(
            f"Position closest to window end for HT is smaller than window end position in {check_end} cases!"
        )

    # Create new HT with obs, exp, and OE values for post-window positions
    # This is to get the values for the section of the transcript after the window
    next_ht = ht.select("post_window_pos")
    # Add new_locus annotation to pull obs, exp, OE values for post window position
    next_ht = next_ht.annotate(
        new_locus=hl.locus(next_ht.locus.contig, next_ht.post_window_pos)
    )
    next_ht = next_ht.annotate(
        next_values=hl.struct(
            cum_obs=ht[next_ht.new_locus, next_ht.transcript].cumulative_obs[
                next_ht.transcript
            ],
            obs=ht[next_ht.new_locus, next_ht.transcript].observed,
            exp=ht[next_ht.new_locus, next_ht.transcript].cumulative_exp,
            mu_snp=ht[next_ht.new_locus, next_ht.transcript].mu_snp,
            oe=ht[next_ht.new_locus, next_ht.transcript].forward_oe,
            reverse_obs=ht[next_ht.new_locus, next_ht.transcript].reverse.obs,
            reverse_exp=ht[next_ht.new_locus, next_ht.transcript].reverse.exp,
        )
    )
    ht = ht.annotate(next_values=next_ht[ht.key].next_values)

    # Annotate expected variant count at current position
    ht = ht.annotate(
        # Translate mu_snp at site to expected at window start site
        # can't use translate_mu_to_exp_expr because that is expecting
        # cumulative mu, and we want to use the value for this site only
        exp_at_site=(ht.mu_snp / ht.total_mu)
        * ht.total_exp,
    )

    # Annotate observed and expected counts for section of transcript pre-window
    # These values are the current cumulative values minus the value for the site
    ht = ht.annotate(
        # Annotate observed count for section of transcript pre-window
        # = current cumulative obs minus the obs at position
        # Use hl.max to keep this value positive
        pre_obs=hl.max(ht.cumulative_obs[ht.transcript] - ht.observed, 0),
        pre_exp=hl.max(ht.cumulative_exp - ht.exp_at_site, 0),
    )
    # Make sure prev exp value isn't 0
    # Otherwise, the chisq calculation will return NaN
    ht = ht.annotate(pre_exp=hl.if_else(ht.pre_exp > 0, ht.pre_exp, 1e-09))

    # Annotate OE value for section of transcript pre-window
    ht = ht.annotate(
        pre_oe=get_obs_exp_expr(
            cond_expr=True, obs_expr=ht.pre_obs, exp_expr=ht.pre_exp,
        )
    )

    # Annotate obsered and expected counts for window of constraint
    # These values are the cumulative values at the first position post-window minus
    # the cumulative values at the current window
    # Also minus the value of the observed or expected variants at the first position post-window
    # This subtraction is for both obs and exp in case the first position post window has an observed variant
    # (otherwise you'd overcount the number of obs within the window)
    ht = ht.annotate(
        window_obs=(ht.next_values.cum_obs - ht.cumulative_obs[ht.transcript])
        - ht.next_values.obs,
        window_exp=(ht.next_values.exp - ht.cumulative_exp)
        - ((ht.next_values.mu_snp / ht.total_mu) * ht.total_exp),
    )

    # Annotate OE value for section of transcript within window
    ht = ht.annotate(
        window_oe=get_obs_exp_expr(
            cond_expr=True, obs_expr=ht.window_obs, exp_expr=ht.window_exp,
        )
    )

    # Annotate observed and expected counts for section of transcript post-window
    # These values are the reverse values at the first position post-window plus
    # the observed and expected values at that first position post-window
    # (need to add here, otherwise would undercount the number post-window)
    ht = ht.annotate(
        post_obs=ht.next_values.reverse_obs + ht.next_values.obs,
        post_exp=ht.next_values.reverse_exp
        + ((ht.next_values.mu_snp / ht.total_mu) * ht.total_exp),
    )

    # Annotate OE value for section of transcript post-window
    ht = ht.annotate(
        post_oe=get_obs_exp_expr(
            cond_expr=True, obs_expr=ht.post_obs, exp_expr=ht.post_exp,
        )
    )

    logger.info("Searching for two breaks...")
    return search_for_break(
        ht, "transcript", simul_break=True, chisq_threshold=chisq_threshold
    )


def expand_two_break_window(
    ht: hl.Table, transcript_percentage: float, chisq_threshold: float = 13.8,
) -> hl.Table:
    """
    Search for larger windows of constraint in transcripts with simultaneous breaks.

    This function searches for breaks for break sizes of `min_window_size` to `max_window_size`,
    where `min_window_size` is the number of base pairs needed, on average, to see 10 missense variants, and
    `max_window_size` is the number of base pairs equal to `transcript_percentage` of the largest transcript size.
    For example, for gnomAD v2.1, `min_window_size` is 100bp, and `max_window_size` is 90% of the largest transcript.

    :param hl.Table ht: Input Table filtered to contain only transcripts with simultaneous breaks.
    :param float transcript_percentage: Maximum percentage of the transcript that can be included within a window of constraint.
    :param float chisq_threshold:  Chi-square significance threshold. Default is 13.8.
    :return: Table with largest simultaneous break window size annotated per transcript.
    :rtype: hl.Table
    """
    # Get smallest possible window size
    # NOTE: Hardcoding 10 here as smallest number of obs mis, not sure if this will change in the future
    min_window_size = hl.eval(ht.obs_mis_10_window_size)
    logger.info("Minimum window size: %i", min_window_size)

    logger.info("Annotating each transcript with max window size...")
    ht = ht.annotate(max_window_size=ht.transcript_size * transcript_percentage)
    max_window_size = round(ht.aggregate(hl.agg.max(ht.max_window_size)))
    logger.info("Maximum window size: %i", max_window_size)

    logger.info("Annotating each transcript with current window size...")
    # Also drop any unnecessary annotations
    ht = ht.select(
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
        "max_window_size",
        break_sizes=(
            hl.case()
            .when(ht.obs_mis_50.contains(ht.transcript), [ht.obs_mis_50_window_size])
            .when(ht.obs_mis_20.contains(ht.transcript), [ht.obs_mis_20_window_size])
            .default([ht.obs_mis_10_window_size])
        ),
        break_chisqs=[ht.max_chisq],
        window_ends=[ht.window_end],
        post_window_pos=[ht.post_window_pos],
    )
    ht.describe()
    ht = ht.checkpoint(f"{temp_path}/simul_break_expand_ready.ht", overwrite=True)

    logger.info("Expanding window sizes...")
    window_size = min_window_size
    while window_size <= max_window_size:
        break_ht = search_for_two_breaks(
            ht=ht,
            break_size=window_size,
            max_break_size=max_window_size,
            chisq_threshold=chisq_threshold,
        )

        # Re-add annotations that are dropped in search_for_two_breaks function
        annot_ht = ht.select(
            "break_sizes", "break_chisqs", "window_ends", "post_window_pos"
        )
        break_ht = break_ht.join(annot_ht, how="left")
        break_ht = break_ht.annotate(
            break_sizes=break_ht.break_sizes.append(window_size),
            break_chisqs=break_ht.break_chisqs.append(break_ht.max_chisq),
            window_ends=break_ht.window_ends.append(break_ht.window_end),
        )
        # This method will checkpoint a LOT of temporary tables...not sure if there is a better way
        break_ht = break_ht.checkpoint(
            f"{temp_path}/simul_break_{window_size}_window.ht", overwrite=True
        )
        break_ht.describe()
        ht = ht.annotate(**break_ht[ht.key])
        window_size += 1

    logger.info("Getting best window sizes for each transcript and returning...")
    # Extract max chi square value
    ht = ht.annotate(max_chisq=hl.sorted(ht.break_chisqs, reverse=True)[0])
    # Get index associated with max chi square
    ht = ht.annotate(chisq_index=ht.break_chisqs.index(ht.max_chisq))
    return ht.annotate(best_window_size=ht.break_sizes[ht.chisq_index]).drop(
        "chisq_index"
    )


def calculate_section_chisq(
    obs_expr: hl.expr.Int64Expression, exp_expr: hl.expr.Float64Expression,
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


def fix_xg(
    context_ht: hl.Table,
    exome_ht: hl.Table,
    xg_transcript: str = "ENST00000419513",
    groupings: List[str] = GROUPINGS,
) -> hl.Table:
    """
    Fix observed and expected counts for XG (gene that spans PAR and non-PAR regions on chrX).

    Expects that context HT is annotated with all of the fields in `groupings`.

    :param hl.Table context_ht: Context Table.
    :param hl.Table exome_ht: Table containing variants from gnomAD exomes.
    :param str xg_transcript: XG transcript string. Default is 'ENST00000419513'.
    :param List[str] groupings: List of Table fields used to group Table to adjust mutation rate.
        Default is GROUPINGS.
    :return: Table filtered to XG and annotated with RMC annotations (forward scans, total obs/exp/mu, overall OE).
    :rtype: hl.Table
    """

    def _fix_xg_exp(
        xg: hl.Table, groupings: List[str] = groupings,
    ) -> hl.expr.StructExpression:
        """
        Fix total expected and total mu counts for XG.

        :param hl.Table xg: Context Table filtered to XG.
        :param List[str] groupings: List of Table fields used to group Table to adjust mutation rate.
            Default is GROUPINGS.
        :return: StructExpression with total mu and total expected values.
        :rtype: hl.expr.StructExpression
        """
        xg = xg.annotate(par=xg.locus.in_x_par())
        par = xg.filter(xg.par)
        nonpar = xg.filter(~xg.par)
        par = calculate_exp_per_transcript(
            par, locus_type="autosomes", groupings=groupings
        )
        nonpar = calculate_exp_per_transcript(
            nonpar, locus_type="X", groupings=groupings
        )
        exp_ht = par.union(nonpar)
        return exp_ht.aggregate(
            hl.struct(
                total_exp=hl.agg.sum(exp_ht.expected),
                total_mu=hl.agg.sum(exp_ht.mu_agg),
            )
        )

    def _fix_xg_obs(xg: hl.Table, exome_ht: hl.Table) -> hl.Table:
        """
        Fix total observed counts for XG.

        :param hl.Table xg: Context Table filtered to XG.
        :param hl.Table exome_ht: Table containing variants from gnomAD exomes.
        :return: Context Table filtered to XG, annotated with total observed and observed per position values.
        :rtype: hl.Table
        """
        xg = xg.annotate(_obs=exome_ht.index(xg.key))
        xg = xg.transmute(observed=hl.int(hl.is_defined(xg._obs)))
        return xg.annotate(total_obs=obs_ht[xg.transcript].observed)

    logger.info(f"Filtering context HT to XG (transcript: {xg_transcript})...")
    xg = context_ht.filter(context_ht.transcript == xg_transcript)

    logger.info("Fixing expected counts for XG...")
    exp_struct = _fix_xg_exp(xg, groupings)

    logger.info("Fixing observed counts for XG...")
    exome_ht = exome_ht.filter(
        exome_ht.transcript_consequences.transcript_id == xg_transcript
    )
    obs_ht = calculate_observed(exome_ht)
    xg = _fix_xg_obs(xg, exome_ht)

    logger.info("Collecting by key...")
    # Context HT is keyed by locus and allele, which means there is one row for every possible missense variant
    # This means that any locus could be present up to three times (once for each possible missense)
    # Collect by key here to ensure all loci are unique
    xg = xg.key_by("locus", "transcript").collect_by_key()
    xg = xg.annotate(
        # Collect the mutation rate probabilities at each locus
        mu_snp=hl.sum(xg.values.mu_snp),
        # Collect the observed counts for each locus
        # (this includes counts for each possible missense at the locus)
        observed=hl.sum(xg.values.observed),
        # Take just the first coverage value, since the locus should have the same coverage across the possible variants
        coverage=xg.values.exome_coverage[0],
    )

    logger.info(
        "Annotating overall observed/expected value (capped at 1) and forward scans..."
    )
    xg = xg.annotate(
        total_exp=exp_struct.total_exp,
        total_mu=exp_struct.total_mu,
        total_obs=obs_ht[xg.transcript].observed,
    )
    xg.describe()
    xg = xg.annotate(overall_oe=hl.min(xg.total_obs / xg.total_exp, 1))
    xg = get_fwd_exprs(
        ht=xg,
        transcript_str="transcript",
        obs_str="observed",
        mu_str="mu_snp",
        total_mu_str="total_mu",
        total_exp_str="total_exp",
    )
    xg = xg.checkpoint(f"{temp_path}/XG.ht", overwrite=True)
    return xg
