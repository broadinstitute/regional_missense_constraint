import logging
from typing import Dict, List, Union

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
    Groups input Table by transcript, filters based on `keep_criteria`,
    and aggregates observed variants count per transcript.

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
    Returns annotation with the cumulative mutation rate probability, shifted by one.

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
    Adjusts the scan with the cumulative number mutation rate probability.

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
    Translates cumulative mutation rate probability into cumulative expected count per base.

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
    Returns the total number of expected variants and aggregate mutation rate per transcript.

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
    Returns observed/expected annotation based on inputs.

    Caps observed/expected value at 1.

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
    Returns annotation with the cumulative number of observed variants, shifted by one.

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
    Adjusts the scan with the cumulative number of observed variants.

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
    Returns the "reverse" section observed and expected variant counts.

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
    Annotates input Table with the forward section cumulative observed, expected, and observed/expected values.

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
    Calls `get_reverse_obs_exp_expr` and `get_obs_exp_expr` to add the reverse section cumulative observed, expected, and observed/expected values.

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
    :param hl.expr.Float64Expression section_oe_expr: Expression of section observed/expected value.
    :param Union[Dict[hl.expr.StringExpression, hl.expr.Int64Expression], hl.expr.Int64Expression] obs_expr: Expression containing observed variants count.
    :param Union[Dict[hl.expr.StringExpression, hl.expr.Float64Expression], hl.expr.Float64Expression] exp_expr: Expression containing expected variants count.
    :return: Struct containing forward or reverse null and alt values (either when searching for first or second break).
    :rtype: hl.expr.StructExpression
    """
    return hl.or_missing(cond_expr, hl.dpois(obs_expr, exp_expr * section_oe_expr))


def get_section_expr(dpois_expr: hl.expr.ArrayExpression,) -> hl.expr.Float64Expression:
    """
    Builds null or alt model by multiplying all section null or alt distributions.

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
    prev: bool = False,
    chisq_threshold: float = 10.8,
) -> hl.Table:
    """
    Searches for breakpoints in a transcript or within a transcript subsection.

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
        If prev is True:
            - pre_window_pos and prev_values
        Otherwise:
            - post_window_pos and next_values
    Also expects:
        - multiallelic variants in input HT have been split.
        - Input HT is autosomes/PAR only, X non-PAR only, or Y non-PAR only.

    Returns HT filtered to lines with maximum chisq if chisq >= max_value, otherwise returns None.

    :param hl.Table ht: Input context Table.
    :param hl.expr.StringExpression search_field: Field of table to search. Value should be either 'transcript' or 'section'. 
    :param bool simul_break: Whether this function is searching for simultaneous breaks. Default is False.
    :param bool prev: Whether simultaneous breaks were calculated with code that uses prev_nonnull. Default if False.
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
        if prev:
            pre_window_expr = hl.format("%s_%s", ht.transcript, ht.pre_window_pos)
            ht = ht.annotate(
                section_null=[
                    # Get null expression for section of transcript pre-window
                    get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=ht.overall_oe,
                        obs_expr=ht.prev_values[pre_window_expr].prev_obs,
                        exp_expr=ht.prev_values[pre_window_expr].prev_exp,
                    ),
                    # Get null expression for window of constraint
                    get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=ht.overall_oe,
                        obs_expr=ht.cumulative_obs[ht.transcript],
                        exp_expr=ht.cumulative_exp,
                    ),
                    # Get null expression for section of transcript post-window
                    get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=ht.overall_oe,
                        obs_expr=ht.reverse.obs,
                        exp_expr=ht.reverse.exp,
                    ),
                ]
            )
            ht = ht.annotate(
                section_alt=[
                    # Get alt expression for section of transcript pre-window
                    get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=ht.prev_values[pre_window_expr].prev_oe,
                        obs_expr=ht.prev_values[pre_window_expr].prev_obs,
                        exp_expr=ht.prev_values[pre_window_expr].prev_exp,
                    ),
                    # Get alt expression for window of constraint
                    get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=get_obs_exp_expr(
                            cond_expr=True,
                            obs_expr=ht.cumulative_obs[ht.transcript],
                            exp_expr=ht.cumulative_exp,
                        ),
                        obs_expr=ht.cumulative_obs[ht.transcript],
                        exp_expr=ht.cumulative_exp,
                    ),
                    # Get alt expression for section of transcript post-window
                    get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=get_obs_exp_expr(
                            cond_expr=True,
                            obs_expr=ht.reverse.obs,
                            exp_expr=ht.reverse.exp,
                        ),
                        obs_expr=ht.reverse.obs,
                        exp_expr=ht.reverse.exp,
                    ),
                ]
            )
        else:
            # Annotate expected variant count at current position
            ht = ht.annotate(
                # Translate mu_snp at site to expected at window start site
                # can't use translate_mu_to_exp_expr because that is expecting
                # cumulative mu, and we want to use the value for this site only
                exp_at_site=(ht.mu_snp / ht.total_mu)
                * ht.total_exp,
            )
            # Annotate observed and expected counts for section of transcript pre-window
            ht = ht.annotate(
                # Annotate observed count for section of transcript pre-window
                # = current cumulative obs minus the obs at position
                # Use hl.max to keep this value positive
                prev_obs=hl.max(ht.cumulative_obs[ht.transcript] - ht.observed, 0),
                prev_exp=hl.max(ht.cumulative_exp - ht.exp_at_site, 0),
            )
            # Make sure prev exp value isn't 0
            # Otherwise, the chisq calculation will return NaN
            ht = ht.annotate(prev_exp=hl.if_else(ht.prev_exp > 0, ht.prev_exp, 1e-09))
            # Annotate OE value for section of transcript pre-window
            ht = ht.annotate(
                pre_oe=get_obs_exp_expr(
                    cond_expr=True, obs_expr=ht.prev_obs, exp_expr=ht.prev_exp,
                )
            )

            ht = ht.annotate(
                section_nulls=[
                    # Get null expression for section of transcript pre-window
                    # These values are the current cumulative values minus the value for the site
                    get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=ht.overall_oe,
                        obs_expr=ht.prev_obs,
                        exp_expr=ht.prev_exp,
                    ),
                    # Get null expression for window of constraint
                    # These values are the cumulative values at the first position post-window minus
                    # the cumulative values at the current window
                    # Also minus the value of the observed or expected variants at the first position post-window
                    # This subtraction is for both obs and exp in case the first position post window has an observed variant
                    # (otherwise you'd overcount the number of obs within the window)
                    get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=ht.overall_oe,
                        obs_expr=(
                            ht.next_values.cum_obs - ht.cumulative_obs[ht.transcript]
                        )
                        - ht.next_values.obs,
                        exp_expr=(ht.next_values.exp - ht.cumulative_exp)
                        - ((ht.next_values.mu_snp / ht.total_mu) * ht.total_exp),
                    ),
                    # Get null expression for section of transcript post-window
                    # These values are the reverse values at the first position post-window plus
                    # the observed and expected values at that first position post-window
                    # (need to add here, otherwise would undercount the number post-window)
                    get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=ht.overall_oe,
                        obs_expr=ht.next_values.reverse_obs + ht.next_values.obs,
                        exp_expr=ht.next_values.reverse_exp
                        + ((ht.next_values.mu_snp / ht.total_mu) * ht.total_exp),
                    ),
                ]
            )
            ht = ht.annotate(
                section_alts=[
                    # Get alt expression for section of transcript pre-window
                    get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=ht.pre_oe,
                        obs_expr=ht.prev_obs,
                        exp_expr=ht.prev_exp,
                    ),
                    # Get alt expression for window of constraint
                    get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=get_obs_exp_expr(
                            cond_expr=True,
                            obs_expr=(
                                ht.next_values.cum_obs
                                - ht.cumulative_obs[ht.transcript]
                            )
                            - ht.next_values.obs,
                            exp_expr=(ht.next_values.exp - ht.cumulative_exp)
                            - ((ht.next_values.mu_snp / ht.total_mu) * ht.total_exp),
                        ),
                        obs_expr=(
                            ht.next_values.cum_obs - ht.cumulative_obs[ht.transcript]
                        )
                        - ht.next_values.obs,
                        exp_expr=(ht.next_values.exp - ht.cumulative_exp)
                        - ((ht.next_values.mu_snp / ht.total_mu) * ht.total_exp),
                    ),
                    # Get alt expression for section of transcript post-window
                    get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=get_obs_exp_expr(
                            cond_expr=True,
                            obs_expr=ht.next_values.reverse_obs + ht.next_values.obs,
                            exp_expr=ht.next_values.reverse_exp
                            + ((ht.next_values.mu_snp / ht.total_mu) * ht.total_exp),
                        ),
                        obs_expr=ht.next_values.reverse_obs + ht.next_values.obs,
                        exp_expr=ht.next_values.reverse_exp
                        + ((ht.next_values.mu_snp / ht.total_mu) * ht.total_exp),
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
    ht = ht.annotate(chisq=(2 * (hl.log(ht.total_alt) - hl.log(ht.total_null))))

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
    Annotates each position in Table with whether that position is a significant breakpoint.

    Also annotates input Table with cumulative observed, expected, and observed/expected values
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
    Annotates total observed, expected, and observed/expected (OE) counts for each section of a transcript.

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
    Returns average number of bases between observed missense variation.

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
    logger.info("Getting total number of bases in the exome (based on GENCODE)...")
    total_bases = get_exome_bases(build=get_reference_genome(ht.locus).name)
    total_variants = ht.count()
    logger.info(f"Total number of bases in the exome: {total_bases}")
    logger.info(f"Total number of missense variants in gnomAD exomes: {total_variants}")

    logger.info("Getting average bases between missense variants and returning...")
    return round(total_bases / total_variants)


# NOTE: renaming this function because don't want to delete yet
def search_for_two_breaks_prev(
    ht: hl.Table,
    exome_ht: hl.Table,
    chisq_threshold: float = 13.8,
    num_obs_var: int = 10,
) -> hl.Table:
    """
    Searches for evidence of constraint within a set window size/number of base pairs.

    Function is designed to search in transcripts that didn't have one single significant break.

    Assumes that:
        - Input Table has a field named 'transcript'.

    :param hl.Table ht: Input Table.
    :param hl.Table exome_ht: Table containing variants from gnomAD exomes.
    :param float chisq_threshold: Chi-square significance threshold. 
        Value should be 10.8 (single break) and 13.8 (two breaks) (values from ExAC RMC code).
    :param int num_obs_var: Number of observed variants. Used when determining the window size for simultaneous breaks. 
        Default is 10, meaning that the window size for simultaneous breaks is the average number of base pairs required to see 10 observed variants.
    :return: Table annotated with is_break at the *end* position of a simultaneous break window.
    :rtype: hl.Table
    """
    break_size = get_avg_bases_between_mis(exome_ht) * num_obs_var
    logger.info(
        f"Number of bases to search for constraint (size for simultaneous breaks): {break_size}"
    )

    logger.info(
        "Annotating each row with previous row's forward observed, expected, and OE values..."
    )
    # Previous row's forward obs, exp, OE values reflect the previous transcript section's total obs, exp, OE values
    ht = ht.annotate(position=hl.format("%s_%s", ht.transcript, ht.locus.position))
    ht = ht.annotate(
        prev_values=hl.scan.group_by(
            ht.position,
            hl.struct(
                prev_obs=hl.scan._prev_nonnull(ht.cumulative_obs[ht.transcript]),
                prev_exp=hl.scan._prev_nonnull(ht.cumulative_exp),
                prev_oe=hl.scan._prev_nonnull(ht.forward_oe),
                prev_pos=hl.scan._prev_nonnull(ht.locus.position),
            ),
        )
    )

    logger.info(
        f"Annotating each position with start position for window \
        if position - {break_size - 1} bases is greater than or equal to the transcript start pos..."
    )
    ht = ht.annotate(
        window_start=hl.or_missing(
            # Window start is defined only if position - break_size + 1 is greater than or
            # equal to the transcript start position
            # e.g., if position is 5, break size is 3, and start position is 1, then window_start is 5 - (3-1) = 3
            (ht.locus.position - (break_size - 1) >= ht.start_pos),
            ht.locus.position - (break_size - 1),
        )
    )

    # Update window start, as start position isn't necessarily present in context HT
    # (context HT filtered to loci with possible missense variants only)
    # Get list of all positions
    all_pos = ht.aggregate(hl.agg.collect(ht.locus.position), _localize=False)

    def _get_closest_pos(
        pos: hl.expr.Int64Expression, all_pos: hl.expr.ArrayExpression
    ) -> int:
        """
        Checks list of all window start positions and finds closest position in HT for each window start.

        This function assumes all_pos is sorted 
        (and all_pos should be sorted, if it is created from a table keyed by locus and alleles).
        
        In the event of a tie, this function will return the *smallest* position that is 
        closest to the position of interest.
        For example, if the position of interest is 4, and the list of positions is [1, 7, 10, 14],
        the closest positions are either 1 or 7. This function will return 1.

        :param hl.expr.Int64Expression pos: Position of interest.
        :param hl.expr.ArrayExpression all_pos: ArrayExpression containing list of all positions.      
        """

        def _return_closest(
            l: hl.expr.Int64Expression, r: hl.expr.Int64Expression
        ) -> hl.expr.Int64Expression:
            """
            Function to return whether left number (l) or right number (r) is closer to input position.

            `l` and `r` are both numbers from an ArrayExpression. 
            `l` is the number with the smaller index (on the left).
            `r` is the number with the larger index (on the right).

            :param hl.expr.Int64Expression l: Left number.
            :param hl.expr.Int64Expression r: Right number.
            :return: `l` or `r`, whichever is closest to the input position `pos`.
            :rtype: hl.expr.Int64Expression
            """
            # NOTE: this has to be less than or equal to to return the *smaller* position in the event of a tie
            # Just less than will return the *larger* number
            # Solution adapted from https://discuss.hail.is/t/find-value-closest-to-input-in-list/1903
            return hl.if_else(hl.abs(l - pos) <= hl.abs(r - pos), l, r)

        return hl.fold(_return_closest, all_pos[0], all_pos[1:])

    logger.info("Annotating HT with position closest to each window start...")
    ht = ht.annotate(
        pre_window_pos=hl.or_missing(
            hl.is_defined(ht.window_start), _get_closest_pos(ht.locus.position, all_pos)
        )
    )
    # Check if pre window pos is ever larger than window_start
    check_start = ht.aggregate(hl.agg.count_where(ht.pre_window_pos > ht.window_start))
    if check_start > 0:
        raise DataException(
            f"Position closest to window start for HT is larger than window start position in {check_start} cases!"
        )
    ht = ht.checkpoint(f"{temp_path}/simul_break_prep.ht", overwrite=True)

    return search_for_break(
        ht, "transcript", simul_break=True, chisq_threshold=chisq_threshold
    )


def search_for_two_breaks(
    ht: hl.Table,
    exome_ht: hl.Table,
    chisq_threshold: float = 13.8,
    num_obs_var: int = 10,
) -> hl.Table:
    """
    Searches for evidence of constraint within a set window size/number of base pairs.

    Function is designed to search in transcripts that didn't have one single significant break.

    Assumes that:
        - Input Table has a field named 'transcript'.

    :param hl.Table ht: Input Table.
    :param hl.Table exome_ht: Table containing variants from gnomAD exomes.
    :param float chisq_threshold: Chi-square significance threshold. 
        Value should be 10.8 (single break) and 13.8 (two breaks) (values from ExAC RMC code).
    :param int num_obs_var: Number of observed variants. Used when determining the window size for simultaneous breaks. 
        Default is 10, meaning that the window size for simultaneous breaks is the average number of base pairs required to see 10 observed variants.
    :return: Table annotated with is_break at the *end* position of a simultaneous break window.
    :rtype: hl.Table
    """
    break_size = get_avg_bases_between_mis(exome_ht) * num_obs_var
    logger.info(
        f"Number of bases to search for constraint (size for simultaneous breaks): {break_size}"
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

    logger.info(
        "Checking how many window end positions are actually defined in the input HT..."
    )
    # The input HT contains every possible missense variant
    # and will be missing positions if no missense variant was possible at that position
    # Duplicate HT to check window ends
    window_ht = ht.select()
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
    ht = ht.select("window_end", "end_pos")

    logger.info(
        "Checkpointing HT with sites that have their window ends defined in the HT..."
    )
    end_ht = ht.filter(
        hl.is_defined(window_ht[hl.locus(ht.locus.contig, ht.window_end)])
    )
    end_ht = end_ht.annotate(post_window_pos=end_ht.window_end)
    end_ht = end_ht.checkpoint(
        f"{temp_path}/simul_break_prep_end_def.ht", overwrite=True
    )
    logger.info(f"Sites with defined window ends: {end_ht.count()}")

    logger.info("Working on sites that don't have their window ends defined in the HT")
    no_end_ht = ht.filter(
        hl.is_missing(window_ht[hl.locus(ht.locus.contig, ht.window_end)])
    )
    no_end_ht = no_end_ht.checkpoint(
        f"{temp_path}/simul_break_prep_no_end.ht", overwrite=True
    )
    no_end_transcripts = no_end_ht.aggregate(
        hl.agg.collect_as_set(no_end_ht.transcript), _localize=False
    )
    logger.info(f"Sites without defined window ends: {no_end_ht.count()}")

    logger.info("Gathering all positions in each transcript...")
    pos_ht = window_ht.filter(no_end_transcripts.contains(window_ht.transcript))
    pos_ht = pos_ht.group_by(transcript=pos_ht.transcript).aggregate(
        positions=hl.agg.collect(pos_ht.locus.position)
    )
    pos_ht = pos_ht.checkpoint(f"{temp_path}/pos_per_transcript.ht", overwrite=True)

    # Annotate positions per transcript back onto no end HT
    no_end_ht = no_end_ht.annotate(
        pos_per_transcript=pos_ht[no_end_ht.transcript].positions
    )
    no_end_ht = no_end_ht.annotate(
        post_window_index=hl.binary_search(
            no_end_ht.pos_per_transcript, no_end_ht.window_end
        ),
        n_pos_per_transcript=hl.len(no_end_ht.pos_per_transcript),
    )
    no_end_ht = no_end_ht.annotate(
        post_window_pos=hl.case()
        .when(
            hl.is_defined(no_end_ht.post_window_index),
            hl.if_else(
                no_end_ht.post_window_index == no_end_ht.n_pos_per_transcript,
                no_end_ht.end_pos,
                no_end_ht.pos_per_transcript[no_end_ht.post_window_index],
            ),
        )
        .or_missing()
    )
    no_end_ht = no_end_ht.checkpoint(
        f"{temp_path}/simul_break_prep_no_end_ready.ht", overwrite=True
    )

    logger.info("Joining no end HT with end HT...")
    ht = end_ht.join(no_end_ht, how="outer")
    ht = ht.annotate(
        post_window_pos=hl.if_else(
            hl.is_defined(ht.post_window_pos), ht.post_window_pos, ht.post_window_pos_1,
        ),
        **annotation_ht[ht.key],
    )
    ht = ht.checkpoint(f"{temp_path}/simul_break_prep.ht", overwrite=True)
    logger.info(f"HT count: {ht.count()}")

    # Check if post window pos is ever smaller than window_end
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

    logger.info("Searching for two breaks...")
    return search_for_break(
        ht, "transcript", simul_break=True, chisq_threshold=chisq_threshold
    )


def calculate_section_chisq(
    obs_expr: hl.expr.Int64Expression, exp_expr: hl.expr.Float64Expression,
) -> hl.expr.Float64Expression:
    """
    Creates expression checking if transcript section is significantly different than the null model (no evidence of regional missense constraint).

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
    Returns struct with constraint flags.

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
    Fixes observed and expected counts for XG (gene that spans PAR and non-PAR regions on chrX).

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
        Fixes total expected and total mu counts for XG.

        :param hl.Table xg: Context Table filtered to XG.
        :param str xg_transcript: XG transcript string. Default is 'ENST00000419513'.
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

    def _fix_xg_obs(xg: hl.Table, exome_ht: hl.Table):
        """
        Fixes total observed counts for XG.

        :param hl.Table xg: Context Table filtered to XG.
        :param hl.Table exome_ht: Table containing variants from gnomAD exomes.
        :return: Context Table filtered to XG, annotated with total observed and observed per position values.
        :rtype: hl.Table
        """
        exome_ht = exome_ht.filter(exome_ht.transcript == xg_transcript)
        obs_ht = calculate_observed(exome_ht)
        xg = xg.annotate(_obs=exome_ht.index(xg.key))
        xg = xg.transmute(observed=hl.int(hl.is_defined(xg._obs)))
        return xg.annotate(total_obs=obs_ht[xg.transcript].observed)

    logger.info(f"Filtering context HT to XG (transcript: {xg_transcript})...")
    xg = context_ht.filter(context_ht.transcript == xg_transcript)

    logger.info("Fixing expected counts for XG...")
    exp_struct = _fix_xg_exp(xg, groupings)

    logger.info("Fixing observed counts for XG...")
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
    xg = xg.annotate(total_exp=exp_struct.total_exp, total_mu=exp_struct.total_mu)
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
