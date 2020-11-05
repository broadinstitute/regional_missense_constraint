import logging
from typing import Dict, List, Union

import hail as hl

from gnomad.utils.reference_genome import get_reference_genome

from gnomad_lof.constraint_utils.generic import annotate_variant_types

from rmc.resources.resource_utils import MISSENSE
from rmc.utils.generic import (
    get_coverage_correction_expr,
    get_exome_bases,
    process_vep,
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
    total_mu_expr: hl.expr.Int32Expression,
    total_exp_expr: hl.expr.Int32Expression,
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
        obs=total_obs_expr - scan_obs_expr, exp=total_exp_expr - cumulative_exp_expr,
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
    ht: hl.Table, search_field: hl.str, chisq_threshold: float = 10.8,
) -> hl.Table:
    """
    Searches for breakpoints in a transcript or within a transcript subsection.

    Expects input context HT to contain the following fields:
        - locus
        - alleles
        - transcript or section
        - coverage (median)
        - mu
        - scan_counts struct
        - overall_obs_exp
        - forward_oe
        - reverse struct
        - reverse_obs_exp
    Also expects:
        - multiallelic variants in input HT have been split.
        - Input HT is autosomes/PAR only, X non-PAR only, or Y non-PAR only.

    Returns HT filtered to lines with maximum chisq if chisq >= max_value, otherwise returns None.

    :param hl.Table ht: Input context Table.
    :param hl.expr.StringExpression search_field: Field of table to search. Value should be either 'transcript' or 'section'. 
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
    ht = ht.annotate(
        section_nulls=[
            # Add forwards section null (going through positions from smaller to larger)
            # section_null = stats.dpois(section_obs, section_exp*overall_obs_exp)[0]
            get_dpois_expr(
                cond_expr=hl.len(ht.cumulative_obs) != 0,
                section_oe_expr=ht.overall_obs_exp,
                obs_expr=ht.cumulative_obs[ht[search_field]],
                exp_expr=ht.cumulative_exp,
            ),
            # Add reverse section null (going through positions larger to smaller)
            get_dpois_expr(
                cond_expr=hl.is_defined(ht.reverse.obs),
                section_oe_expr=ht.overall_obs_exp,
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

    return search_for_break(ht, "transcript", chisq_threshold)


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
    logger.info(
        "Getting total observed and expected counts for each transcript section..."
    )
    # Get total obs per section
    section_obs = ht.group_by(ht.section).aggregate(obs=hl.agg.sum(ht.observed))

    # Get total mu per section -- translate this to total expected per section
    section_mu = ht.group_by(ht.section).aggregate(total=hl.agg.sum(ht.mu_snp))
    ht = ht.annotate(section_total_mu=section_mu[ht.section].total,)
    ht = ht.annotate(
        section_exp=(ht.section_total_mu / ht.total_mu) * ht.total_exp,
        section_obs=section_obs[ht.section].obs,
    )

    logger.info("Getting observed/expected value for each transcript section...")
    ht = ht.annotate(
        break_oe=get_obs_exp_expr(
            cond_expr=hl.is_defined(ht.section),
            obs_expr=ht.section_obs,
            exp_expr=ht.section_exp,
        )
    )

    logger.info("Getting section chi-squared values...")
    ht = ht.annotate(
        section_chisq=calculate_section_chisq(
            obs_expr=ht.break_obs, exp_expr=ht.break_exp,
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
        pre_ht, search_field="section", chisq_threshold=chisq_threshold
    )
    post_ht = search_for_break(
        post_ht, search_field="section", chisq_threshold=chisq_threshold
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
    annot_expr = {f"is_break_{break_num - 1}": ht.is_break}
    ht = ht.annotate(**annot_expr)

    logger.info(
        "Splitting each transcript into two sections: pre-first breakpoint and post..."
    )
    ht = ht.annotate(
        section=hl.if_else(
            # NOTE: Initially added breakpoint position to second section ("post") with logic that
            # first position of scan is always missing
            # Have fixed scans to not be one line behind -- should breakpoint pos still get added to
            # second section?
            ht.locus.position >= break_ht[ht.transcript].locus.position,
            hl.format("%s_%s", ht.transcript, "post"),
            hl.format("%s_%s", ht.transcript, "pre"),
        ),
    )
    return process_sections(ht, chisq_threshold)


def get_avg_bases_between_mis(ht: hl.Table, missense: str = MISSENSE) -> int:
    """
    Returns average number of bases between observed missense variation.

    For example, if the total number of bases is 30, and the total number of missense variants is 10,
    this function will return 3.

    This function is used to determine the minimum size window to check for significant missense depletion
    when searching for two simultaneous breaks.

    :param hl.Table ht: Input gnomAD Table.
    :param str missense: String representing missense variant VEP annotation. Default is MISSENSE.
    :return: Average number of bases between observed missense variants, rounded to the nearest integer,
    :rtype: int
    """
    logger.info("Getting total number of bases in the exome (based on GENCODE)...")
    total_bases = get_exome_bases(build=get_reference_genome(ht.locus).name)

    logger.info(
        "Filtering to missense variants in canonical protein coding transcripts..."
    )
    ht = process_vep(ht, filter_csq=True, csq=missense)
    total_variants = ht.count()
    return round(total_bases / total_variants)


def search_for_two_breaks(
    ht: hl.Table,
    exome_ht: hl.Table,
    transcript: str,
    chisq_threshold: float,
    num_obs_var: int = 10,
) -> Union[hl.Table, None]:
    """
    Searches for evidence of constraint within a set window size/number of base pairs.

    Function is designed to search in transcripts that didn't have one single significant break.
    Currently designed to one run one transcript at a time.

    Assumes that:
        - Input Table has a field named 'transcript'.

    :param hl.Table ht: Input Table.
    :param hl.Table exome_ht: Table containing variants from gnomAD exomes.
    :param str transcript: Transcript of interest.
    :param float chisq_threshold: Chi-square significance threshold. 
        Value should be 10.8 (single break) and 13.8 (two breaks) (values from ExAC RMC code).
    :param int num_obs_var: Number of observed variants. Used when determining the window size for simultaneous breaks. 
        Default is 10, meaning that the window size for simultaneous breaks is the average number of base pairs required to see 10 observed variants.
    :return: Tuple of largest chi-square value and breakpoint positions if significant break was found. Otherwise, None.
    :rtype: Union[hl.Table, None]
    """
    break_size = (
        get_avg_bases_between_mis(
            exome_ht, get_reference_genome(exome_ht.head(1).locus).name
        )
        * num_obs_var
    )
    logger.info(
        f"Number of bases to search for constraint (size for simultaneous breaks): {break_size}"
    )

    # I don't think there's any way to search for simultaneous breaks without a loop?
    ht = ht.filter(ht.transcript == transcript)
    start_pos = ht.head(1).take(1).locus.position
    end_pos = ht.tail(1).take(1).locus.position
    best_chisq = 0
    breakpoints = ()

    breaks = {}
    while start_pos < end_pos:
        ht = ht.annotate(
            section=hl.case()
            # Position is within window if pos is larger than start_pos and
            # less than start_pos + window size, then pos is within window
            .when(
                (ht.locus.position >= start_pos)
                & (ht.locus.position < (start_pos + break_size)),
                hl.format("%s_%s", ht.transcript, "window"),
            )
            # If pos < start_pos, pos is outside of window and pre breaks
            .when(
                ht.locus.position < start_pos, hl.format("%s_%s", ht.transcript, "pre")
            )
            # Otherwise, pos is outside window and post breaks
            .default(hl.format("%s_%s", ht.transcript, "post"))
        )
        new_ht = process_sections(ht, chisq_threshold)

        # Check if found break
        if new_ht.aggregate(hl.agg.counter(new_ht.is_break) > 0):
            max_chisq = new_ht.aggregate(hl.agg.max(new_ht.max_chisq))
            if (max_chisq > best_chisq) and (max_chisq >= chisq_threshold):
                breakpoints = (start_pos, (start_pos + break_size) - 1)
                best_chisq = max_chisq
                annot_expr = {f"{transcript}_breakpoints": breakpoints}
                new_ht = new_ht.annotate_globals(**annot_expr)
                breaks[max_chisq] = [new_ht]

        start_pos += 1

    if best_chisq != 0:
        return breaks[best_chisq]
    return None


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
