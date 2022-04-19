import logging
from typing import Dict, List, Set, Tuple, Union

import hail as hl

from gnomad.utils.file_utils import file_exists
from gnomad.resources.resource_utils import DataException

from gnomad_lof.constraint_utils.generic import annotate_variant_types

from rmc.resources.basics import (
    constraint_prep,
    multiple_breaks,
    oe_bin_counts_tsv,
    rmc_browser,
    rmc_results,
    simul_break,
    temp_path,
    TOTAL_EXOME_BASES,
)
from rmc.resources.grch37.reference_data import (
    clinvar_path_mis,
    de_novo,
)
from rmc.utils.generic import (
    get_coverage_correction_expr,
    get_exome_bases,
    import_clinvar_hi_variants,
    import_de_novo_variants,
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

MULTI_BREAK_ANNOTATIONS = [
    "mu_snp",
    "observed",
    "total_exp",
    "total_mu",
    "total_obs",
    "cumulative_obs",
    "cumulative_exp",
    "break_list",
    "start_pos",
    "end_pos",
]
"""
List of annotations to keep when finalizing first/additional (multiple) breaks HT.
"""

SIMUL_BREAK_ANNOTATIONS = [
    "max_chisq",
    "break_pos",
    "start_pos",
    "end_pos",
]
"""
List of annotations to keep when finalizing two simultaneous breaks HT.
"""

CONSTRAINT_ANNOTATIONS = [
    "section_start_pos",
    "section_end_pos",
    "section",
    "section_exp",
    "section_obs",
    "section_chisq",
    "section_oe",
]
"""
List of annotations to keep on finalized regional missense constraint results HT.
"""


def calculate_observed(ht: hl.Table) -> hl.Table:
    """
    Group input Table by transcript, filter based on `keep_criteria`, and aggregate observed variants count per transcript.

    .. note::
        Assumes input HT has been filtered using `keep_criteria`.

    :param hl.Table ht: Input Table.
    :return: Table annotated with observed variant counts.
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
    """
    return (cumulative_mu_expr[transcript_expr] / total_mu_expr) * total_exp_expr


def calculate_exp_per_transcript(
    context_ht: hl.Table,
    locus_type: str,
    groupings: List[str] = GROUPINGS,
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
    """
    logger.info("Grouping by %s...", groupings)
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
    """
    return hl.fold(lambda i, j: i * j, 1, dpois_expr)


def search_for_break(
    ht: hl.Table,
    search_field: hl.str,
    chisq_threshold: float = 6.6,
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
        Default is 6.6 (corresponds to p-value of 0.01).
        Default in ExAC was 10.8 (corresponds to p-value of 0.001).
        Reference: https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm
    :return: Table annotated with whether position is a breakpoint.
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
    group_ht = group_ht.checkpoint(f"{temp_path}/max_chisq.ht", overwrite=True)
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
        Value should be 6.6 (corresponds to p-value of 0.01).
    :return: Table with cumulative observed, expected, and observed/expected values annotated for forward and reverse directions.
        Table also annotated with boolean for whether each position is a breakpoint.
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
        section_oe=get_obs_exp_expr(
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
        Value should be 6.6 (corresponds to p-value of 0.01).
    :return: Table annotated with whether position is a breakpoint.
    """
    ht = get_subsection_exprs(ht)

    # Rename section_oe (each section's observed/expected value) to be overall_oe
    # This is because the overall OE ratio used in searching for additional breaks should be
    # transcript section OE and not transcript overall OE
    ht = ht.transmute(overall_oe=ht.section_oe)

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
        chisq_threshold=chisq_threshold,
    )
    post_ht = search_for_break(
        post_ht,
        search_field="transcript",
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
        Value should be 6.6 (corresponds to p-value of 0.01).
    :return: Table annotated with whether position is a breakpoint.
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
    """
    return ((obs_expr - exp_expr) ** 2) / exp_expr


def get_section_info(
    ht: hl.Table,
    section_num: int,
    section_type: str,
    indices: Tuple[int],
    simul_breaks: bool,
) -> hl.Table:
    """
    Get the number of observed variants, number of expected variants, and chi square value for transcript section.

    Section observed and expected missense calculations are slightly different between the first/addtional breaks results and the
    two simultaneous breaks results. This is because the breakpoints are treated differently:

    For first/additional breaks, breakpoints are treated as part of the previous transcript subsection.
    For example, for a transcript with 200 positions and a single breakpoint between 99 and 100,
    the transcript would be divided into the sections [1,99], (99,200].

    For two simultaneous breaks, the two breakpoints are treated as a section.
    For example, for a transcript with 200 positions, and two breaks at 50 and 100,
    the transcript would be divided into the sections [1, 50), [50, 100], (100, 200].

    .. note::
        Assumes that the input Table is annotated with a list of breakpoint positions (`break_pos`) and
        with each transcript's start and end positions (`start_pos`, `end_pos`).

    :param hl.Table ht: Input Table.
    :param int section_num: Transcript section number (e.g., 1 for first section, 2 for second, 3 for third, etc.).
    :param str section_type: Transcript section type. Must be one of 'first', 'middle', or 'last'.
    :param Tuple[int] indices: List of indices pointing to breakpoints.
    :param bool simul_breaks: Whether this function is getting section information for two simultaneous breaks.
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
        if simul_breaks:
            ht = ht.filter(ht.locus.position < ht.break_pos[0])
            ht = ht.annotate(
                section_start_pos=ht.start_pos,
                section_end_pos=ht.break_pos[0] - 1,
            )
        else:
            ht = ht.filter(ht.locus.position <= ht.break_pos[0])
            ht = ht.annotate(
                section_start_pos=ht.start_pos,
                section_end_pos=ht.break_pos[0],
            )
    elif section_type == "middle":
        if simul_breaks:
            ht = ht.filter(
                (ht.locus.position >= ht.break_pos[indices[0]])
                & (ht.locus.position <= ht.break_pos[indices[1]])
            )
            ht = ht.annotate(
                section_start_pos=ht.break_pos[indices[0]],
                section_end_pos=ht.break_pos[indices[1]],
            )
        else:
            ht = ht.filter(
                (ht.locus.position > ht.break_pos[indices[0]])
                & (ht.locus.position <= ht.break_pos[indices[1]])
            )
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
    ht: hl.Table, max_n_breaks: int, simul_breaks: bool
) -> hl.Table:
    """
    Annotate each transcript section with observed, expected, OE, and section chi square values.

    .. note::
        Needs to be run for each break number. For example, this function needs to be run for transcripts with three breaks,
        and it needs to be run again for transcripts with four breaks.

    :param hl.Table ht: Input Table.
    :param int max_n_breaks: Largest number of breaks.
    :param bool simul_breaks: Whether this function is getting section information for two simultaneous breaks.
    :return: Table with section and section values annotated.
    """
    logger.info("Getting section information for first section of each transcript...")
    count = 1
    section_ht = get_section_info(
        ht,
        section_num=count,
        section_type="first",
        indices=None,
        simul_breaks=simul_breaks,
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
                simul_breaks=simul_breaks,
            )
            section_ht = section_ht.union(temp_ht)
            count += 1
    end_ht = get_section_info(
        ht,
        section_num=count,
        section_type="last",
        indices=None,
        simul_breaks=simul_breaks,
    )
    return section_ht.union(end_ht)


def get_unique_transcripts_per_break(
    ht: hl.Table,
) -> Tuple[Dict[int, Union[Set[str], hl.expr.SetExpression]], int]:
    """
    Return a DictExpression with break numbers (key) and set of transcripts unique to each break number (value).

    Key is in format `break_{number}_transcripts`. Value is empty SetExpression if that break number has no unique transcripts.

    Function also checkpoints temporary table containing breakpoint positions and number of breaks for each transcript.

    .. note::
        - This function will only get unique transcripts for transcripts with one or one + additional breaks.
        - Assumes input Table is annotated with list containing booleans for whether that locus is a breakpoint
        (`break_list`).

    :param hl.Table: Input Table.
    :return: Tuple of dictionary with break number (key) and set of transcripts unique to that break number or empty SetExpression (value) and
        maximum number of breaks for any transcript.
    """
    transcripts_per_break = {}
    # Use `break_list` annotation (list of booleans for whether a row is a breakpoint)
    # to filter one/additional breaks ht to rows that are significant breakpoints ONLY
    ht = ht.filter(ht.break_list.any(lambda x: x))

    # Group HT (filtered to breakpoint positions only) by transcript and
    # count the number of breakpoints associated with each transcript
    group_ht = ht.group_by("transcript").aggregate(
        break_pos=hl.sorted(hl.agg.collect(ht.locus.position))
    )
    group_ht = group_ht.annotate(n_breaks=hl.len(group_ht.break_pos))

    # Checkpoint to force hail to finish this group_by computation
    group_ht = group_ht.checkpoint(
        f"{temp_path}/breaks_per_transcript.ht", overwrite=True
    )
    max_n_breaks = hl.eval(group_ht.aggregate(hl.agg.max(group_ht.n_breaks)))

    for i in range(1, max_n_breaks + 1):
        temp_ht = group_ht.filter(group_ht.n_breaks == i)
        transcripts = temp_ht.aggregate(hl.agg.collect_as_set(temp_ht.transcript))
        if len(transcripts > 0):
            transcripts_per_break[f"break_{i}_transcripts"] = transcripts
        else:
            transcripts_per_break[f"break_{i}_transcripts"] = hl.missing(
                hl.tset(hl.tstr)
            )

    # Save transcripts per break DictExpression
    hl.experimental.write_expression(
        transcripts_per_break, f"{temp_path}/transcripts_per_break.he"
    )
    return (transcripts_per_break, max_n_breaks)


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
    """
    ht = ht.annotate(
        regions=hl.struct(
            start=ht.section_start_pos,
            stop=ht.section_end_pos,
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
    multi_breaks_ht: hl.Table,
    annotations: List[str] = MULTI_BREAK_ANNOTATIONS,
) -> hl.Table:
    """
    Organize table of transcripts with multiple breaks.

    Get number of transcripts unique to each break number and drop any extra annotations.
    Also calculate each section's observed, expected, OE, and chi-square values.

    Assumes:
        - Table is annotated with set of transcripts per break (e.g., `break_1_transcripts`)

    :param hl.Table multi_breaks_ht: Input Table. Example schema (truncated at two breaks for space reasons):
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
    :param List[str] annotations: List of annotations to keep from input Table.
        Default is MULTI_BREAK_ANNOTATIONS.
    :return: Table annotated with transcript subsection values and breakpoint positions.
    """
    logger.info("Getting transcripts associated with each break number...")
    # Get number of transcripts unique to each break number
    # Transcript sets in globals of multiple breaks HT are not unique
    # i.e., `break_1_transcripts` could contain transcripts that are also present in `break_2_transcripts`
    multi_breaks_ht = multi_breaks_ht.select_globals(
        multiple_breaks_chisq_threshold=multi_breaks_ht.chisq_threshold
    )
    transcripts_per_break, max_n_breaks = get_unique_transcripts_per_break(
        multi_breaks_ht
    )

    for break_num in transcripts_per_break:
        logger.info(
            "Break number %i has %i transcripts",
            break_num.split("_")[1],
            len(transcripts_per_break[break_num]),
        )

    logger.info("Selecting only relevant annotations from HT and checkpointing...")
    multiple_breaks_transcripts = multi_breaks_ht.aggregate(
        hl.agg.collect_as_set(multi_breaks_ht.transcript)
    )

    # Reading in processed context HT (annotated with cumulative obs/exp missense values)
    # because multiple breaks HT only contains one row per transcript
    # Need the context HT to be able to calculate the observed and expected values for each transcript subsection
    ht = constraint_prep.ht()
    ht = ht.filter(hl.literal(multiple_breaks_transcripts).contains(ht.transcript))
    # This table is checkpointed as part of `get_unique_transcripts_per_break` above
    break_ht = hl.read_table(f"{temp_path}/breaks_per_transcript.ht")
    ht = ht.annotate(break_pos=break_ht[ht.transcript].break_pos)
    ht = ht.select(*annotations)
    ht = ht.checkpoint(f"{temp_path}/multiple_breaks.ht", overwrite=True)

    logger.info("Get transcript section annotations (obs, exp, OE, chisq)...")
    hts = []
    for i in range(1, max_n_breaks + 1):
        transcripts = hl.eval(transcripts_per_break[i])
        if not transcripts:
            logger.info("Break number %i has no transcripts. Continuing...")
            continue

        # Filter HT to transcripts associated with this break only
        # and annotate section information
        temp_ht = ht.filter(hl.literal(transcripts).contains(ht.transcript))
        temp_ht = annotate_transcript_sections(temp_ht, i, simul_breaks=False)

        # Recalculate section oe to remove cap at 1
        # Cap is only necessary when calculating constraint results and should not be included in release
        temp_ht = temp_ht.annotate(section_oe=temp_ht.section_obs / temp_ht.section_exp)
        temp_ht = temp_ht.checkpoint(
            f"{temp_path}/break_{i}_sections.ht",
            overwrite=True,
        )
        hts.append(temp_ht)

    # Merge all HTs together
    ht = hts[0].union(*hts[1:])
    ht = ht.annotate_globals(**transcripts_per_break)
    ht = ht.checkpoint(multiple_breaks.path, overwrite=True)
    return ht


def finalize_simul_breaks(
    simul_ht: hl.Table, annotations: List[str] = SIMUL_BREAK_ANNOTATIONS
) -> hl.Table:
    """
    Process Table containing transcripts with two simultaneous breaks.

    Calculate each transcript subsection observed, expected, OE, and chi-square values.

    :param hl.Table simul_ht: Input Table with simultaneous breaks results.
    :return: Table annotated with transcript subsection values.
    """
    # Rename chisq cutoff in simultaneous breaks HT
    simul_ht = simul_ht.transmute_globals(
        simul_breaks_chisq_threshold=simul_ht.chisq_threshold
    )
    simul_ht = simul_ht.select(*annotations)
    simul_break_transcripts = simul_ht.aggregate(
        hl.agg.collect_as_set(simul_ht.transcript)
    )

    logger.info(
        "Reading in processed context Table, adding simul breaks annotations, and checkpointing..."
    )
    # Reading in processed context HT (annotated with cumulative obs/exp missense values)
    # because simultaneous breaks HT only contains one row per transcript
    # Need the context Table to be able to calculate the observed and expected values for each transcript subsection
    ht = constraint_prep.ht()
    ht = ht.filter(hl.literal(simul_break_transcripts).contains(ht.transcript))
    ht = ht.annotate(
        start_pos=simul_ht[ht.transcript].start,
        end_pos=simul_ht[ht.transcript].stop,
        max_chisq=simul_ht[ht.transcript].max_chisq,
        break_pos=simul_ht[ht.transcript].break_pos,
    )
    ht = ht.checkpoint(f"{temp_path}/simul_breaks.ht", overwrite=True)

    logger.info("Get transcript section annotations (obs, exp, OE, chisq)...")
    ht = annotate_transcript_sections(ht, max_n_breaks=2, simul_breaks=True)
    ht = ht.checkpoint(simul_break.path, overwrite=True)
    return ht


def finalize_all_breaks_results(
    multi_breaks_ht: hl.Table,
    simul_breaks_ht: hl.Table,
    annotations: List[str] = CONSTRAINT_ANNOTATIONS,
) -> None:
    """
    Finalize all breaks results.

    Finalize results for multiple breaks (first/additional breaks) and simultaneous breaks,
    then merge results into single Table.

    Also reformat Table to match desired release HT schema.

    :param hl.Table multi_breaks_ht: Input Table with first/additional (multiple) break results.
    :param hl.Table simul_breaks_ht: Input Table with simultaneous breaks results.
    :param int max_n_breaks: Largest number of breakpoints in any transcript. Used only for multiple breaks results.
    :param List[str] annotations: List of annotations to keep from input Table.
        Default is CONSTRAINT_ANNOTATIONS.
    :return: None; writes Tables to resource path.
    """
    logger.info("Finalizing multiple breaks results...")
    multi_breaks_ht = finalize_multiple_breaks(multi_breaks_ht).select(*annotations)

    logger.info("Finalizing simultaneous breaks results...")
    simul_breaks_ht = finalize_simul_breaks(simul_breaks_ht).select(*annotations)
    ht = multi_breaks_ht.union(simul_breaks_ht, unify=True)
    ht = ht.transmute_globals(
        chisq_threshold=hl.struct(
            multiple_breaks=ht.multiple_breaks_chisq_threshold,
            simultaneous_breaks=ht.simul_breaks_chisq_threshold,
        )
    )
    ht = ht.checkpoint(rmc_results.path, overwrite=True)

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
    """
    return ht1.annotate(**{f"{annot_str}": hl.int(hl.is_defined(ht2[ht1.locus]))})


def get_oe_bins(ht: hl.Table, build: str, get_total_exome_bases: bool = False) -> None:
    """
    Group RMC results HT by obs/exp (OE) bin and annotate.

    Add the following annotations:
        - Proportion coding base pairs
        - Proportion de novo missense from controls/cases
        - Proportion ClinVar pathogenic/likely pathogenic missense variants in haploinsufficient genes.

    Assumes input Table is annotated with the following annotations:
        - `section_start_pos`: Start position for transcript subsection
        - `section_end_pos`: End position for transcript subsection
        - `section_obs`: Number of observed missense variants within transcript subsection
        - `section_exp`: Proportion of expected missense variatns within transcript subsection
        - `section_oe`: Observed/expected missense variation ratio within transcript subsection

    :param hl.Table ht: Input Table containing all breaks results.
    :param str build: Reference genome build.
    :param bool get_total_exome_bases: Boolean for whether to recalculate total number of bases in exome.
        If False, will use value from `TOTAL_EXOME_BASES`. Default is False.
    :return: None; writes TSV with OE bins + annotations to `oe_bin_counts_tsv` resource path.
    """
    if build != "GRCh37":
        raise DataException(
            "ClinVar and de novo files currently only exist for GRCh37!"
        )

    logger.info("Reading in ClinVar, de novo missense, and transcript HTs...")
    if not file_exists(clinvar_path_mis.path):
        import_clinvar_hi_variants(build="GRCh37", overwrite=True)
    if not file_exists(de_novo.path):
        import_de_novo_variants(build="GRCh37", overwrite=True)

    clinvar_ht = clinvar_path_mis.ht()
    dn_ht = de_novo.ht()

    # Split de novo HT into two HTs -- one for controls and one for cases
    dn_controls_ht = dn_ht.filter(dn_ht.case_control == "control")
    dn_case_ht = dn_ht.filter(dn_ht.case_control != "control")

    # Get total number of coding base pairs, also ClinVar and DNM variants
    if get_total_exome_bases:
        total_bp = get_exome_bases(build=build)
    else:
        total_bp = TOTAL_EXOME_BASES
    total_clinvar = clinvar_ht.count()
    total_control = dn_controls_ht.count()
    total_case = dn_case_ht.count()

    # Annotate with control and case DNM, ClinVar P/LP variants,
    ht = check_loci_existence(ht, dn_controls_ht, "dnm_controls")
    ht = check_loci_existence(ht, dn_case_ht, "dnm_cases")
    ht = check_loci_existence(ht, clinvar_ht, "clinvar_path")
    ht = ht.annotate(
        oe_bin=hl.case()
        .when(ht.section_oe <= 0.2, "0-0.2")
        .when((ht.section_oe > 0.2) & (ht.section_oe <= 0.4), "0.2-0.4")
        .when((ht.section_oe > 0.4) & (ht.section_oe <= 0.6), "0.4-0.6")
        .when((ht.section_oe > 0.6) & (ht.section_oe <= 0.8), "0.6-0.8")
        .default("0.8-1.0")
    )
    # Checkpoint HT here because it becomes a couple tables below
    ht = ht.checkpoint(f"{temp_path}/breaks_oe_bin.ht", overwrite=True)

    # Group HT by section to get number of base pairs per section
    # Need to group by to avoid overcounting
    group_ht = ht.group_by("section").aggregate(
        start=hl.agg.take(ht.section_start_pos, 1)[0],
        end=hl.agg.take(ht.section_end_pos, 1)[0],
        obs=hl.agg.take(ht.section_obs, 1)[0],
        exp=hl.agg.take(ht.section_exp, 1)[0],
        chisq=hl.agg.take(ht.section_chisq, 1)[0],
        oe=hl.agg.take(ht.section_oe, 1)[0],
        oe_bin=hl.agg.take(ht.oe_bin, 1)[0],
    )
    group_ht = group_ht.checkpoint(f"{temp_path}/sections.ht", overwrite=True)
    group_ht = group_ht.annotate(bp=group_ht.end - group_ht.start)
    group_ht = group_ht.group_by("oe_bin").aggregate(bp_sum=hl.agg.sum(group_ht.bp))

    # Group HT
    assess_ht = ht.group_by("oe_bin").aggregate(
        dnm_controls=hl.agg.sum(ht.dnm_controls),
        dnm_case=hl.agg.sum(ht.dnm_cases),
        clinvar=hl.agg.sum(ht.clinvar_path),
    )
    assess_ht_count = assess_ht.count()
    if assess_ht_count != 5:
        raise DataException(
            f"Expected 5 OE bins but found {assess_ht_count}. Please double check and rerun!"
        )

    logger.info("Reformatting annotations on assessment HT to be proportions...")
    assess_ht = assess_ht.annotate(bp=group_ht[assess_ht.key])
    assess_ht = assess_ht.transmute(
        bp=assess_ht.bp / total_bp,
        controls=assess_ht.dnm_controls / total_control,
        case=assess_ht.dnm_case / total_case,
        clinvar=assess_ht.clinvar / total_clinvar,
    )
    # Add a show to double check HT
    assess_ht.show()

    logger.info("Exporting assessment HT (in preparation for plotting)...")
    assess_ht.export(oe_bin_counts_tsv)


# TODO: Remove this (expect RMC pipeline to be run after LoF constraint pipeline)
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
