import logging
from typing import Dict, List, Set, Tuple, Union

import hail as hl

from gnomad.resources.grch37.gnomad import coverage, public_release
from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists

from gnomad_lof.constraint_utils.generic import annotate_variant_types

from rmc.resources.basics import (
    multiple_breaks,
    oe_bin_counts_tsv,
    rmc_browser,
    rmc_results,
    temp_path,
)
from rmc.resources.grch37.gnomad import filtered_exomes
from rmc.resources.grch37.reference_data import clinvar_path_mis, de_novo, gene_model
from rmc.utils.generic import (
    get_constraint_transcripts,
    get_coverage_correction_expr,
    keep_criteria,
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
    "mu_scan",
    "total_mu",
    "cumulative_obs",
    "observed",
    "cumulative_exp",
    "total_obs",
    "reverse",
    "forward_oe",
    "section_oe",
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
    "mu_scan",
    "cumulative_exp",
    "break_list",
    "start_pos",
    "end_pos",
]
"""
List of annotations to keep when finalizing release HT.
"""


def add_obs_annotation(
    ht: hl.Table,
    gnomad_data_type: str = "exomes",
    filter_csq: bool = False,
) -> hl.Table:
    """
    Add observed variant count for each variant in input Table.

    Check if locus/allele are present in gnomAD and add as annotation.

    :param ht: Input Table.
    :param gnomad_data_type: gnomAD data type. Used to retrieve public release and coverage resources.
        Must be one of "exomes" or "genomes" (check is done within `public_release`).
        Default is "exomes".
    :param filter_csq: Whether to filter gnomAD data to a specific consequence. Default is False.
        If True, then only variants that match this consequence in the input Table will be annotated as observed.
    :return: Table with observed variant annotation.
    """
    if filter_csq:
        gnomad_ht = filtered_exomes.ht()

    else:
        logger.info("Adding observed annotation...")
        gnomad_ht = public_release(gnomad_data_type).ht()
        gnomad_cov = coverage(gnomad_data_type).ht()
        gnomad_ht = gnomad_ht.select(
            "filters",
            ac=gnomad_ht.freq[0].AC,
            af=gnomad_ht.freq[0].AF,
            gnomad_coverage=gnomad_cov[gnomad_ht.locus].median,
        )
        gnomad_ht = gnomad_ht.filter(
            keep_criteria(
                gnomad_ht.ac, gnomad_ht.af, gnomad_ht.filters, gnomad_ht.gnomad_coverage
            )
        )

    ht = ht.annotate(_obs=gnomad_ht.index(ht.key))
    return ht.transmute(observed=hl.int(hl.is_defined(ht._obs)))


def calculate_observed(ht: hl.Table) -> hl.Table:
    """
    Group input Table by transcript, filter based on `keep_criteria`, and aggregate observed variants count per transcript.

    .. note::
        Assumes input HT has been filtered using `keep_criteria`.

    :param ht: Input Table.
    :return: Table annotated with observed variant counts.
    """
    return ht.group_by(ht.transcript_consequences.transcript_id).aggregate(
        observed=hl.agg.count()
    )


def get_cumulative_obs_expr(
    group_str: hl.expr.StringExpression,
    observed_expr: hl.expr.Int64Expression,
) -> hl.expr.DictExpression:
    """
    Return annotation with the cumulative number of observed variants (non-inclusive).

    Value is non-inclusive (does not include value of row)
    due to the nature of `hl.scan` and needs to be corrected later.

    This function can produce the scan when searching for the first break or when searching for additional break(s).

    :param group_str: StringExpression containing transcript or transcript subsection information.
        Used to group observed and expected values.
    :param observed_expr: Observed variants expression.
    :return: DictExpression containing scan expressions for cumulative observed variant counts for `search_expr`.
    """
    return hl.scan.group_by(group_str, hl.scan.sum(observed_expr))


def adjust_obs_expr(
    cumulative_obs_expr: hl.expr.DictExpression,
    obs_expr: hl.expr.Int64Expression,
    group_str: str = "section",
) -> hl.expr.Int64Expression:
    """
    Adjust the scan with the cumulative number of observed variants.

    This adjustment is necessary because scans are always one line behind, and we want the values to match per line.

    This function can correct the scan created when searching for the first break or when searching for additional break(s).

    .. note::
        This function expects that `cumulative_obs_expr` is a DictExpression keyed by `group_str`.

    :param cumulative_obs_expr: DictExpression containing scan expression with cumulative observed counts per base.
    :param obs_expr: IntExpression with value of either 0 (no observed variant at site) or 1 (variant found in gnomAD).
    :param group_str: Field used to group Table observed and expected values. Default is "section".
    :return: Adjusted cumulative observed counts expression.
    """
    return hl.if_else(
        # Check if the current transcript/section exists in the obs_scan dictionary
        # If it doesn't exist, that means this is the first line in the HT for that particular transcript
        # The first line of a scan is always missing, but we want it to exist
        # Thus, set the cumulative_obs equal to the current observed value
        hl.is_missing(cumulative_obs_expr.get(group_str)),
        obs_expr,
        # Otherwise, add the current obs to the scan to make sure the cumulative value isn't one line behind
        cumulative_obs_expr[group_str] + obs_expr,
    )


def get_cumulative_mu_expr(
    group_str: hl.expr.StringExpression,
    mu_expr: hl.expr.Float64Expression,
) -> hl.expr.DictExpression:
    """
    Return annotation with the cumulative mutation rate probability (non-inclusive).

    Value is non-inclusive due to the nature of `hl.scan` and needs to be corrected later.

    This function can produce the scan when searching for the first break or when searching for additional break(s).

    :param group_str: StringExpression containing transcript or transcript subsection information.
        Used to group observed and expected values.
    :param mu_expr: FloatExpression containing mutation rate probability per variant.
    :return: DictExpression containing scan expressions for cumulative mutation rate probability.
    """
    return hl.scan.group_by(group_str, hl.scan.sum(mu_expr))


def adjust_mu_expr(
    cumulative_mu_expr: hl.expr.DictExpression,
    mu_expr: hl.expr.Float64Expression,
    group_str: hl.expr.StringExpression,
) -> hl.expr.Float64Expression:
    """
    Adjust the scan with the cumulative number mutation rate probability.

    This adjustment is necessary because scans are always one line behind, and we want the values to match per line.

    This function can correct the scan created when searching for the first break or when searching for additional break(s).

    :param cumulative_mu_expr: DictExpression containing scan expressions for cumulative mutation rate probability.
    :param mu_expr: FloatExpression containing mutation rate probability per variant.
    :param group_str: StringExpression containing transcript or transcript subsection information.
        Used to group observed and expected values.
    :return: Adjusted cumulative mutation rate expression.
    """
    return hl.if_else(
        # Check if the current transcript/section exists in the mu_scan dictionary
        # If it doesn't exist, that means this is the first line in the HT for that particular transcript/section
        # The first line of a scan is always missing, but we want it to exist
        # Thus, set the cumulative_mu equal to the current mu_snp value
        hl.is_missing(cumulative_mu_expr.get(group_str)),
        mu_expr,
        # Otherwise, add the current mu_snp to the scan to make sure the cumulative value isn't one line behind
        cumulative_mu_expr[group_str] + mu_expr,
    )


def translate_mu_to_exp_expr(
    cumulative_mu_expr: hl.expr.DictExpression,
    total_mu_expr: hl.expr.Float64Expression,
    total_exp_expr: hl.expr.Float64Expression,
) -> hl.expr.DictExpression:
    """
    Translate cumulative mutation rate probability into cumulative expected count per base.

    Expected variants counts are produced per base by first calculating the fraction of probability of mutation per base,
    then multiplying that fraction by the total expected variants count for a transcript or transcript sub-section.

    The expected variants count for section of interest is mutation rate per SNP adjusted by location in the genome/CpG status
    (plateau model) and coverage (coverage model).

    :param cumulative_mu_expr: DictExpression containing scan expressions for cumulative mutation rate probability.
    :param total_mu_expr: FloatExpression describing total sum of mutation rate probabilities per transcript.
    :param total_exp_expr: FloatExpression describing total expected variant counts per transcript.
    :return: Cumulative expected variants count expression.
    """
    return (cumulative_mu_expr / total_mu_expr) * total_exp_expr


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

    :param context_ht: Context Table.
    :param locus_type: Locus type of input table. One of "X", "Y", or "autosomes".
        NOTE: will treat any input other than "X" or "Y" as autosomes.
    :param groupings: List of Table fields used to group Table to adjust mutation rate.
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

    :param cond_expr: Condition to check prior to adding obs/exp expression.
    :param obs_expr: Expression containing number of observed variants.
    :param exp_expr: Expression containing number of expected variants.
    :return: Observed/expected expression.
    """
    return hl.or_missing(cond_expr, hl.min(obs_expr / exp_expr, 1))


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

    :param total_obs_expr: Expression containing total number of observed variants for transcript.
    :param total_exp_expr: Expression containing total number of expected variants for transcript.
    :param scan_obs_expr: Expression containing cumulative number of observed variants for transcript.
    :param cumulative_exp_expr: Expression containing cumulative number of expected variants for transcript.
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
    obs_str: str,
    mu_str: str,
    total_mu_str: str = "section_mu",
    total_exp_str: str = "section_exp",
    group_str: str = "section",
) -> hl.Table:
    """
    Annotate input Table with the forward section cumulative observed, expected, and observed/expected values.

    .. note::
        'Forward' refers to moving through the transcript from smaller to larger chromosomal positions.

    Expects:
        - Input HT is annotated with transcript, observed, mutation rate, total mutation rate (per section),
        and total expected counts (per section).

    :param hl.Table ht: Input Table.
    :param obs_str: Name of field containing observed variants counts.
    :param mu_str: Name of field containing mutation rate probability per variant.
    :param total_mu_str: Name of field containing total mutation rate per section of interest (transcript or sub-section of transcript).
        Default is 'section_mu'.
    :param total_exp_str: Name of field containing total expected variants count per section of interest (transcript or sub-section of transcript).
        Default is 'section_exp'.
    :param group_str: Field used to group Table observed and expected values. Default is "section".
    :return: Table with forward values (cumulative obs, exp, and forward o/e) annotated.
    """
    logger.info("Getting cumulative observed variant counts...")
    ht = ht.annotate(
        obs_scan=get_cumulative_obs_expr(
            group_str=group_str,
            observed_expr=ht[obs_str],
        )
    )
    ht = ht.annotate(
        cumulative_obs=adjust_obs_expr(
            cumulative_obs_expr=ht.obs_scan,
            obs_expr=ht[obs_str],
            group_str=ht[group_str],
        )
    ).drop("obs_scan")

    logger.info("Getting cumulative expected variant counts...")
    # Get scan of mu_snp
    ht = ht.annotate(mu_scan=get_cumulative_mu_expr(ht[group_str], ht[mu_str]))
    # Adjust scan of mu_snp
    ht = ht.annotate(mu_scan=adjust_mu_expr(ht.mu_scan, ht[mu_str], ht[group_str]))
    ht = ht.annotate(
        cumulative_exp=translate_mu_to_exp_expr(
            ht.mu_scan, ht[total_mu_str], ht[total_exp_str]
        )
    )

    logger.info("Getting forward observed/expected count and returning...")
    # NOTE: adding cond_expr here because get_obs_exp_expr expects it
    # cond_expr is necessary for reverse obs/exp, which is why the function has it
    ht = ht.annotate(cond_expr=True)
    ht = ht.annotate(
        forward_oe=get_obs_exp_expr(
            ht.cond_expr,
            ht.cumulative_obs,
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

    :param ht: Input Table.
    :param total_obs_expr: Expression containing total number of observed variants per transcript (if searching for first break)
        or per section (if searching for additional breaks).
    :param total_exp_expr: Expression containing total number of expected variants per transcript (if searching for first break)
        or per section (if searching for additional breaks).
    :param scan_obs_expr: Expression containing cumulative number of observed variants per transcript
        (if searching for first break) or per section (if searching for additional breaks).
    :param scan_exp_expr: Expression containing cumulative number of expected variants per transcript
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
    Calculate probabilities (in log10 space) of the observed values under a Poisson model with rate given by
    expected * section observed/expected values.

    :param cond_expr: Conditional expression to check before calculating probability.
    :param section_oe_expr: Expression of section observed/expected value.
    :param obs_expr: Expression containing observed variants count.
    :param exp_expr: Expression containing expected variants count.
    :return: log10 of the probability under Poisson model.
    """
    # log_p = True returns the natural logarithm of the probability density
    # Divide this value by hl.log(10) to convert to log base 10
    return hl.or_missing(
        cond_expr,
        hl.dpois(obs_expr, exp_expr * section_oe_expr, log_p=True) / hl.log(10),
    )


def search_for_break(
    ht: hl.Table,
    chisq_threshold: float = 6.6,
    group_str: str = "section",
    min_num_exp_mis: int = 10,
) -> hl.Table:
    """
    Search for breakpoints in a transcript or within a transcript subsection.

    Also checkpoints intermediate table with max chi square values per transcript
    or transcript subsection.

    Expects input HT to contain the following fields:
        - locus
        - section
        - mu_snp
        - cumulative_exp
        - cumulative_obs
        - section_oe
        - forward_oe
        - reverse struct
        - reverse_obs_exp
        - section_obs
        - section_exp

    Also expects:
        - Input HT was created using a VEP context HT.
        - Multiallelic variants in input HT have been split.
        - Input HT is autosomes/PAR only, X non-PAR only, or Y non-PAR only.

    :param ht: Input Table.
    :param chisq_threshold: Chi-square significance threshold.
        Default is 6.6 (single break; p = 0.01).
    :param group_str: Field used to group Table observed and expected values. Default is 'section'.
    :param min_num_exp_mis: Minimum number of expected missense per transcript/transcript section.
        Sections that have fewer than this number of expected missense variants will not
        be computed (chi square will be annotated as -1).
        Default is 10.
    :return: Table annotated with whether position is a breakpoint (`is_break`).
    """
    logger.info(
        "Creating section null (no regional variability in missense depletion)\
        and alt (evidence of domains of missense constraint) expressions..."
    )
    logger.info(
        "Skipping breakpoints that create at least one section\
        that has < %i expected missense variants...",
        min_num_exp_mis,
    )
    # Split transcript or transcript subsection into two sections
    # Split transcript when searching for first break
    # Split transcript subsection when searching for additional breaks
    ht = ht.annotate(
        total_null=hl.or_missing(
            (ht.cumulative_exp >= min_num_exp_mis) & (ht.reverse.exp >= min_num_exp_mis),
            # Add forwards section null (going through positions from smaller to larger)
            # section_null = stats.dpois(section_obs, section_exp*overall_obs_exp)[0]
            get_dpois_expr(
                cond_expr=hl.is_defined(ht.cumulative_obs),
                section_oe_expr=ht.section_oe,
                obs_expr=ht.cumulative_obs,
                exp_expr=ht.cumulative_exp,
            )
            # Add reverse section null (going through positions larger to smaller)
            + get_dpois_expr(
                cond_expr=hl.is_defined(ht.reverse.obs),
                section_oe_expr=ht.section_oe,
                obs_expr=ht.reverse.obs,
                exp_expr=ht.reverse.exp,
            ),
        )
    )
    ht = ht.annotate(
        total_alt=hl.or_missing(
            hl.is_defined(ht.total_null),
            # Add forward section alt
            # section_alt = stats.dpois(section_obs, section_exp*section_obs_exp)[0]
            get_dpois_expr(
                cond_expr=hl.is_defined(ht.cumulative_obs),
                section_oe_expr=ht.forward_oe,
                obs_expr=ht.cumulative_obs,
                exp_expr=ht.cumulative_exp,
            )
            # Add reverse section alt
            + get_dpois_expr(
                cond_expr=hl.is_defined(ht.reverse.obs),
                section_oe_expr=ht.reverse_obs_exp,
                obs_expr=ht.reverse.obs,
                exp_expr=ht.reverse.exp,
            ),
        )
    )

    logger.info("Adding chisq annotation and getting max chisq per section...")
    ht = ht.annotate(
        chisq=hl.or_missing(
            hl.is_defined(ht.total_alt),
            2 * (ht.total_alt - ht.total_null),
        )
    )
    # hl.agg.max ignores NaNs
    group_ht = ht.group_by(group_str).aggregate(max_chisq=hl.agg.max(ht.chisq))
    group_ht = group_ht.checkpoint(f"{temp_path}/max_chisq.ht", overwrite=True)
    ht = ht.annotate(max_chisq=group_ht[ht[group_str]].max_chisq)
    return ht.annotate(
        is_break=((ht.chisq == ht.max_chisq) & (ht.chisq >= chisq_threshold))
    )


def get_subsection_exprs(
    ht: hl.Table,
    section_str: str = "section",
    obs_str: str = "observed",
    mu_str: str = "mu_snp",
    total_mu_str: str = "section_mu",
    total_exp_str: str = "section_exp",
) -> hl.Table:
    """
    Annotate total observed, expected, and observed/expected (OE) counts for each section of a transcript.

    .. note::
        Assumes input Table is annotated with:
            - section
            - observed variants count per site
            - mutation rate probability per site
        Names of annotations must match section_str, obs_str, and mu_str.

    :param ht: Input Table.
    :param section_str: Name of section annotation.
    :param obs_str: Name of observed variant counts annotation.
    :param mu_str: Name of mutation rate probability per site annotation.
    :param total_mu_str: Name of annotation containing sum of mutation rate probabilities per transcript.
    :param total_exp_str: Name of annotation containing total expected variant counts per transcript.
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
        section_mu=section_counts[ht[section_str]].mu,
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


def process_sections(ht: hl.Table, chisq_threshold: float, group_str: str = "section"):
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

    :param ht: Input Table.
    :param chisq_threshold: Chi-square significance threshold.
        Value should be 6.6 (single break) or 9.2 (two breaks) (p = 0.01).
    :param group_str: Field used to group observed and expected values. Default is 'section'.
    :return: Table annotated with whether position is a breakpoint.
    """
    # TODO: When re-running, make sure `get_subsection_exprs`,
    # `get_fwd_exprs` don't run again for first break search only
    # Also rename total to section for this run
    # TODO: update code to stop continually finding first break (we still need a break_list annotation or something similar)
    ht = get_subsection_exprs(ht)

    logger.info(
        "Annotating cumulative observed and expected counts..."
    )
    ht = get_fwd_exprs(
        ht=ht,
        group_str="section",
        obs_str="observed",
        mu_str="mu_snp",
        total_mu_str="section_mu",
        total_exp_str="section_exp",
    )

    logger.info(
        "Annotating reverse observed and expected counts..."
    )
    ht = get_reverse_exprs(
        ht=ht,
        total_obs_expr=ht.section_obs,
        total_exp_expr=ht.section_exp,
        scan_obs_expr=ht.cumulative_obs,
        scan_exp_expr=ht.cumulative_exp,
    )

    logger.info("Searching for a break in each section and returning...")
    ht = search_for_break(
        ht,
        chisq_threshold=chisq_threshold,
    )
    return ht


def calculate_section_chisq(
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
) -> hl.expr.Float64Expression:
    """
    Create expression checking if transcript section is significantly different than the null model (no evidence of regional missense constraint).

    Formula is: (section obs - section exp)^2 / section exp. Taken from ExAC RMC code.

    :param obs_expr: Total number of observed missense variants in section.
    :param exp_expr: Total number of expected missense variants in section.
    :return: Transcript section chi-squared value.
    """
    return ((obs_expr - exp_expr) ** 2) / exp_expr


def get_all_breakpoint_pos(ht: hl.Table) -> hl.GroupedTable:
    """
    Get all breakpoint positions per transcript.

    :param ht: Input Table.
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

    :param ht: Input Table.
    :param section_num: Transcript section number (e.g., 1 for first section, 2 for second, 3 for third, etc.).
    :param section_type: Transcript section type. Must be one of 'first', 'middle', or 'end'.
    :param indices: List of indices pointing to breakpoints.
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

    :param ht: Input Table.
    :param max_n_breaks: Largest number of breaks.
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

    :param ht: Input Table.
    :param max_n_breaks: Largest number of breaks.
    :return: Dictionary with break number (key) and set of transcripts unique to that break number or empty SetExpression (value).
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

    :param ht: Input Table.
    :return: Table with schema described above.
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
    :param max_n_breaks: Largest number of breakpoints in any transcript.
    :param annotations: List of annotations to keep from input Table.
        Default is FINAL_ANNOTATIONS.
    :return: Table annotated with transcript subsection values and breakpoint positions.
    """
    # TODO: Update with new flow
    logger.info("Removing outlier transcripts...")
    outlier_transcripts = get_constraint_transcripts(outlier=True)
    ht = ht.filter(~outlier_transcripts.contains(ht.transcript))

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

    :param ht: Input Table with simultaneous breaks results.
    :return: Table annotated with transcript subsection values.
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

    :param breaks_ht: Input Table with multiple break results.
    :param simul_breaks_ht: Input Table with simultaneous breaks results.
    :param max_n_breaks: Largest number of breakpoints in any transcript. Used only for multiple breaks results.
    :param annotations: List of annotations to keep from input Table.
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

    :param ht1: Table to be annotated.
    :para ht2: Table to check for loci from `ht1`.
    :param str annot_str: Name of annotation to be added to `ht1` designating whether locus is present in `ht2`.
    :return: Annotated version of `ht1`.
    """
    return ht1.annotate(**{f"{annot_str}": hl.int(hl.is_defined(ht2[ht1.locus]))})


def get_oe_bins(ht: hl.Table, build: str) -> None:
    """
    Group RMC results HT by obs/exp (OE) bin and annotate.

    Add the following annotations:
        - Proportion coding base pairs
        - Proportion de novo missense from controls/cases
        - Proportion ClinVar pathogenic/likely pathogenic missense variants in haploinsufficient genes.

    Assumes input Table is annotated with the following annotations:
        - `section_start`: Start position for transcript subsection
        - `section_end`: End position for transcript subsection
        - `section_obs`: Number of observed missense variants within transcript subsection
        - `section_exp`: Proportion of expected missense variatns within transcript subsection
        - `section_oe`: Observed/expected missense variation ratio within transcript subsection

    :param ht: Input Table containing all breaks results.
    :param build: Reference genome build.
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
    transcript_ht = gene_model.ht()

    # Split de novo HT into two HTs -- one for controls and one for cases
    dn_controls_ht = dn_ht.filter(dn_ht.case_control == "control")
    dn_case_ht = dn_ht.filter(dn_ht.case_control != "control")

    # Get total number of coding base pairs, also ClinVar and DNM variants
    # TODO: use exon field in gene model HT to get only coding bases (rather than using transcript end - start)
    transcript_ht = transcript_ht.annotate(
        bp=transcript_ht.end_pos - transcript_ht.start_pos
    )
    total_bp = transcript_ht.aggregate(hl.agg.sum(transcript_ht.bp))
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
        start=hl.agg.take(ht.section_start, 1)[0],
        end=hl.agg.take(ht.section_end, 1)[0],
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

    :param obs_syn_expr: Expression containing number of observed synonymous variants in gnomAD.
    :param obs_mis_expr: Expression containing number of observed missense variants in gnomAD.
    :param obs_lof_expr: Expression containing number of observed loss-of-function (LoF) variants in gnomAD.
    :param exp_syn_expr: Expression containing number of expected synonymous variants.
    :param exp_mis_expr: Expression containing number of expected missense variants.
    :param exp_lof_expr: Expression containing number of expected LoF variants.
    :param raw_syn_z_expr: Expression containing number of Z score for synonymous variants.
    :param raw_mis_z_expr: Expression containing number of Z score for missense variants.
    :param raw_lof_z_expr: Expression containing number of Z score for LoF variants.
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


def group_rmc_ht_by_section(overwrite: bool = False) -> hl.Table:
    """
    Group RMC results Table by transcript subsection and return interval and section missense o/e.

    .. note::
        - Function reads RMC results Table from resource path.
        - Assumes RMC HT is annotated with `locus`, `transcript`, `section`, `section_start_pos`,
        `section_end_pos`, and `section_oe`.
        - Assumes `transcript` is one of RMC HT's key fields.

    :param overwrite: Whether to overwrite temporary checkpointed Table if it exists.
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
