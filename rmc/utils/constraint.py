import logging
import re
import subprocess
from typing import Dict, List, Tuple, Union

import hail as hl

from gnomad.resources.grch37.gnomad import coverage, public_release
from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists

from gnomad_lof.constraint_utils.generic import annotate_variant_types

from rmc.resources.basics import (
    SIMUL_BREAK_TEMP_PATH,
    TEMP_PATH,
    TEMP_PATH_WITH_DEL,
)
from rmc.resources.gnomad import filtered_exomes
from rmc.resources.reference_data import clinvar_path_mis, de_novo, gene_model
from rmc.utils.generic import (
    get_coverage_correction_expr,
    keep_criteria,
    import_clinvar_hi_variants,
    import_de_novo_variants,
)
from rmc.resources.rmc import (
    no_breaks,
    oe_bin_counts_tsv,
    simul_search_bucket_path,
    simul_search_round_bucket_path,
    single_search_bucket_path,
    single_search_round_ht_path,
)
from rmc.utils.simultaneous_breaks import merge_simul_break_temp_hts


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
            group_str=ht[group_str],
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
    Calculate probabilities (in log10 space) of the observed values under a Poisson model.

    Rate in model given by expected * section observed/expected values.

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
    search_num: int,
    chisq_threshold: float = 6.6,
    group_str: str = "section",
    min_num_exp_mis: float = 10.0,
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

    :param ht: Input Table.
    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param chisq_threshold: Chi-square significance threshold.
        Default is 6.6 (single break; p = 0.01).
    :param group_str: Field used to group Table observed and expected values. Default is 'section'.
    :param min_num_exp_mis: Minimum number of expected missense per transcript/transcript section.
        Sections that have fewer than this number of expected missense variants will not
        be computed (chi square will be annotated as a missing value).
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
            (ht.cumulative_exp >= min_num_exp_mis)
            & (ht.reverse.exp >= min_num_exp_mis),
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

    ht = ht.checkpoint(
        f"{TEMP_PATH_WITH_DEL}/round{search_num}_all_loci_chisq.ht", overwrite=True
    )

    # hl.agg.max ignores NaNs
    group_ht = ht.group_by(group_str).aggregate(max_chisq=hl.agg.max(ht.chisq))
    group_ht = group_ht.checkpoint(f"{TEMP_PATH_WITH_DEL}/max_chisq.ht", overwrite=True)
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


def process_sections(
    ht: hl.Table,
    search_num: int,
    chisq_threshold: float,
    group_str: str = "section",
):
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
    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
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

    logger.info("Annotating cumulative observed and expected counts...")
    ht = get_fwd_exprs(
        ht=ht,
        group_str=group_str,
        obs_str="observed",
        mu_str="mu_snp",
        total_mu_str="section_mu",
        total_exp_str="section_exp",
    )

    logger.info("Annotating reverse observed and expected counts...")
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
        search_num,
        chisq_threshold=chisq_threshold,
    )
    return ht


def get_rescue_1break_transcripts(
    overwrite: bool, rescue_threshold: float = 5.0
) -> Tuple[hl.expr.SetExpression]:
    """
    Get transcripts that have a single breakpoint in the first round of the 'rescue' search.

    These transcripts did not have a single or simultaneous breakpoint that was over
    the initial search threshold but do have a single breakpoint that is above
    the lower 'rescue' search threshold.

    This function performs the single break search for the first 'rescue' round
    using the statistics already computed during the single break search for the
    first 'initial' round and writes out the resulting tables:
    a section (transcript)-level breakpoint table and a locus-level table
    for transcripts where a 'rescue' breakpoint was found.

    Function uses the no-break results from the first round of the single and
    simultaneous break searches.

    :param overwrite: Whether to overwrite output data if it exists.
    :param rescue_threshold: Lower chi square significance threshold associated with
        the 'rescue' search pathway.
    :return: Tuple set of sections (transcript_start_stop) with two simultaneous breaks in
        the initial search and set of sections with single breakpoint above 'rescue' threshold.
    """
    # Read in the simultaneous breaks results from initial search round 1
    simul_results_path = simul_search_round_bucket_path(
        is_rescue=False,
        search_num=1,
        bucket_type="final_results",
    )
    simul_ht = hl.read_table(
        f"{simul_results_path}/merged.ht",
    )
    simul_ht = simul_ht.annotate(transcript=simul_ht.section.split("_")[0])
    simul_sections = simul_ht.aggregate(hl.agg.collect_as_set(simul_ht.section))

    # Get merged no-break table of round 1 of initial search
    ht = merge_round_no_break_ht(is_rescue=False, search_num=1)
    ht = ht.checkpoint(
        f"{TEMP_PATH_WITH_DEL}/tmp_single_rescue_transcripts.ht", overwrite=True
    )

    # Filter to transcripts with candidate breakpoints over the rescue threshold
    ht = ht.filter(ht.max_chisq >= rescue_threshold)

    # Annotate breakpoints of transcripts found with rescue threshold and checkpoint
    breakpoint_ht = ht.annotate(is_break=(ht.chisq == ht.max_chisq))
    breakpoint_ht = breakpoint_ht.filter(breakpoint_ht.is_break)
    breakpoint_ht = breakpoint_ht.annotate_globals(chisq_threshold=rescue_threshold)
    breakpoint_ht = breakpoint_ht.key_by("section")
    breakpoint_ht = breakpoint_ht.checkpoint(
        single_search_round_ht_path(
            is_rescue=True,
            search_num=1,
            is_break_found=True,
            is_breakpoint_only=True,
        ),
        overwrite=overwrite,
    )

    # Write break found HT for rescue search
    ht = ht.annotate(breakpoint=breakpoint_ht[ht.section].locus.position)
    ht = ht.filter(hl.is_defined(ht.breakpoint))
    ht = ht.checkpoint(
        single_search_round_ht_path(
            is_rescue=True,
            search_num=1,
            is_break_found=True,
            is_breakpoint_only=False,
        ),
        overwrite=overwrite,
    )
    rescue_single_sections = ht.aggregate(hl.agg.collect_as_set(ht.section))
    return simul_sections, rescue_single_sections


def get_rescue_2breaks_transcripts(
    overwrite: bool, initial_threshold: float = 6.6, rescue_threshold: float = 5.0
) -> hl.expr.SetExpression:
    """
    Get transcripts that have two simultaneous breakpoints above the 'rescue' threshold.

    These transcripts did not have a single break/simultaneous breaks over the initial search threshold
    but have simultaneous breakpoints that are above the lower 'rescue' search threshold.

    .. note::
        This function assumes that `process_sections` was run with `save_full_chisq_ht`
        set to True during the initial simultaneous breaks search (round 1).

    :param overwrite: Whether to overwrite output data if it exists.
    :param initial_threshold: Chi square significance threshold associated with
        the initial search pathway.
    :param rescue_threshold: Lower chi square significance threshold associated with
        the 'rescue' search pathway.
    :return: SetExpression of sections (transcript_start_stop)
        above 'rescue' threshold.
    """
    # Read in the transcripts found in the single break rescue search
    single_ht = hl.read_table(
        single_search_round_ht_path(
            is_rescue=True,
            search_num=1,
            is_break_found=True,
            is_breakpoint_only=True,
        )
    )
    single_break_rescue_sections = single_ht.aggregate(
        hl.agg.collect_as_set(single_ht.section)
    )

    # Merge all of the temporary chi square HTs saved in round 1 of
    # initial simul breaks search
    merge_simul_break_temp_hts(
        input_hts_path=SIMUL_BREAK_TEMP_PATH,
        batch_phrase="batch_temp_chisq",
        query_phrase="dataproc_temp_chisq",
        output_ht_path=f"{TEMP_PATH_WITH_DEL}/rescue_simul_chisq.ht",
        overwrite=overwrite,
    )
    ht = hl.read_table(f"{TEMP_PATH_WITH_DEL}/rescue_simul_chisq.ht")
    ht = ht.filter(
        (ht.max_chisq < initial_threshold)
        & (ht.max_chisq >= rescue_threshold)
        & ~hl.literal(single_break_rescue_sections).contains(ht.section)
    )
    ht = ht.drop("transcript")
    results_path = simul_search_round_bucket_path(
        is_rescue=True,
        search_num=1,
        bucket_type="final_results",
    )
    ht = ht.checkpoint(f"{results_path}/merged.ht", overwrite=overwrite)
    return ht.aggregate(hl.agg.collect_as_set(ht.section))


def get_rescue_transcripts_and_create_no_breaks_ht(overwrite: bool) -> None:
    """
    Get transcripts found in rescue pipeline and write final no breaks HT.

    Function gets transcripts found in rescue single and simultaneous breaks searches
    and creates final HT with all transcripts that do not have evidence of RMC.

    :param overwrite: Whether to overwrite output data if it exists.
    :return: None; function writes HT to resource path.
    """
    # Get the sections (transcript_start_stop) found in initial simultaneous breaks search,
    # rescue single break search, and rescue simultaneous breaks search
    init_simul_sections, rescue_single_sections = get_rescue_1break_transcripts(
        overwrite=overwrite
    )
    rescue_simul_sections = get_rescue_2breaks_transcripts(overwrite=overwrite)
    sections_with_breaks = init_simul_sections.union(rescue_simul_sections).union(
        rescue_single_sections
    )

    # Read in the no break found HT from initial search round 1
    ht = hl.read_table(
        single_search_round_ht_path(
            is_rescue=False,
            search_num=1,
            is_break_found=False,
            is_breakpoint_only=False,
        )
    )
    ht = ht.filter(~hl.literal(sections_with_breaks).contains(ht.section))
    no_break_sections = ht.aggregate(hl.agg.collect_as_set(ht.section))
    no_break_transcripts = hl.map(lambda x: x.split("_")[0], no_break_sections)
    hl.experimental.write_expression(
        no_break_transcripts, no_breaks, overwrite=overwrite
    )


def get_break_search_round_nums(
    rounds_path: str,
    round_num_regex: str = r"round(\d+)/$",
    google_project: str = "broad-mpg-gnomad",
) -> List[str]:
    r"""
    Get round numbers for a particular type of break search, e.g. single break in initial search.

    Function returns all round numbers for a particular type of break search
    by matching the round paths in a top-level bucket to a regex pattern.
    Regex matches from all capture groups and match instances in a given path are merged.
    Every round number, regardless of whether breaks were discovered, is returned.

    :param rounds_path: Path to top-level bucket containing break search round buckets.
    :param round_num_regex: Regex pattern to match the round number
        in a round bucket path. Default is r'round(\d+)/$'.
    :param google_project: Google project to use to read data from requester-pays buckets.
        Default is 'broad-mpg-gnomad'.
    :return: Sorted list of round numbers.
    """
    r = re.compile(round_num_regex)
    round_paths = (
        subprocess.check_output(
            ["gsutil", "-u", f"{google_project}", "ls", f"{rounds_path}"]
        )
        .decode("utf8")
        .strip()
        .split("\n")
    )
    round_nums = []
    for path in round_paths:
        m = r.findall(path)
        if len(m) > 0:
            # Merge all capture groups and match instances
            round_nums.append(int("".join("".join(t) for t in m)))
    return sorted(round_nums)


def check_break_search_round_nums(is_rescue: bool) -> List[int]:
    """
    Check for valid single and simultaneous break search round number outputs.

    Function checks that single and simultaneous search round numbers match, are
    increasing consecutive integers starting at 1, and that at least one search round was run.

    .. note::
        Assumes there is a folder for each search round run, regardless of whether there were breaks discovered

    :param is_rescue: Whether to operate on searches in rescue pathway.
    :return: Sorted list of round numbers.
    """
    single_search_round_nums = get_break_search_round_nums(
        single_search_bucket_path(is_rescue=is_rescue)
    )
    simul_search_round_nums = get_break_search_round_nums(
        simul_search_bucket_path(is_rescue=is_rescue)
    )
    logger.info(
        "Single search round numbers: %s\nSimultaneous search round numbers: %s",
        ",".join(map(str, single_search_round_nums)),
        ",".join(map(str, simul_search_round_nums)),
    )
    if single_search_round_nums != list(range(1, max(single_search_round_nums))):
        raise DataException(
            "Single search round numbers are not consecutive and increasing from 1, please double-check!"
        )
    if simul_search_round_nums != list(range(1, max(simul_search_round_nums))):
        raise DataException(
            "Simultaneous search round numbers are not consecutive and increasing from 1, please double-check!"
        )
    if single_search_round_nums != simul_search_round_nums:
        raise DataException(
            "Round numbers from single and simultaneous searches do not match, please double-check!"
        )
    if len(single_search_round_nums) == 0:
        raise DataException("No search rounds recorded, please double-check!")
    if len(single_search_round_nums) == 1:
        logger.warning(
            "Only one break search round recorded. Either no breaks were found or break search is not complete, please double-check!"
        )
    return single_search_round_nums


def merge_round_no_break_ht(is_rescue: bool, search_num: int) -> hl.Table:
    """
    Get merged single and simultaneous search no-break table from a given round of break search.

    :param is_rescue: Whether to operate on search in rescue pathway.
    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :return: Table of loci in sections where no breaks were found in the break search round.
    """
    single_no_break_path = single_search_round_ht_path(
        is_rescue=is_rescue,
        search_num=search_num,
        is_break_found=False,
        is_breakpoint_only=False,
    )
    if file_exists(single_no_break_path):
        raise DataException(
            f"No table found at {single_no_break_path}. Please double-check!"
        )
    ht = hl.read_table(single_no_break_path)
    simul_results_path = simul_search_round_bucket_path(
        is_rescue=is_rescue,
        search_num=search_num,
        bucket_type="final_results",
    )
    simul_no_break_path = f"{simul_results_path}/merged.ht"
    if file_exists(simul_no_break_path):
        simul_by_section_ht = hl.read_table(simul_no_break_path)
        ht = ht.annotate(breakpoints=simul_by_section_ht[ht.section].breakpoints)
        ht = ht.filter(hl.is_missing(ht.breakpoints)).drop("breakpoints")
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


def merge_rmc_hts(round_nums: List[int], is_rescue: bool) -> hl.Table:
    """
    Get table of final RMC sections in a given pathway (initial or rescue) after all break searches are complete.

    :param round_nums: List of round numbers to merge results across.
    :param is_rescue: Whether to operate on search in rescue pathway.
    :return: Table of final RMC sections. Schema:
        ----------------------------------------
        Row fields:
            'section_obs': int64
            'section_exp': float64
            'section_oe': float64
            'section_chisq': float64
            'search_type': str
            'transcript': str
            'interval': interval<locus<GRCh37>>
        ----------------------------------------
        Key: ['transcript', 'interval']
        ----------------------------------------
    """
    if len(round_nums) < 2:
        raise DataException(
            "At least two rounds of break search are needed if evidence of RMC is found, please double-check!"
        )

    logger.info(
        "Collecting and merging no-break HTs from each search round starting at round 2..."
    )
    hts = []
    for search_num in round_nums[1:]:
        # For each search round:
        # Get locus-level merged no-break table (no breaks in both single and simultaneous searches)
        ht = merge_round_no_break_ht(is_rescue=is_rescue, search_num=search_num)
        # Group to section-level and retain section obs- and exp-related annotations
        ht = ht.group_by("section").aggregate(
            section_obs=hl.agg.take(ht.section_obs, 1)[0],
            section_exp=hl.agg.take(ht.section_exp, 1)[0],
            section_oe=hl.agg.take(ht.section_oe, 1)[0],
        )
        hts.append(ht)
    rmc_ht = hts[0].union(*hts[1:])
    # Calculate chi-square value for each section
    rmc_ht = rmc_ht.annotate(
        section_chisq=calculate_section_chisq(rmc_ht.section_obs, rmc_ht.section_exp)
    )
    # Annotate search pathway type (initial or rescue)
    rmc_ht = rmc_ht.annotate(search_type="rescue" if is_rescue else "initial")

    # Convert section label to transcript and start and end positions
    rmc_ht = rmc_ht.key_by()
    rmc_ht = rmc_ht.transmute(
        transcript=rmc_ht.section.split("_")[0],
        start_pos=hl.int(rmc_ht.section.split("_")[1]),
        end_pos=hl.int(rmc_ht.section.split("_")[2]),
    )
    # Convert start and end positions to interval
    rmc_ht = rmc_ht.annotate(chr=gene_model.ht()[rmc_ht.transcript].chrom)
    rmc_ht = rmc_ht.transmute(
        interval=hl.parse_locus_interval(
            hl.format(
                "[%s:%s-%s]",
                rmc_ht.chr,
                rmc_ht.start_pos,
                rmc_ht.end_pos,
            )
        ),
    )
    rmc_ht = rmc_ht.key_by("transcript", "interval")
    return rmc_ht


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


def get_oe_bins(ht: hl.Table) -> None:
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
        - `section_exp`: Proportion of expected missense variants within transcript subsection
        - `section_oe`: Observed/expected missense variation ratio within transcript subsection

    :param hl.Table ht: Input Table containing all breaks results.
    :return: None; writes TSV with OE bins + annotations to `oe_bin_counts_tsv` resource path.
    """
    logger.info("Reading in ClinVar, de novo missense, and transcript HTs...")
    if not file_exists(clinvar_path_mis.path):
        import_clinvar_hi_variants(overwrite=True)
    if not file_exists(de_novo.path):
        import_de_novo_variants(overwrite=True)

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
    ht = ht.checkpoint(f"{TEMP_PATH}/breaks_oe_bin.ht", overwrite=True)

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
    group_ht = group_ht.checkpoint(f"{TEMP_PATH}/sections.ht", overwrite=True)
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
