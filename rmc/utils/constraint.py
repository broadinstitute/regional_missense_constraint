import logging
import re
import subprocess
from typing import Dict, List, Set, Union

import hail as hl
import scipy
from gnomad.resources.grch37.gnomad import coverage, public_release
from gnomad.resources.grch37.reference_data import vep_context
from gnomad.resources.resource_utils import DataException
from gnomad.utils.constraint import annotate_mutation_type
from gnomad.utils.file_utils import file_exists

from rmc.resources.basics import (
    SINGLE_BREAK_TEMP_PATH,
    TEMP_PATH,
    TEMP_PATH_WITH_FAST_DEL,
    TEMP_PATH_WITH_SLOW_DEL,
)
from rmc.resources.gnomad import constraint_ht, filtered_exomes
from rmc.resources.reference_data import (
    clinvar_plp_mis_haplo,
    gene_model,
    ndd_de_novo,
    transcript_cds,
    transcript_ref,
)
from rmc.resources.resource_utils import MISSENSE
from rmc.resources.rmc import (
    CONSTRAINT_ANNOTATIONS,
    CURRENT_FREEZE,
    FINAL_ANNOTATIONS,
    MIN_EXP_MIS,
    P_VALUE,
    SIMUL_SEARCH_ANNOTATIONS,
    constraint_prep,
    context_with_oe,
    context_with_oe_dedup,
    no_breaks_he_path,
    oe_bin_counts_tsv,
    rmc_browser,
    rmc_results,
    simul_search_bucket_path,
    simul_search_round_bucket_path,
    single_search_bucket_path,
    single_search_round_ht_path,
)
from rmc.utils.generic import (
    get_aa_from_context,
    get_annotations_from_context_ht_vep,
    get_constraint_transcripts,
    get_coverage_correction_expr,
    get_ref_aa,
    import_clinvar,
    import_de_novo_variants,
    keep_criteria,
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
    group_ht = annotate_mutation_type(group_ht)

    logger.info("Adjusting aggregated mutation rate with plateau model...")
    if locus_type == "X":
        model = group_ht.plateau_x_models["total"][group_ht.cpg]
    elif locus_type == "Y":
        model = group_ht.plateau_y_models["total"][group_ht.cpg]
    else:
        model = group_ht.plateau_models["total"][group_ht.cpg]

    group_ht = group_ht.annotate(mu_adj=group_ht.mu_agg * model[1] + model[0])

    logger.info(
        "Adjusting aggregated mutation rate with coverage correction to get expected"
        " counts..."
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
    Calculate probabilities (natural log) of the observed values under a Poisson model.

    Rate in model given by expected * section observed/expected values.

    :param cond_expr: Conditional expression to check before calculating probability.
    :param section_oe_expr: Expression of section observed/expected value.
    :param obs_expr: Expression containing observed variants count.
    :param exp_expr: Expression containing expected variants count.
    :return: natural log of the probability under Poisson model.
    """
    # log_p = True returns the natural logarithm of the probability density
    return hl.or_missing(
        cond_expr,
        hl.dpois(obs_expr, exp_expr * section_oe_expr, log_p=True),
    )


def get_max_chisq_per_group(
    ht: hl.Table,
    group_str: str,
    chisq_str: str,
    freeze: int,
) -> hl.Table:
    """
    Group input Table by given field and return maximum chi square value per group.

    'Group' in this context refers to either a transcript or transcript subsection.

    :param ht: Input Table.
    :param group_str: String of field containing transcript or transcript subsection information.
        Used to group observed and expected values.
    :param chisq_str: String of field containing chi square values to be checked.
    :param freeze: RMC data freeze number.
    :return: Table annotated with maximum chi square value per group
    """
    group_ht = ht.group_by(group_str).aggregate(
        # hl.agg.max ignores NaNs
        section_max_chisq=hl.agg.max(ht[chisq_str])
    )
    group_ht = group_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/freeze{freeze}_group_max_chisq.ht", overwrite=True
    )
    ht = ht.annotate(section_max_chisq=group_ht[ht.section].section_max_chisq)
    return ht


def search_for_break(
    ht: hl.Table,
    search_num: int,
    freeze: int,
    chisq_threshold: float = scipy.stats.chi2.ppf(1 - P_VALUE, 1),
    group_str: str = "section",
    min_num_exp_mis: float = MIN_EXP_MIS,
    save_chisq_ht: bool = False,
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
    :param freeze: RMC data freeze number.
    :param chisq_threshold: Chi-square significance threshold.
        Default is `scipy`.stats.chi2.ppf(1 - P_VALUE, 1)`.
        Default value used in ExAC was 10.8, which corresponds to a p-value of 0.001
        with 1 degree of freedom.
        (https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm)
    :param group_str: Field used to group Table observed and expected values. Default is 'section'.
    :param min_num_exp_mis: Minimum number of expected missense per transcript/transcript section.
        Sections that have fewer than this number of expected missense variants will not
        be computed (chi square will be annotated as a missing value).
        Default is MIN_EXP_MIS.
    :param save_chisq_ht: Whether to save HT with chi square values annotated for every locus.
        This saves a lot of extra data and should only occur once.
        Default is False.
    :return: Table annotated with whether position is a breakpoint (`is_break`).
    """
    logger.info(
        "Creating section null (no regional variability in missense depletion)       "
        " and alt (evidence of domains of missense constraint) expressions..."
    )
    logger.info(
        "Skipping breakpoints that create at least one section        that has < %i"
        " expected missense variants...",
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
    all_loci_chisq_ht_path = (
        f"{TEMP_PATH_WITH_FAST_DEL}/freeze{freeze}_round{search_num}_all_loci_chisq.ht"
    )
    if save_chisq_ht:
        all_loci_chisq_ht_path = f"{SINGLE_BREAK_TEMP_PATH}/all_loci_chisq.ht"
    ht = ht.checkpoint(all_loci_chisq_ht_path, overwrite=True)

    ht = get_max_chisq_per_group(ht, group_str, "chisq", freeze)
    return ht.annotate(
        is_break=((ht.chisq == ht.section_max_chisq) & (ht.chisq >= chisq_threshold))
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
    freeze: int,
    chisq_threshold: float = scipy.stats.chi2.ppf(1 - P_VALUE, 1),
    group_str: str = "section",
    save_chisq_ht: bool = False,
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
    :param freeze: RMC data freeze number.
    :param chisq_threshold: Chi-square significance threshold.
        Default is `scipy.stats.chi2.ppf(1 - P_VALUE, 1)`.
        Default value used in ExAC was 10.8, which corresponds to a p-value of 0.001
        with 1 degree of freedom.
        (https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm)
    :param group_str: Field used to group observed and expected values. Default is 'section'.
    :param save_chisq_ht: Whether to save HT with chi square values annotated for every locus.
        This saves a lot of extra data and should only occur during the initial search round.
        Default is False.
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
        freeze=freeze,
        chisq_threshold=chisq_threshold,
        save_chisq_ht=save_chisq_ht,
    )
    return ht


def create_no_breaks_he(freeze: int, overwrite: bool) -> None:
    """
    Write final no breaks HailExpression.

    :param freeze: RMC freeze number.
    :param overwrite: Whether to overwrite output data if it exists.
    :return: None; function writes HailExpression to resource path.
    """
    # Get the sections (transcript_start_stop) found in first round of simultaneous breaks search
    simul_results_path = simul_search_round_bucket_path(
        search_num=1,
        bucket_type="final_results",
        freeze=freeze,
    )
    simul_ht = hl.read_table(
        f"{simul_results_path}/merged.ht",
    )
    simul_sections = simul_ht.aggregate(hl.agg.collect_as_set(simul_ht.section))

    # Read in the no break found HT from the first round of single search
    ht = hl.read_table(
        single_search_round_ht_path(
            search_num=1,
            is_break_found=False,
            is_breakpoint_only=False,
            freeze=freeze,
        )
    )
    ht = ht.filter(~hl.literal(simul_sections).contains(ht.section))
    no_break_sections = ht.aggregate(hl.agg.collect_as_set(ht.section))
    logger.info(
        "%i transcripts did not have any evidence of RMC", len(no_break_sections)
    )
    no_break_transcripts = hl.map(lambda x: x.split("_")[0], no_break_sections)
    hl.experimental.write_expression(
        no_break_transcripts, no_breaks_he_path(freeze), overwrite=overwrite
    )


def merge_simul_break_temp_hts(
    input_hts_path: str,
    batch_phrase: str,
    query_phrase: str,
    output_ht_path: str,
    overwrite: bool,
    google_project: str = "broad-mpg-gnomad",
) -> None:
    """
    Read in simultaneous breaks temporary HTs at specified path and merge.

    .. note::
        - Assumes temp HTs are keyed by "section" or ("section", "i", "j")
        - Assumes temp HTs contain all annotations in SIMUL_BREAK_ANNOTATIONS

    :param input_hts_path: Path to bucket containing input temporary HTs.
    :param batch_phrase: String indicating name associated with HTs created
        using Hail Batch jobs, e.g. "under", or "batch_temp_chisq".
    :param query_phrase: String indicating name associated with HTs created
        using Hail Query jobs via Dataproc, e.g. "dataproc", or "dataproc_temp_chisq".
    :param output_ht_path: Desired path for output HT.
    :param overwrite: Whether to overwrite output HT if it exists.
    :param google_project: Google project to use to read data from requester-pays buckets.
        Default is 'broad-mpg-gnomad'.
    :return: None; function writes HT to specified path.
    """
    logger.info("Collecting all HT paths...")
    temp_ht_paths = (
        subprocess.check_output(
            ["gsutil", "-u", f"{google_project}", "ls", f"{input_hts_path}"]
        )
        .decode("utf8")
        .strip()
        .split("\n")
    )

    intermediate_hts = []
    ht_count = 0
    for ht_path in temp_ht_paths:
        ht_path = ht_path.strip("/")
        if (batch_phrase in ht_path or query_phrase in ht_path) and ht_path.endswith(
            "ht"
        ):
            ht_count += 1
            logger.info("Working on %s", ht_path)
            temp = hl.read_table(ht_path)
            if temp.count() > 0:
                # Tables containing transcripts/transcript sections that are over the transcript/section length threshold
                # are keyed by section, i, j
                # Tables containing transcripts/transcript sections that are under the length threshold are keyed
                # only by section
                # Rekey all tables here and select only the required fields to ensure the union on line 83 is able to work
                # This `key_by` should not shuffle because `section` is already the first key for both Tables
                temp = temp.key_by("section")
                row_fields = set(temp.row)
                if len(SIMUL_SEARCH_ANNOTATIONS.intersection(row_fields)) < len(
                    SIMUL_SEARCH_ANNOTATIONS
                ):
                    raise DataException(
                        "The following fields are missing from the temp table:"
                        f" {SIMUL_SEARCH_ANNOTATIONS.difference(row_fields)}!"
                    )
                temp = temp.select(*SIMUL_SEARCH_ANNOTATIONS)
                intermediate_hts.append(temp)

    logger.info("Found %i HTs and appended %i", ht_count, len(intermediate_hts))
    if len(intermediate_hts) == 0:
        raise DataException(
            "All temp tables had 0 rows. Please double check the temp tables!"
        )
    ht = intermediate_hts[0].union(*intermediate_hts[1:])
    ht.write(output_ht_path, overwrite=overwrite)


def get_break_search_round_nums(
    rounds_path: str,
    round_num_regex: str = r"round(\d+)/$",
    google_project: str = "broad-mpg-gnomad",
) -> List[int]:
    r"""
    Get round numbers for a particular type of break search, e.g. single break search.

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


def check_break_search_round_nums(freeze: int = CURRENT_FREEZE) -> List[int]:
    """
    Check for valid single and simultaneous break search round number outputs.

    Function checks that single and simultaneous search round numbers match, are
    increasing consecutive integers starting at 1, and that at least one search round was run.

    .. note::
        Assumes there is a folder for each search round run, regardless of whether there were breaks discovered

    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: Sorted list of round numbers.
    """
    # Get sorted round numbers
    single_search_round_nums = get_break_search_round_nums(
        single_search_bucket_path(freeze=freeze)
    )
    simul_search_round_nums = get_break_search_round_nums(
        simul_search_bucket_path(freeze=freeze)
    )
    logger.info(
        "Single search round numbers: %s\nSimultaneous search round numbers: %s",
        ",".join(map(str, single_search_round_nums)),
        ",".join(map(str, simul_search_round_nums)),
    )
    if len(single_search_round_nums) == 0 or len(simul_search_round_nums) == 0:
        raise DataException(
            "No rounds recorded for at least one of single and simultaneous search,"
            " please double-check!"
        )
    if single_search_round_nums != list(range(1, max(single_search_round_nums) + 1)):
        raise DataException(
            "Single search round numbers are not consecutive and increasing from 1,"
            " please double-check!"
        )
    if simul_search_round_nums != list(range(1, max(simul_search_round_nums) + 1)):
        raise DataException(
            "Simultaneous search round numbers are not consecutive and increasing from"
            " 1, please double-check!"
        )
    if single_search_round_nums != simul_search_round_nums:
        raise DataException(
            "Round numbers from single and simultaneous searches do not match, please"
            " double-check!"
        )
    if len(single_search_round_nums) == 1:
        logger.warning(
            "Only one break search round recorded. Either no breaks were found or break"
            " search is not complete, please double-check!"
        )
    return single_search_round_nums


def merge_round_no_break_ht(
    search_num: int,
    freeze: int,
    keep_annotations: Set[str] = CONSTRAINT_ANNOTATIONS,
) -> hl.Table:
    """
    Get merged single and simultaneous search no-break table from a given round of break search.

    Function starts with the round-specific single search no-breaks table
    and removes all sections in the round-specific simultaneous search breaks table.

    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param freeze: RMC data freeze number.
    :param keep_annotations: Fields to keep in the table. Default is `CONSTRAINT_ANNOTATIONS`.
    :return: Table of loci in sections where no breaks were found in the break search round. Schema:
        ----------------------------------------
        Row fields:
            'locus': locus<GRCh37>
            'section': str
            ...
        ----------------------------------------
        Key: ['locus', 'section']
        ----------------------------------------
        Note that there may be additional row fields depending on `keep_annotations`.
    """
    single_no_break_path = single_search_round_ht_path(
        search_num=search_num,
        is_break_found=False,
        is_breakpoint_only=False,
        freeze=freeze,
    )
    if not file_exists(single_no_break_path):
        raise DataException(
            f"No table found at {single_no_break_path}. Please double-check!"
        )
    ht = hl.read_table(single_no_break_path)

    # Check that all fields in `keep_annotations` are in the table
    row_fields = set(ht.row)
    if len(keep_annotations.intersection(row_fields)) < len(keep_annotations):
        raise DataException(
            "The following fields are missing from the break search result table:"
            f" {keep_annotations.difference(row_fields)}!"
        )
    ht = ht.select(*keep_annotations)

    simul_results_path = simul_search_round_bucket_path(
        search_num=search_num,
        bucket_type="final_results",
        freeze=freeze,
    )
    simul_break_path = f"{simul_results_path}/merged.ht"
    if file_exists(simul_break_path):
        simul_ht = hl.read_table(simul_break_path)
        simul_sections = simul_ht.aggregate(hl.agg.collect_as_set(simul_ht.section))
        ht = ht.filter(~hl.literal(simul_sections).contains(ht.section))
    else:
        logger.warning(
            "Simul breaks results HT (round %i) did not exist. Please double check that"
            " this was expected!",
            search_num,
        )
    return ht


def calculate_section_chisq(
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
) -> hl.expr.Float64Expression:
    """
    Check for significance of regional missense constraint within transcript section.

    Function calculates chi square expression that assess whether observed and expected
    missense counts in a given transcript session are significantly different than the null
    model (no evidence of regional missense constraint).

    Formula is: (section obs - section exp)^2 / section exp. Taken from ExAC RMC code.

    :param obs_expr: Total number of observed missense variants in section.
    :param exp_expr: Total number of expected missense variants in section.
    :return: Transcript section chi-squared value.
    """
    return ((obs_expr - exp_expr) ** 2) / exp_expr


def merge_rmc_hts(round_nums: List[int], freeze: int) -> hl.Table:
    """
    Get table of final RMC sections after all break searches are complete.

    :param round_nums: List of round numbers to merge results across.
    :param freeze: RMC data freeze number.
    :return: Table of final RMC sections. Schema:
        ----------------------------------------
        Row fields:
            'section_obs': int64
            'section_exp': float64
            'section_oe': float64
            'section_chisq': float64
            'transcript': str
            'interval': interval<locus<GRCh37>>
        ----------------------------------------
        Key: ['interval', 'transcript']
        ----------------------------------------
    """
    logger.warning(
        "This function performs a join followed by a rekey, which will trigger a"
        " shuffle!"
    )
    if len(round_nums) < 2:
        raise DataException(
            "At least two rounds of break search are needed if evidence of RMC is"
            " found, please double-check!"
        )

    logger.info(
        "Collecting and merging no-break HTs from each search round starting at round"
        " 2..."
    )
    hts = []
    for search_num in round_nums[1:]:
        # For each search round:
        # Get locus-level merged no-break table (no breaks in both single and simultaneous searches)
        ht = merge_round_no_break_ht(
            search_num=search_num,
            freeze=freeze,
            keep_annotations=FINAL_ANNOTATIONS,
        )
        # Group to section-level and retain section obs- and exp-related annotations
        ht = ht.group_by("section").aggregate(
            section_obs=hl.agg.take(ht.section_obs, 1)[0],
            section_exp=hl.agg.take(ht.section_exp, 1)[0],
            section_oe=hl.agg.take(ht.section_oe, 1)[0],
            chrom=hl.agg.take(ht.locus.contig, 1)[0],
        )
        hts.append(ht)
    rmc_ht = hts[0].union(*hts[1:])
    rmc_ht = rmc_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/freeze{freeze}_search_union.ht",
        overwrite=True,
    )
    # Calculate chi-square value for each section
    rmc_ht = rmc_ht.annotate(
        section_chisq=calculate_section_chisq(rmc_ht.section_obs, rmc_ht.section_exp)
    )

    # Convert section label to transcript and start and end positions
    rmc_ht = rmc_ht.key_by()
    rmc_ht = rmc_ht.transmute(
        transcript=rmc_ht.section.split("_")[0],
        start_pos=hl.int(rmc_ht.section.split("_")[1]),
        end_pos=hl.int(rmc_ht.section.split("_")[2]),
    )

    # TODO: Check that transcripts are fully covered
    # (Check that all section start and end positions
    # line up to cover transcript start/end positions)
    # TODO: Check that section start/ends reflect breakpoints found in searches

    # Convert start and end positions to interval
    rmc_ht = rmc_ht.transmute(
        interval=hl.parse_locus_interval(
            hl.format(
                "[%s:%s-%s]",
                rmc_ht.chrom,
                rmc_ht.start_pos,
                rmc_ht.end_pos,
            )
        ),
    )
    rmc_ht = rmc_ht.key_by("interval", "transcript")
    return rmc_ht


def get_oe_annotation(ht: hl.Table, freeze: int) -> hl.Table:
    """
    Annotate input Table with observed to expected missense (OE) ratio per transcript.

    Use regional OE value if available, otherwise use transcript OE value.

    .. note::
        - Assumes input Table has `locus` and `trancript` annotations
        - OE values are transcript specific
        - Assumes merged RMC results HT exists
        - Assumes merged RMC results HT is annotated per transcript section with:
            - `section_oe`: Missense observed/expected ratio
            - `interval`: Transcript section start position to end position

    :param hl.Table ht: Input Table.
    :param freeze: RMC data freeze number.
    :return: Table with `oe` annotation.
    """
    overall_oe_ht = (
        constraint_prep.ht().select_globals().select("total_obs", "total_exp")
    )
    group_ht = overall_oe_ht.group_by("transcript").aggregate(
        obs=hl.agg.take(overall_oe_ht.total_obs, 1)[0],
        exp=hl.agg.take(overall_oe_ht.total_exp, 1)[0],
    )
    # Recalculating transcript level OE ratio because previous OE ratio (`overall_oe`)
    # is capped at 1 for regional missense constraint calculation purposes
    group_ht = group_ht.annotate(transcript_oe=group_ht.obs / group_ht.exp)

    # Read in LoF constraint HT to get OE ratio for five transcripts missing in v2 RMC results
    # # 'ENST00000304270', 'ENST00000344415', 'ENST00000373521', 'ENST00000381708', 'ENST00000596936'
    # All 5 of these transcripts have extremely low coverage in gnomAD
    # Will keep for consistency with v2 LoF results but they look terrible
    # NOTE: LoF HT is keyed by gene and transcript, but `_key_by_assert_sorted` doesn't work here for v2 version
    # Throws this error: hail.utils.java.FatalError: IllegalArgumentException
    lof_ht = constraint_ht.ht().select("oe_mis").key_by("transcript")
    ht = ht.annotate(
        gnomad_transcript_oe=lof_ht[ht.transcript].oe_mis,
        rmc_transcript_oe=group_ht[ht.transcript].transcript_oe,
    )
    ht = ht.transmute(
        transcript_oe=hl.coalesce(ht.rmc_transcript_oe, ht.gnomad_transcript_oe)
    )

    if not file_exists(rmc_results.versions[freeze].path):
        raise DataException("Merged RMC results table does not exist!")
    rmc_ht = rmc_results.versions[freeze].ht().key_by("interval")
    ht = ht.annotate(
        section_oe=rmc_ht.index(ht.locus, all_matches=True)
        .filter(lambda x: x.transcript == ht.transcript)
        .section_oe
    )
    ht = ht.annotate(
        section_oe=hl.or_missing(
            hl.len(ht.section_oe) > 0,
            ht.section_oe[0],
        ),
    )
    return ht.transmute(oe=hl.coalesce(ht.section_oe, ht.transcript_oe))


def create_context_with_oe(
    freeze: int,
    missense_str: str = MISSENSE,
    n_partitions: int = 30000,
    overwrite_temp: bool = False,
) -> None:
    """
    Filter VEP context Table to missense variants in canonical transcripts, and add missense observed/expected.

    Function writes two Tables:
        - `context_with_oe`: Table of all missense variants in canonical transcripts + all annotations
            (includes duplicate variants, i.e. those in multiple transcripts).
            Annotations are: transcript, most severe consequence, codons, reference and alternate amino acids,
            Polyphen-2, and SIFT.
        - `context_with_oe_dedup`: Deduplicated version of `context_with_oe` that only contains missense o/e and transcript annotations.

    :param freeze: RMC data freeze number.
    :param str missense_str: String representing missense variant consequence. Default is MISSENSE.
    :param int n_partitions: Number of desired partitions for the VEP context Table.
        Repartition VEP context Table to this number on read.
        Default is 30000.
    :param bool overwrite_temp: Whether to overwrite intermediate temporary (OE-independent) data if it already exists.
        If False, will read existing intermediate temporary data rather than overwriting.
        Default is False.
    :return: None; function writes Table to resource path.
    """
    logger.info("Importing set of transcripts to keep...")
    transcripts = get_constraint_transcripts(outlier=False)

    logger.info(
        "Reading in VEP context HT and filtering to missense variants in canonical"
        " transcripts..."
    )
    # Using vep_context.path to read in table with fewer partitions
    # VEP context resource has 62164 partitions
    ht = (
        hl.read_table(vep_context.path, _n_partitions=n_partitions)
        .select_globals()
        .select("vep", "was_split")
    )
    ht = process_vep(ht, filter_csq=True, csq=missense_str)
    ht = ht.filter(transcripts.contains(ht.transcript_consequences.transcript_id))
    ht = get_annotations_from_context_ht_vep(ht)
    # Save context Table to temporary path with specified deletion policy because this is a very large file
    # and relevant information will be saved at `context_with_oe` (written below)
    # Checkpointing here to force this expensive computation to compute
    # before joining with RMC tables in `get_oe_annotation`
    # This computation is expensive, so overwrite only if specified (otherwise, read existing file)
    ht = ht.checkpoint(
        f"{TEMP_PATH_WITH_SLOW_DEL}/vep_context_mis_only_annot.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    logger.info(
        "Adding regional missense constraint missense o/e annotation and writing to"
        " resource path..."
    )
    ht = get_oe_annotation(ht, freeze)
    ht = ht.key_by("locus", "alleles", "transcript")
    ht = ht.checkpoint(
        context_with_oe.versions[freeze].path,
        overwrite=True,
    )
    logger.info("Output OE-annotated context HT fields: %s", set(ht.row))

    logger.info(
        "Creating dedup context with oe (with oe and transcript annotations only)..."
    )
    ht = ht.select("oe").key_by("locus", "alleles")
    ht = ht.collect_by_key()
    ht = ht.annotate(**{"oe": hl.nanmin(ht.values.oe)})
    ht = ht.annotate(
        **{"transcript": ht.values.find(lambda x: x.oe == ht.oe).transcript}
    )
    ht = ht.drop("values")
    ht.write(
        context_with_oe_dedup.versions[freeze].path,
        overwrite=True,
    )
    logger.info("Output OE-annotated dedup context HT fields: %s", set(ht.row))


def annot_rmc_with_aa(ht: hl.Table, overwrite_temp: bool):
    """
    Annotate RMC results HT with amino acid information.

    :param ht: Input RMC results HT.
    :param overwrite_temp: Whether to overwrite temporary data.
        If False, will read existing temp data rather than overwriting.
        If True, will overwrite temp data.
    :return: RMC results HT annotated with amino acid information.
    """
    logger.info("Getting amino acid information from context HT...")
    rmc_transcripts = ht.aggregate(hl.agg.collect_as_set(ht.transcript))
    # Get amino acid information from context table for variants in chosen transcripts
    context_ht = get_aa_from_context(
        overwrite_temp=overwrite_temp, keep_transcripts=rmc_transcripts
    )
    # Get reference AA label for each locus-transcript combination in context HT
    context_ht = get_ref_aa(
        ht=context_ht,
        overwrite_temp=overwrite_temp,
    )
    # NOTE: For transcripts on the negative strand, `start_aa` is the larger AA number and
    # `stop_aa` is the smaller number
    ht = ht.annotate(
        start_aa=context_ht[ht.start_coordinate, ht.transcript].ref_aa,
        stop_aa=context_ht[ht.stop_coordinate, ht.transcript].ref_aa,
    )
    ht = ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/rmc_results_aa_annot.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )
    return check_and_fix_missing_aa(ht, context_ht, overwrite_temp)


def join_and_fix_aa(
    ht: hl.Table, fix_ht: hl.Table, fix_ts_boundaries: bool
) -> hl.Table:
    """
    Add start and stop amino acid information from fix_ht to ht.

    :param ht: Input HT. Should be RMC results HT annotated with amino acid information.
    :param fix_ht: HT start and stop amino acids only for RMC regions that were
        previously missing amino acid information.
    :param fix_ts_boundaries: Whether to fix missing amino acid annotations at transcript boundaries
        (CDS starts/stops) or at RMC region start/stops. If True, adjusts AAs at transcript boundaries only.
    :return: RMC results HT annotated with amino acid information from fix HT.
    """
    if fix_ts_boundaries:
        return ht.annotate(
            start_aa=hl.if_else(
                hl.is_missing(ht.start_aa) & ht.is_transcript_start,
                fix_ht[ht.interval, ht.transcript].start_aa,
                ht.start_aa,
            ),
            stop_aa=hl.if_else(
                hl.is_missing(ht.stop_aa) & ht.is_transcript_stop,
                fix_ht[ht.interval, ht.transcript].stop_aa,
                ht.stop_aa,
            ),
        )
    return ht.annotate(
        start_aa=hl.if_else(
            hl.is_missing(ht.start_aa),
            fix_ht[ht.interval, ht.transcript].start_aa,
            ht.start_aa,
        ),
        stop_aa=hl.if_else(
            hl.is_missing(ht.stop_aa),
            fix_ht[ht.interval, ht.transcript].stop_aa,
            ht.stop_aa,
        ),
    )


def fix_transcript_start_stop_aas(
    ht: hl.Table, context_ht: hl.Table, overwrite_temp: bool
) -> hl.Table:
    """
    Fix missing start or stop amino acid at transcript start or stop coordinates.

    Fix for transcript starts missing amino acid annotations:
        - Met1 (if + strand)
        - Closest amino acid (if - strand)

    Fix for transcript ends missing amino acid annotations:
        - Closest amino acid (if + strand)
        - Met1 (if - strand)

    :param ht: Input HT. Should be RMC results HT annotated with amino acid information.
    :param context_ht: VEP context HT filtered to keep only transcript ID, protein number, and amino acid information.
        Keys are ["locus", "transcript"].
    :param overwrite_temp: Whether to overwrite temporary data.
        If False, will read existing temp data rather than overwriting.
        If True, will overwrite temp data.
    :return: HT containing start and stop amino acids only for RMC regions that were
        previously missing amino acid information at transcript stops/ends.
    """
    # Filter to regions missing AA at transcript starts/ends
    miss_start_stop_ht = ht.filter(
        (hl.is_missing(ht.start_aa) & ht.is_transcript_start)
        | (hl.is_missing(ht.stop_aa) & ht.is_transcript_stop)
    )
    miss_start_stop_ht = miss_start_stop_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/transcript_start_stop_missing_aa.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    # Get all positions with defined amino acids in context HT
    miss_start_stop_transcripts = miss_start_stop_ht.aggregate(
        hl.agg.collect_as_set(miss_start_stop_ht.transcript)
    )
    context_ht = context_ht.filter(
        hl.literal(miss_start_stop_transcripts).contains(context_ht.transcript)
    )
    # Keep only key fields: locus and transcript
    pos_ht = context_ht.select()
    pos_ht = pos_ht.group_by("transcript").aggregate(
        positions=hl.sorted(hl.agg.collect(pos_ht.locus.position))
    )
    pos_ht = pos_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/pos_per_transcript.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )
    miss_start_stop_ht = miss_start_stop_ht.annotate(
        # Find the closest position with defined AA to the start coordinate
        closest_start_pos=hl.bind(
            lambda x: x[
                hl.binary_search(x, miss_start_stop_ht.start_coordinate.position)
            ],
            pos_ht[miss_start_stop_ht.transcript].positions,
        ),
        # Find the closest position with defined AA to the stop coordinate
        closest_stop_pos=hl.bind(
            lambda x: x[
                hl.binary_search(x, miss_start_stop_ht.stop_coordinate.position)
            ],
            pos_ht[miss_start_stop_ht.transcript].positions,
        ),
    )
    # Checkpoint to make sure the binary search/join only runs once,
    # especially since there is another join immediately afterwards
    miss_start_stop_ht = miss_start_stop_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/ts_missing_start_stop_pos.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    # Get amino acid annotations
    miss_start_stop_ht = miss_start_stop_ht.annotate(
        start_aa=hl.if_else(
            hl.is_missing(miss_start_stop_ht.start_aa),
            hl.if_else(
                # Set all missing + strand transcript starts to Met1
                miss_start_stop_ht.strand == "+",
                "Met1",
                context_ht[
                    hl.locus(
                        miss_start_stop_ht.start_coordinate.contig,
                        miss_start_stop_ht.closest_start_pos,
                    ),
                    miss_start_stop_ht.transcript,
                ].ref_aa,
            ),
            miss_start_stop_ht.start_aa,
        ),
        stop_aa=hl.if_else(
            hl.is_missing(miss_start_stop_ht.stop_aa),
            hl.if_else(
                # Set all missing - strand transcript stops to Met1
                miss_start_stop_ht.strand == "-",
                "Met1",
                context_ht[
                    hl.locus(
                        miss_start_stop_ht.stop_coordinate.contig,
                        miss_start_stop_ht.closest_stop_pos,
                    ),
                    miss_start_stop_ht.transcript,
                ].ref_aa,
            ),
            miss_start_stop_ht.stop_aa,
        ),
    )
    return miss_start_stop_ht.select("start_aa", "stop_aa")


def fix_region_start_stop_aas(
    ht: hl.Table, context_ht: hl.Table, overwrite_temp: bool
) -> hl.Table:
    """
    Fix missing start or stop amino acid information for non-transcript starts and stops.

    .. note::
        - This function assumes that all start coordinates missing AA annotations are one position larger than
            the previous exon stop position. (This was the case in RMC freeze 5.)
        - This function assumes that all stop coordinates missing AA annotations are one position smaller than
            the following exon start position. (This was the case in RMC freeze 5).
        - See this ticket for more information about RMC freeze 5:
            https://github.com/broadinstitute/rmc_production/issues/149

    Fix for region start coordinates missing AA annotations:
        - Function adds amino acid from following exon start to fix missing AA annotation

    Fix for region stop coordinates missing AA annotations:
        - Function adds amino acid from previous exon stop to fix missing AA annotation

    :param ht: Input HT. Should be RMC results HT annotated with amino acid information.
    :param context_ht: VEP context HT filtered to keep only transcript ID, protein number, and amino acid information.
    :param overwrite_temp: Whether to overwrite temporary data.
        If False, will read existing temp data rather than overwriting.
        If True, will overwrite temp data.
    :return: HT containing start and stop amino acids only for RMC regions that were
        previously missing amino acid information.
    """
    missing_ht = ht.filter(hl.is_missing(ht.start_aa) | hl.is_missing(ht.stop_aa))
    missing_ht = missing_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/region_start_stop_missing_aa.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )
    missing_transcripts = hl.literal(
        missing_ht.aggregate(hl.agg.collect_as_set(missing_ht.transcript))
    )

    # Read in CDS HT resource
    # Filter CDS and context HT to transcripts with internal missing start or stop AAs
    context_ht = context_ht.filter(missing_transcripts.contains(context_ht.transcript))
    cds_ht = transcript_cds.ht()
    cds_ht = cds_ht.filter(missing_transcripts.contains(cds_ht.transcript))
    cds_ht = cds_ht.annotate(
        exon_start=cds_ht.interval.start,
        exon_end=cds_ht.interval.end,
    )

    # Create two versions of CDS HT to fix missing region start and stop AAs
    def _create_exon_position_ht(
        cds_ht: hl.Table, fix_missing_start: bool, overwrite_temp: bool
    ) -> hl.Table:
        """
        Create Table keyed by either exon start or stop position.

        :param cds_ht: CDS HT resource filtered to transcripts missing region start or stop
            AA annotations only.
        :param fix_missing_start: Whether the output HT is designed to fix region starts
            that are missing AA annotations.
        :param overwrite_temp: Whether to overwrite temporary data.
            If False, will read existing temp data rather than overwriting.
            If True, will overwrite temp data.
        :return: CDS HT keyed by either exon start or stop position and transcript.
        """
        if fix_missing_start:
            # Region starts missing AAs need to get adjusted to use the AA from the previous
            # exon start
            # Reorder HT so that _prev_nonnull returns the next exon start
            checkpoint_path = f"{TEMP_PATH_WITH_FAST_DEL}/cds_start_fix.ht"
            cds_ht = cds_ht.order_by(
                hl.asc(cds_ht.transcript), hl.desc(cds_ht.exon_start)
            )
            cds_ht = cds_ht.annotate(
                next_exon_start=hl.scan._prev_nonnull(cds_ht.exon_start)
            )
            cds_ht = cds_ht.key_by(cds_ht.exon_stop, cds_ht.transcript)
        else:
            checkpoint_path = f"{TEMP_PATH_WITH_FAST_DEL}/cds_stop_fix.ht"
            cds_ht = cds_ht.annotate(
                prev_exon_stop=hl.scan._prev_nonnull(cds_ht.exon_stop)
            )
            cds_ht = cds_ht.key_by(cds_ht.exon_start, cds_ht.transcript)

        cds_ht = cds_ht.checkpoint(
            checkpoint_path,
            _read_if_exists=not overwrite_temp,
            overwrite=overwrite_temp,
        )
        return cds_ht

    # Get relevant exon start and stop positions
    start_fix_ht = _create_exon_position_ht(
        cds_ht=cds_ht, fix_missing_start=True, overwrite_temp=overwrite_temp
    )
    stop_fix_ht = _create_exon_position_ht(
        cds_ht=cds_ht, fix_missing_start=False, overwrite_temp=overwrite_temp
    )
    missing_ht = missing_ht.annotate(
        next_exon_start=hl.or_missing(
            hl.is_missing(missing_ht.start_aa),
            stop_fix_ht[
                hl.locus(
                    missing_ht.start_coordinate.contig,
                    missing_ht.start_coordinate.locus - 1,
                ),
                missing_ht.transcript,
            ].next_exon_start,
        ),
        prev_exon_stop=hl.or_missing(
            hl.is_missing(missing_ht.stop_aa),
            start_fix_ht[
                hl.locus(
                    missing_ht.stop_coordinate.contig,
                    missing_ht.stop_coordinate.locus + 1,
                ),
                missing_ht.transcript,
            ].prev_exon_stop,
        ),
    )
    missing_ht = missing_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/region_start_stop_missing_aa_exon_annot.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )
    missing_ht = missing_ht.annotate(
        start_aa=hl.if_else(
            hl.is_missing(missing_ht.start_aa),
            context_ht[missing_ht.next_exon_start, missing_ht.transcript].ref_aa,
            missing_ht.start_aa,
        ),
        stop_aa=hl.if_else(
            hl.is_missing(missing_ht.stop_aa),
            context_ht[missing_ht.prev_exon_stop, missing_ht.transcript].ref_aa,
            missing_ht.stop_aa,
        ),
    )
    return missing_ht.select("start_aa", "stop_aa")


def check_and_fix_missing_aa(
    ht: hl.Table, context_ht: hl.Table, overwrite_temp: bool
) -> hl.Table:
    """
    Check for and fix any missing amino acid information in RMC results HT.

    :param ht: Input Table. Should be RMC results HT annotated with amino acid information.
    :param context_ht: VEP context HT filtered to keep only transcript ID, protein number, and amino acid information.
    :param overwrite_temp: Whether to overwrite temporary data.
        If False, will read existing temp data rather than overwriting.
        If True, will overwrite temp data.
    :return: RMC results HT annotated with amino acids for all region start/stops.
    """
    missing_ht = ht.filter(hl.is_missing(ht.start_aa) | hl.is_missing(ht.stop_aa))
    missing_start_or_stop = missing_ht.count()
    if missing_start_or_stop == 0:
        return ht

    logger.info(
        "%i regions are missing either their stop or start AA!", missing_start_or_stop
    )

    # Get CDS start/stops from CDS HT
    transcript_ht = transcript_ref.ht().key_by()

    # Annotate whether RMC region is at the beginning or end of the transcript in coordinate space
    # Also add strand annotation from `transcript_ref`
    ht = ht.annotate(**transcript_ht[ht.transcript])
    ht = ht.annotate(
        # NOTE: This is actually the transcript end for transcripts on the negative strand
        is_transcript_start=ht.start_coordinate == ht.gnomad_cds_start,
        is_transcript_stop=ht.stop_coordinate == ht.gnomad_cds_end,
    )
    transcript_start_stop_fix_ht = fix_transcript_start_stop_aas(
        ht, context_ht, overwrite_temp
    )
    ht = join_and_fix_aa(ht, transcript_start_stop_fix_ht, fix_ts_boundaries=True)
    ht = ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/rmc_results_transcript_start_stop_aa_fix.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )
    region_fix_ht = fix_region_start_stop_aas(ht, context_ht, overwrite_temp)
    ht = join_and_fix_aa(ht, region_fix_ht, fix_ts_boundaries=False)
    ht = ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/rmc_results_all_aa_fix.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )
    return ht


def add_globals_rmc_browser(ht: hl.Table) -> hl.Table:
    """
    Annotate HT globals with RMC transcript information.

    Function is used when reformatting RMC results for browser release.
    Annotates:
        - transcripts not searched for RMC
        - transcripts without evidence of RMC
        - outlier transcripts
    :param HT: Input Table. Should be RMC results HT.
    :return: RMC results HT with updated globals annotations.
    """
    # Get transcripts with evidence of RMC
    rmc_transcripts = hl.literal(ht.aggregate(hl.agg.collect_as_set(ht.transcript)))

    # Get all QC pass transcripts and outlier transcripts
    qc_pass_transcripts = get_constraint_transcripts(outlier=False)
    outlier_transcripts = get_constraint_transcripts(outlier=True)

    # Get all transcripts displayed on browser
    transcript_ht = gene_model.ht()
    all_transcripts = transcript_ht.aggregate(
        hl.agg.collect_as_set(transcript_ht.transcript)
    )
    transcripts_no_rmc = qc_pass_transcripts.difference(rmc_transcripts)
    transcripts_not_searched = all_transcripts.difference(qc_pass_transcripts)
    ht = ht.select_globals()
    return ht.annotate_globals(
        transcripts_not_searched=transcripts_not_searched,
        transcripts_no_rmc=transcripts_no_rmc,
        outlier_transcripts=outlier_transcripts,
    )


def format_rmc_browser_ht(freeze: int, overwrite_temp: bool) -> None:
    """
    Reformat annotations in input HT for release.

    This reformatting is necessary to load data into the gnomAD browser.

    Desired schema:
    ----------------------------------------
    Global fields:
        'transcripts_not_searched': set<str>
        'transcripts_no_rmc': set<str>
        'outlier_transcripts': set<str>
    ----------------------------------------
    Row fields:
        'transcript': str
        'regions': array<struct {
            start_coordinate: locus<GRCh37>,
            stop_coordinate: locus<GRCh37>,
            start_aa: str,
            stop_aa: str,
            obs: int64,
            exp: float64,
            oe: float64,
            chisq: float64
        }>
    ----------------------------------------
    Key: ['transcript']
    ----------------------------------------

    :param freeze: RMC data freeze number.
    :param overwrite_temp: Whether to overwrite temporary data.
        If False, will read existing temp data rather than overwriting.
        If True, will overwrite temp data.
    :return: None; writes Table with desired schema to resource path.
    """
    ht = rmc_results.versions[freeze].ht()

    # Annotate start and stop coordinates per region
    ht = ht.annotate(
        start_coordinate=ht.interval.start,
        stop_coordinate=ht.interval.end,
    )

    # Annotate amino acids
    ht = annot_rmc_with_aa(ht, overwrite_temp)

    # Add region struct
    ht = ht.annotate(
        region=hl.struct(
            start_coordinate=ht.start_coordinate,
            stop_coordinate=ht.stop_coordinate,
            start_aa=ht.start_aa,
            stop_aa=ht.stop_aa,
            obs=ht.section_obs,
            exp=ht.section_exp,
            oe=ht.section_oe,
            chisq=ht.section_chisq,
        )
    )
    # Group Table by transcript
    ht = ht.group_by("transcript").aggregate(regions=hl.agg.collect(ht.region))

    # Annotate globals and write
    ht = add_globals_rmc_browser(ht)
    ht.write(rmc_browser.versions[freeze].path, overwrite=True)


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
    if not file_exists(clinvar_plp_mis_haplo.path):
        import_clinvar(overwrite=True)
    if not file_exists(ndd_de_novo.path):
        import_de_novo_variants(overwrite=True)

    clinvar_ht = clinvar_plp_mis_haplo.ht()
    dn_ht = ndd_de_novo.ht()
    transcript_ht = transcript_ref.ht()

    # Split de novo HT into two HTs -- one for controls and one for cases
    dn_controls_ht = dn_ht.filter(dn_ht.case_control == "control")
    dn_case_ht = dn_ht.filter(dn_ht.case_control != "control")

    # Get total number of coding base pairs, also ClinVar and DNM variants
    # TODO: use exon field in gene model HT to get only coding bases (rather than using transcript end - start)
    transcript_ht = transcript_ht.annotate(
        bp=transcript_ht.gnomad_cds_end - transcript_ht.gnomad_cds_start
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
            f"Expected 5 OE bins but found {assess_ht_count}. Please double check and"
            " rerun!"
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
