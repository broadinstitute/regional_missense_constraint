import logging
import re
import subprocess
from typing import List, Set

import hail as hl
import scipy
from gnomad.resources.grch37.gnomad import coverage, public_release
from gnomad.resources.grch37.reference_data import vep_context
from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists, check_file_exists_raise_error

from rmc.resources.basics import (
    SINGLE_BREAK_TEMP_PATH,
    TEMP_PATH,
    TEMP_PATH_WITH_FAST_DEL,
    TEMP_PATH_WITH_SLOW_DEL,
)
from rmc.resources.gnomad import constraint_ht, prop_obs_coverage
from rmc.resources.reference_data import clinvar_plp_mis_haplo, gene_model, ndd_de_novo
from rmc.resources.resource_utils import KEEP_CODING_CSQ, MISSENSE
from rmc.resources.rmc import (
    CURRENT_FREEZE,
    FINAL_ANNOTATIONS,
    MIN_EXP_MIS,
    P_VALUE,
    SIMUL_SEARCH_ANNOTATIONS,
    constraint_prep,
    context_with_oe,
    context_with_oe_dedup,
    coverage_plateau_models_path,
    filtered_context,
    no_breaks_he_path,
    oe_bin_counts_tsv,
    rmc_results,
    simul_search_bucket_path,
    simul_search_round_bucket_path,
    single_search_bucket_path,
    single_search_round_ht_path,
)
from rmc.utils.generic import (
    filter_context_using_gnomad,
    filter_to_region_type,
    generate_models,
    get_annotations_from_context_ht_vep,
    get_constraint_transcripts,
    get_coverage_correction_expr,
    import_clinvar,
    import_de_novo_variants,
    keep_criteria,
    process_context_ht,
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
    "coverage",
]
"""
Fields to group by when calculating expected variants per variant type.

Fields based off of gnomAD LoF repo.
"""


def add_obs_annotation(
    ht: hl.Table,
    gnomad_data_type: str = "exomes",
) -> hl.Table:
    """
    Add observed variant count for each variant in input Table.

    Check if locus/allele are present in gnomAD and add as annotation.

    :param ht: Input Table.
    :param gnomad_data_type: gnomAD data type. Used to retrieve public release and coverage resources.
        Must be one of 'exomes' or 'genomes' (check is done within `public_release`).
        Default is 'exomes'.
    :return: Table with observed variant annotation.
    """
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


def get_fwd_cumulative_count_expr(
    group_expr: hl.expr.StringExpression,
    count_expr: hl.expr.Int64Expression,
) -> hl.expr.DictExpression:
    """
    Return annotation with the cumulative number of variant counts (non-inclusive).

    Function can return either the cumulative expected and observed counts.

    Value returned is non-inclusive (does not include value of row) due to the nature of `hl.scan`
    and needs to be corrected later.

    This function can produce the scan when searching for the first break or when searching for additional break(s).

    :param group_expr: Transcript or transcript subsection expression. Used to group count values.
    :param count_expr: Variant count expression.
    :return: DictExpression containing scan expressions for cumulative variant counts for `search_expr`.
        Keys of DictExpression are the elements in `group_expr`.
    """
    return hl.scan.group_by(group_expr, hl.scan.sum(count_expr))


def adjust_fwd_cumulative_count_expr(
    cumulative_count_expr: hl.expr.DictExpression,
    count_expr: hl.expr.Int64Expression,
    group_expr: hl.expr.StringExpression,
) -> hl.expr.Int64Expression:
    """
    Adjust the scan with the cumulative number of variant counts.

    This adjustment is necessary because scans are always one line behind, and we want the values to match per line.

    This function can correct the scan created when searching for the first break or when searching for additional break(s).

    .. note::
        This function expects that `cumulative_count_expr` is a DictExpression keyed by elements in `group_expr`.

    :param cumulative_count_expr: DictExpression containing scan expression with cumulative variant counts per site.
    :param count_expr: IntExpression with variant counts at site.
    :param group_expr: StringExpression containing transcript or transcript subsection information.
    :return: Adjusted cumulative variant counts expression.
    """
    return hl.if_else(
        # Check if the current transcript/section exists in the cumulative count dictionary
        # If it doesn't exist, that means this is the first line in the HT for that particular transcript
        # The first line of a scan is always missing, but we want it to exist
        # Thus, set the cumulative count equal to the current observed value
        hl.is_missing(cumulative_count_expr.get(group_expr)),
        count_expr,
        # Otherwise, add the current count to the scan to make sure the cumulative value isn't one line behind
        cumulative_count_expr[group_expr] + count_expr,
    )


def calculate_exp_from_mu(
    context_ht: hl.Table,
    locus_type: str,
    groupings: List[str] = GROUPINGS,
) -> hl.Table:
    """
    Annotate Table with the per-variant (locus-allele) expected counts based on the per-variant mu.

    .. note::
        - Assumes that input Table is annotated with all of the fields in `groupings` and that
            the names match exactly.
        - Assumes that input Table is annotated with `cpg` and `mu_snp` (raw mutation rate probability
            without coverage correction).
        - Assumes that input Table is filtered to autosomes/PAR only, X nonPAR only, or Y nonPAR only.
        - Assumes that input Table contains coverage and plateau models in its global annotations
            (`coverage_model`, `plateau_models`).
        - Adds `expected` and `coverage_correction` annotations.

    :param context_ht: Variant-level input context Table.
    :param locus_type: Locus type of input Table. One of "X", "Y", or "autosomes".
        NOTE: Will treat any input other than "X" or "Y" as autosomes.
    :param groupings: List of Table fields used to group Table to adjust mutation rate.
        Table must be annotated with these fields. Default is `GROUPINGS`.
    :return: Table annotated with per-variant expected counts and coverage correction.
    """
    logger.info(
        "Grouping by %s and aggregating mutation rates within groupings...", groupings
    )
    group_ht = context_ht.group_by(*groupings).aggregate(
        mu_agg=hl.agg.sum(context_ht.mu_snp),
        n_grouping=hl.agg.count(),
        # `cpg` is a function of `context`, `ref``, `alt` which are part of `groupings`
        # so will be the same for all alleles in a grouping
        cpg=hl.agg.take(context_ht.cpg, 1)[0],
    )

    logger.info("Adjusting aggregated mutation rates with plateau model...")
    if locus_type == "X":
        model = group_ht.plateau_x_models["total"][group_ht.cpg]
    elif locus_type == "Y":
        model = group_ht.plateau_y_models["total"][group_ht.cpg]
    else:
        model = group_ht.plateau_models["total"][group_ht.cpg]

    group_ht = group_ht.annotate(mu_adj=group_ht.mu_agg * model[1] + model[0])

    logger.info(
        "Adjusting aggregated mutation rates with coverage correction to get expected"
        " counts for each grouping..."
    )
    group_ht = group_ht.annotate(
        coverage_correction=get_coverage_correction_expr(
            group_ht.coverage, group_ht.coverage_model
        ),
    )
    group_ht = group_ht.annotate(exp_agg=group_ht.mu_adj * group_ht.coverage_correction)

    logger.info(
        "Annotating expected counts per allele by distributing expected counts equally"
        " among all alleles in each grouping..."
    )
    group_ht = group_ht.annotate(
        expected=group_ht.exp_agg / group_ht.n_grouping
    ).select("expected", "coverage_correction")
    context_ht = context_ht.annotate(
        **group_ht.index(*[context_ht[g] for g in groupings])
    )
    return context_ht


def create_filtered_context_ht(
    csq: Set[str] = KEEP_CODING_CSQ,
    n_partitions: int = 30000,
    overwrite: bool = False,
    build_models_from_scratch: bool = False,
) -> None:
    """
    Create allele-level VEP context Table with constraint annotations including expected variant counts.

    This Table is used to create the constraint prep Table.

    Table contains only missense, nonsense, read-through, and synonymous alleles in all canonical
    protein-coding transcripts as annotated by VEP. Table is filtered to alleles not found
    or rare in gnomAD exomes at covered sites.

    :param csq: Variant consequences to filter Table to. Default is `KEEP_CODING_CSQ`.
    :param n_partitions: Number of desired partitions for the Table. Default is 30000.
    :param overwrite: Whether to overwrite temporary data. Default is False.
    :param build_models_from_scratch: Whether to build plateau and coverage models from scratch. Default is False.
    :return: None; writes Table to path.
    """
    logger.info(
        "Preprocessing VEP context HT to filter to missense, nonsense, and"
        " synonymous variants in all canonical transcripts and add constraint"
        " annotations..."
    )
    # NOTE: Constraint outlier transcripts are not removed
    ht = process_context_ht(filter_csq=csq)

    logger.info(
        "Filtering context HT to all covered sites not found or rare in gnomAD"
        " exomes and checkpointing..."
    )
    ht = filter_context_using_gnomad(ht, "exomes")
    # Reducing number of partitions as the VEP context table has 62k
    ht = ht.naive_coalesce(n_partitions)
    ht = ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/processed_context.ht",
        _read_if_exists=not overwrite,
        overwrite=overwrite,
    )

    # TODO: Import models built in gnomad-constraint rather than rebuilding here
    if build_models_from_scratch:
        logger.info("Building plateau and coverage models...")
        coverage_ht = prop_obs_coverage.ht()
        coverage_x_ht = hl.read_table(prop_obs_coverage.path.replace(".ht", "_x.ht"))
        coverage_y_ht = hl.read_table(prop_obs_coverage.path.replace(".ht", "_y.ht"))
        (
            coverage_model,
            plateau_models,
            plateau_x_models,
            plateau_y_models,
        ) = generate_models(
            coverage_ht,
            coverage_x_ht,
            coverage_y_ht,
        )
        # Write out models to HailExpression to save
        hl.experimental.write_expression(
            hl.struct(
                coverage=coverage_model,
                plateau_models=plateau_models,
                plateau_X=plateau_x_models,
                plateau_Y=plateau_y_models,
            ),
            coverage_plateau_models_path,
            overwrite=overwrite,
        )
    check_file_exists_raise_error(
        coverage_plateau_models_path,
        error_if_not_exists=True,
        error_if_not_exists_msg=(
            "Coverage and plateau models HailExpression does not exist!"
            " Please double check and/or rerun with"
            " `build_models_from_scratch` = True"
        ),
    )
    models = hl.experimental.read_expression(coverage_plateau_models_path)

    # Also annotate as HT globals
    ht = ht.annotate_globals(
        plateau_models=models.plateau_models,
        plateau_x_models=models.plateau_X,
        plateau_y_models=models.plateau_Y,
        coverage_model=models.coverage,
    )

    logger.info("Calculating expected values per allele and checkpointing...")
    ht = (
        calculate_exp_from_mu(
            filter_to_region_type(ht, "autosomes"),
            locus_type="autosomes",
        )
        .union(calculate_exp_from_mu(filter_to_region_type(ht, "chrX"), locus_type="X"))
        .union(calculate_exp_from_mu(filter_to_region_type(ht, "chrY"), locus_type="Y"))
    )
    ht = ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/context_exp.ht",
        _read_if_exists=not overwrite,
        overwrite=overwrite,
    )

    logger.info("Removing alleles with negative expected values...")
    # Negative expected values happen for a handful of alleles on chrY due to
    # a negative coefficient in the model on this chr
    # We decided to just remove these sites
    ht = ht.filter(ht.expected >= 0)

    logger.info(
        "Annotating context HT with number of observed variants and writing out..."
    )
    ht = add_obs_annotation(ht)
    ht.write(filtered_context.path, overwrite=True)


def create_constraint_prep_ht(
    filter_csq: Set[str] = {MISSENSE}, n_partitions: int = 15000, overwrite: bool = True
) -> None:
    """
    Create locus-level constraint prep Table from filtered context Table for all canonical protein-coding transcripts.

    This Table is used in the first step of regional constraint breakpoint search.

    Table can be filtered to variants with specific consequences before aggregation by locus.

    :param filter_csq: Whether to filter Table to specific consequences. Default is True.
    :param n_partitions: Number of desired partitions for the Table. Default is 15000.
    :param csq: Desired consequences. Default is {`MISSENSE`}. Must be specified if filter is True.
    :param overwrite: Whether to overwrite Table. Default is True.
    :return: None; writes Table to path.
    """
    ht = filtered_context.ht()
    if filter_csq:
        logger.info("Filtering to %s...", filter_csq)
        ht = ht.filter(hl.literal(filter_csq).contains(ht.annotation))

    logger.info("Aggregating by locus...")
    # Context HT is keyed by locus and allele, which means there is one row for every possible missense variant
    # This means that any locus could be present up to three times (once for each possible missense)
    ht = ht.group_by("locus", "transcript").aggregate(
        observed=hl.agg.sum(ht.observed),
        expected=hl.agg.sum(ht.expected),
    )

    logger.info("Adding section annotation...")
    # Add transcript start and stop positions from browser HT
    transcript_ht = gene_model.ht().select("start", "stop")
    ht = ht.annotate(**transcript_ht[ht.transcript])
    ht = ht.annotate(section=hl.format("%s_%s_%s", ht.transcript, ht.start, ht.stop))
    ht = ht.key_by("locus", "section").drop("start", "stop", "transcript")

    ht = ht.naive_coalesce(n_partitions)
    ht.write(constraint_prep.path, overwrite=overwrite)


def get_obs_exp_expr(
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
) -> hl.expr.Float64Expression:
    """
    Return observed/expected annotation based on inputs.

    Cap observed/expected (OE) value at 1. This is to avoid pulling out regions that are
    enriched for missense variation. Code in this pipeline is looking for missense constraint,
    so regions with an OE >= 1.0 can be grouped together.

    Function can generate observed/expected values across the entire transcript or
    section of a transcript depending on inputs.
    Function can also generate 'forward' (moving from smaller to larger positions") or
    'reverse' (moving from larger to smaller positions) section obs/exp values.

    :param obs_expr: Expression containing number of observed variants.
    :param exp_expr: Expression containing number of expected variants.
    :return: Observed/expected expression.
    """
    return hl.nanmin(obs_expr / exp_expr, 1)


def get_reverse_cumulative_obs_exp_expr(
    section_obs_expr: hl.expr.Int64Expression,
    section_exp_expr: hl.expr.Float64Expression,
    fwd_cumulative_obs_expr: hl.expr.Int64Expression,
    fwd_cumulative_exp_expr: hl.expr.Float64Expression,
) -> hl.expr.StructExpression:
    """
    Return the "reverse" cumulative section observed and expected variant counts in a struct.

    The reverse counts are the counts moving from larger to smaller positions
    (backwards from the end of the transcript back to the beginning of the transcript).
    Reverse cumulative value = total section value - forward cumulative value.

    :param section_obs_expr: Expression containing total observed variant counts for transcripts or transcript subsections.
    :param section_exp_expr: Expression containing total expected variant counts for transcripts or transcript subsections.
    :param fwd_cumulative_obs_expr: Expression containing cumulative observed variant counts per site.
    :param fwd_cumulative_exp_expr: Expression containing cumulative expected variant counts per site.
    :return: Struct with reverse cumulative observed and expected variant counts per site.
    """
    return hl.struct(
        # NOTE: Adding hl.max to exp expression to make sure reverse exp is never negative
        # Without this, ran into errors where reverse exp was -5e-14
        # Picked 1e-09 here as tiny number that is not 0
        # ExAC code also did not allow reverse exp to be zero, as this breaks the likelihood ratio tests
        reverse_cumulative_obs=section_obs_expr - fwd_cumulative_obs_expr,
        reverse_cumulative_exp=hl.nanmax(
            section_exp_expr - fwd_cumulative_exp_expr, 1e-09
        ),
    )


def annotate_fwd_exprs(ht: hl.Table) -> hl.Table:
    """
    Annotate input Table with the forward section cumulative observed, expected, and observed/expected values.

    .. note::
        'Forward' refers to moving through the transcript from smaller to larger chromosomal positions.

    Expects input HT to contain the following fields:
        - section
        - observed
        - expected

    Adds the following fields:
        - fwd_cumulative_obs
        - fwd_cumulative_exp
        - fwd_oe

    :param ht: Input Table.
    :return: Table with forward values (cumulative obs, exp, and forward o/e) annotated.
    """
    logger.info("Getting forward cumulative observed and expected variant counts...")
    ht = ht.annotate(
        fwd_cumulative_obs=adjust_fwd_cumulative_count_expr(
            cumulative_count_expr=get_fwd_cumulative_count_expr(
                group_expr=ht.section,
                count_expr=ht.observed,
            ),
            count_expr=ht.observed,
            group_expr=ht.section,
        ),
        fwd_cumulative_exp=adjust_fwd_cumulative_count_expr(
            cumulative_count_expr=get_fwd_cumulative_count_expr(
                group_expr=ht.section,
                count_expr=ht.expected,
            ),
            count_expr=ht.expected,
            group_expr=ht.section,
        ),
    )

    logger.info("Getting forward observed/expected counts and returning...")
    return ht.annotate(
        fwd_oe=get_obs_exp_expr(
            obs_expr=ht.fwd_cumulative_obs, exp_expr=ht.fwd_cumulative_exp
        )
    )


def annotate_reverse_exprs(ht: hl.Table) -> hl.Table:
    """
    Annotate input Table with the reverse section cumulative observed, expected, and observed/expected values.

    .. note::
        'Reverse' refers to moving through the transcript from larger to smaller chromosomal positions.

    Expects input HT to contain the following fields:
        - section_obs
        - section_exp
        - fwd_cumulative_obs
        - fwd_cumulative_exp

    Adds the following fields:
        - reverse_cumulative_obs
        - reverse_cumulative_exp
        - reverse_oe

    :param ht: Input Table.
    :return: Table with reverse values annotated.
    """
    logger.info("Getting reverse cumulative observed and expected variant counts...")
    # Reverse cumulative value = total value - forward cumulative value
    ht = ht.annotate(
        **get_reverse_cumulative_obs_exp_expr(
            section_obs_expr=ht.section_obs,
            section_exp_expr=ht.section_exp,
            fwd_cumulative_obs_expr=ht.fwd_cumulative_obs,
            fwd_cumulative_exp_expr=ht.fwd_cumulative_exp,
        )
    )
    logger.info("Getting reverse observed/expected counts and returning...")
    # Set reverse o/e to missing if reverse expected value is 0 (to avoid NaNs)
    return ht.annotate(
        reverse_oe=get_obs_exp_expr(
            ht.reverse_cumulative_obs, ht.reverse_cumulative_exp
        )
    )


def get_dpois_expr(
    oe_expr: hl.expr.Float64Expression,
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
) -> hl.expr.StructExpression:
    """
    Calculate probability densities (natural log) of the observed values under a Poisson model.

    Poisson rate for each region is defined as expected * given observed/expected values.

    :param oe_expr: Expression of observed/expected value.
    :param obs_expr: Expression containing observed variants count.
    :param exp_expr: Expression containing expected variants count.
    :return: Natural log of the probability density under Poisson model.
    """
    # log_p = True returns the natural logarithm of the probability density
    return hl.dpois(obs_expr, exp_expr * oe_expr, log_p=True)


def annotate_max_chisq_per_section(
    ht: hl.Table,
    freeze: int,
) -> hl.Table:
    """
    Get maximum chi square value per transcript or transcript subsection.

    Expects input HT to contain the following fields:
        - section
        - chisq

    Adds the following field:
        - section_max_chisq

    :param ht: Input Table.
    :param freeze: RMC data freeze number.
    :return: Table annotated with maximum chi square value per transcript or transcript subsection.
    """
    group_ht = ht.group_by("section").aggregate(
        # hl.agg.max ignores NaNs
        section_max_chisq=hl.agg.max(ht.chisq)
    )
    group_ht = group_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/freeze{freeze}_group_max_chisq.ht", overwrite=True
    )
    return ht.annotate(section_max_chisq=group_ht[ht.section].section_max_chisq)


def search_for_break(
    ht: hl.Table,
    search_num: int,
    freeze: int,
    chisq_threshold: float = scipy.stats.chi2.ppf(1 - P_VALUE, 1),
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
        - fwd_cumulative_exp
        - fwd_cumulative_obs
        - fwd_oe
        - reverse_cumulative_obs
        - reverse_cumulative_exp
        - reverse_oe
        - section_obs
        - section_exp
        - section_oe

    Also expects:
        - Input HT was created using a VEP context HT.
        - Multiallelic variants in input HT have been split.

    :param ht: Input Table.
    :param search_num: Search iteration number (e.g., second round of searching for single break would be 2).
    :param freeze: RMC data freeze number.
    :param chisq_threshold: Chi-square significance threshold.
        Default is `scipy.stats.chi2.ppf(1 - P_VALUE, 1)`.
        Default value used in ExAC was 10.8, which corresponds to a p-value of 0.001
        with 1 degree of freedom.
        (https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm)
    :param min_num_exp_mis: Minimum number of expected missense per transcript/transcript section.
        Sections that have fewer than this number of expected missense variants will not
        be computed (chi square will be annotated as a missing value).
        Default is `MIN_EXP_MIS`.
    :param save_chisq_ht: Whether to save HT with chi square values annotated for every locus.
        This saves a lot of extra data and should only occur once.
        Default is False.
    :return: Table annotated with whether position is a breakpoint (`is_break`).
    """
    logger.info(
        "Creating section null (no regional variability in missense depletion)"
        " and alt (evidence of domains of missense constraint) probability densities..."
    )
    logger.info(
        (
            "Skipping breakpoints that create at least one section that has < %i"
            " expected missense variants..."
        ),
        min_num_exp_mis,
    )
    ht = ht.annotate(
        total_null=hl.or_missing(
            (ht.fwd_cumulative_exp >= min_num_exp_mis)
            & (ht.reverse_cumulative_exp >= min_num_exp_mis),
            # Add forwards section null (going through positions from smaller to larger)
            # section_null = stats.dpois(section_obs, section_exp*overall_obs_exp)[0]
            get_dpois_expr(
                oe_expr=ht.section_oe,
                obs_expr=ht.fwd_cumulative_obs,
                exp_expr=ht.fwd_cumulative_exp,
            )
            # Add reverse section null (going through positions from larger to smaller)
            + get_dpois_expr(
                oe_expr=ht.section_oe,
                obs_expr=ht.reverse_cumulative_obs,
                exp_expr=ht.reverse_cumulative_exp,
            ),
        )
    )
    ht = ht.annotate(
        total_alt=hl.or_missing(
            hl.is_defined(ht.total_null),
            # Add forward section alt
            # section_alt = stats.dpois(section_obs, section_exp*section_obs_exp)[0]
            get_dpois_expr(
                oe_expr=ht.fwd_oe,
                obs_expr=ht.fwd_cumulative_obs,
                exp_expr=ht.fwd_cumulative_exp,
            )
            # Add reverse section alt
            + get_dpois_expr(
                oe_expr=ht.reverse_oe,
                obs_expr=ht.reverse_cumulative_obs,
                exp_expr=ht.reverse_cumulative_exp,
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

    ht = annotate_max_chisq_per_section(ht, freeze)
    return ht.annotate(
        is_break=((ht.chisq == ht.section_max_chisq) & (ht.chisq >= chisq_threshold))
    )


def annotate_subsection_exprs(ht: hl.Table) -> hl.Table:
    """
    Annotate total observed, expected, and observed/expected (OE) counts for each section of a transcript.

    Expects input HT to contain the following fields:
        - section
        - observed
        - expected

    Adds the following fields:
        - section_obs
        - section_exp
        - section_oe

    :param ht: Input Table.
    :return: Table annotated with section observed, expected, and OE counts.
    """
    logger.info(
        "Getting total observed and expected counts for each transcript or transcript"
        " subsection..."
    )
    group_ht = ht.group_by("section").aggregate(
        section_obs=hl.agg.sum(ht.observed),
        section_exp=hl.agg.sum(ht.expected),
    )
    ht = ht.annotate(**group_ht[ht.section])
    logger.info(
        "Getting observed/expected value for each transcript or transcript"
        " subsection..."
    )
    return ht.annotate(
        section_oe=get_obs_exp_expr(obs_expr=ht.section_obs, exp_expr=ht.section_exp)
    )


def process_sections(
    ht: hl.Table,
    search_num: int,
    freeze: int,
    chisq_threshold: float = scipy.stats.chi2.ppf(1 - P_VALUE, 1),
    save_chisq_ht: bool = False,
):
    """
    Search for breaks within given transcripts or transcript subsections using per-site observed and expected variant counts.

    Expects that input Table has the following annotations:
        - section
        - observed
        - expected

    :param ht: Input Table.
    :param search_num: Search iteration number (e.g., second round of searching for single break would be 2).
    :param freeze: RMC data freeze number.
    :param chisq_threshold: Chi-square significance threshold. See docstring for `search_for_break` for details.
    :param save_chisq_ht: Whether to save HT with chi square values annotated for every locus.
        This saves a lot of extra data and should only occur during the initial search round.
        Default is False.
    :return: Table annotated with whether position is a breakpoint.
    """
    logger.info("Annotating section total observed, expected, and obs/exp...")
    ht = annotate_subsection_exprs(ht)

    logger.info("Annotating forward cumulative observed, expected, and obs/exp...")
    ht = annotate_fwd_exprs(ht)

    logger.info("Annotating reverse cumulative observed, expected, and obs/exp...")
    ht = annotate_reverse_exprs(ht)
    tmp_obs_exp_annot_path = f"{TEMP_PATH_WITH_FAST_DEL}/rmc/freeze{freeze}_single_search_prep_round{search_num}_chisq{chisq_threshold}.ht"
    ht = ht.checkpoint(tmp_obs_exp_annot_path, overwrite=True)

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
    keep_annotations: Set[str] = FINAL_ANNOTATIONS,
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
            (
                "Simul breaks results HT (round %i) did not exist. Please double check"
                " that this was expected!"
            ),
            search_num,
        )
    return ht


def calculate_oe_neq_1_chisq(
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
) -> hl.expr.Float64Expression:
    """
    Check for significance that observed/expected values for regions are different from 1.

    Formula is: (obs - exp)^2 / exp.

    :param obs_expr: Observed variant counts.
    :param exp_expr: Expected variant counts.
    :return: Chi-squared value.
    """
    return ((obs_expr - exp_expr) ** 2) / exp_expr


def merge_rmc_hts(
    round_nums: List[int],
    freeze: int,
    overwrite_temp: bool = False,
) -> hl.Table:
    """
    Get table of final RMC sections after all break searches are complete.

    :param round_nums: List of round numbers to merge results across.
    :param freeze: RMC data freeze number.
    :param overwrite_temp: Whether to overwrite intermediate temporary data if it already exists.
        If False, will read existing intermediate temporary data rather than overwriting.
        Default is False.
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
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )
    # Calculate chi-square value for each section having O/E different from 1
    rmc_ht = rmc_ht.annotate(
        section_chisq=calculate_oe_neq_1_chisq(rmc_ht.section_obs, rmc_ht.section_exp)
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
        - Assumes `constraint_prep` Table is missense-specific
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
    rmc_prep_ht = constraint_prep.ht().select_globals()
    # Add transcript annotation from section field as this is required for joins to other tables
    rmc_prep_ht = rmc_prep_ht.annotate(transcript=rmc_prep_ht.section.split("_")[0])
    group_rmc_prep_ht = rmc_prep_ht.group_by("transcript").aggregate(
        obs=hl.agg.sum(rmc_prep_ht.observed),
        exp=hl.agg.sum(rmc_prep_ht.expected),
    )
    # Recalculating transcript level OE ratio because previous OE ratio (`overall_oe`)
    # is capped at 1 for regional missense constraint calculation purposes
    group_rmc_prep_ht = group_rmc_prep_ht.annotate(
        transcript_oe=group_rmc_prep_ht.obs / group_rmc_prep_ht.exp
    )

    # Read in LoF constraint HT to get OE ratio for five transcripts missing in v2 RMC results
    # # 'ENST00000304270', 'ENST00000344415', 'ENST00000373521', 'ENST00000381708', 'ENST00000596936'
    # All 5 of these transcripts have extremely low coverage in gnomAD
    # Will keep for consistency with v2 LoF results but they look terrible
    # NOTE: LoF HT is keyed by gene and transcript, but `_key_by_assert_sorted` doesn't work here for v2 version
    # Throws this error: hail.utils.java.FatalError: IllegalArgumentException
    lof_ht = constraint_ht.ht().select("oe_mis").key_by("transcript")

    ht = ht.annotate(
        gnomad_transcript_oe=lof_ht[ht.transcript].oe_mis,
        rmc_transcript_oe=group_rmc_prep_ht[ht.transcript].transcript_oe,
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
        "Reading in SNPs-only, VEP-annotated context HT and filtering to missense"
        " variants in canonical transcripts..."
    )
    # Using vep_context.path to read in table with fewer partitions
    # VEP context resource has 62164 partitions
    ht = (
        hl.read_table(vep_context.path, _n_partitions=n_partitions)
        .select_globals()
        .select("vep", "was_split")
    )
    ht = process_vep(ht, filter_csq={missense_str}, filter_outlier_transcripts=True)
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
    if not file_exists(clinvar_plp_mis_haplo.path):
        import_clinvar(overwrite=True)
    if not file_exists(ndd_de_novo.path):
        import_de_novo_variants(overwrite=True)

    clinvar_ht = clinvar_plp_mis_haplo.ht()
    dn_ht = ndd_de_novo.ht()
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
