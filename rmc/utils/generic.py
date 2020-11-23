import logging
from typing import Dict, Tuple

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad_lof.constraint_utils.constraint_basics import (
    add_most_severe_csq_to_tc_within_ht,
    annotate_constraint_groupings,
    annotate_with_mu,
    build_models,
    prepare_ht,
)
from gnomad_lof.constraint_utils.generic import fast_filter_vep
from rmc.resources.basics import (
    ACID_NAMES_PATH,
    CODON_TABLE_PATH,
    DIVERGENCE_SCORES_TSV_PATH,
    mutation_rate,
    MUTATION_RATE_TABLE_PATH,
)
import rmc.resources.grch37.reference_data as grch37
import rmc.resources.grch38.reference_data as grch38
from rmc.resources.resource_utils import BUILDS, MISSENSE


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_missense_constraint_generic")
logger.setLevel(logging.INFO)


## Resources from Kaitlin
def get_codon_lookup() -> Dict:
    """
    Reads in codon lookup table and returns as dictionary (key: codon, value: amino acid)

    .. note:: 
        This is only necessary for testing on ExAC and should be replaced with VEP annotations.

    :return: Dictionary of codon translation
    :rtype: dict
    """
    codon_lookup = {}
    with hl.hadoop_open(CODON_TABLE_PATH) as c:
        c.readline()
        for line in c:
            line = line.strip().split(" ")
            codon_lookup[line[0]] = line[1]
    return codon_lookup


def get_acid_names() -> Dict:
    """
    Reads in amino acid table and stores as dict (key: 3 letter name, value: (long name, one letter name)

    :return: Dictionary of amino acid names
    :rtype: dict
    """
    acid_map = {}
    with hl.hadoop_open(ACID_NAMES_PATH) as a:
        a.readline()
        for line in a:
            line = line.strip().split("\t")
            acid_map[line[1]] = (line[0], line[2])
    return acid_map


def get_mutation_rate() -> Dict:
    """
    Reads in mutation rate table and stores as dict

    :return: Dictionary of mutation rate information (key: context, value: (alt, mu_snp))
    :rtype: dict
    """
    mu = {}
    # from    n_kmer  p_any_snp_given_kmer    mu_kmer to      count_snp       p_snp_given_kmer        mu_snp
    with hl.hadoop_open(MUTATION_RATE_TABLE_PATH) as m:
        for line in m:
            context, _, _, _, new_kmer, _, _, mu_snp = line.strip().split("\t")
            mu[context] = (new_kmer[1], mu_snp)
    return mu


def get_divergence_scores() -> Dict:
    """
    Reads in divergence score file and stores as dict (key: transcript, value: score)

    :return: Divergence score dict
    :rtype: dict
    """
    div_scores = {}
    with hl.hadoop_open(DIVERGENCE_SCORES_TSV_PATH) as d:
        d.readline()
        for line in d:
            transcript, score = line.strip().split("\t")
            try:
                div_scores[transcript.split(".")[0]] = float(score)
            except:
                continue
    return div_scores


## Functions to process reference genome related resources
def process_context_ht(
    build: str, trimers: bool = True, missense_str: str = MISSENSE,
) -> hl.Table:
    """
    Prepares context HT (SNPs only, annotated with VEP) for regional missense constraint calculations.

    Filters to missense variants in canonical protein coding transcripts. 
    Also annotates with probability of mutation for each variant, CpG status, and methylation level.

    .. note::
        `trimers` needs to be True for gnomAD v2.

    :param str build: Reference genome build; must be one of BUILDS.
    :param bool trimers: Whether to filter to trimers or heptamers. Default is True.
    :return: Context HT filtered to missense variants in canonical transcripts and annotated with mutation rate, CpG status, and methylation level.
    :rtype: hl.Table
    """
    if build not in BUILDS:
        raise DataException(f"Build must be one of {BUILDS}.")

    logger.info("Reading in SNPs-only, VEP-annotated context ht...")

    if build == "GRCh37":
        ht = grch37.full_context.ht()
    else:
        ht = grch38.full_context.ht()

    # `prepare_ht` annotates HT with: ref, alt, methylation_level, exome_coverage, cpg, transition, variant_type
    ht = prepare_ht(ht, trimers)

    logger.info(
        f"Filtering to canonical transcripts, annotating with most severe consequence, and filtering to {missense_str}..."
    )
    ht = process_vep(ht, filter_csq=True, csq=missense_str)

    logger.info("Annotating with mutation rate...")
    # Mutation rate HT is keyed by context, ref, alt, methylation level
    mu_ht = mutation_rate.ht().select("mu_snp")
    ht, grouping = annotate_constraint_groupings(ht)
    ht = ht.select(
        "context",
        "ref",
        "alt",
        "methylation_level",
        "exome_coverage",
        "cpg",
        "transition",
        "variant_type",
        *grouping,
    )
    return annotate_with_mu(ht, mu_ht)


## Functions for obs/exp related resources
def get_exome_bases(build: str) -> int:
    """
    Imports gencode gtf into a hail Table, filters to coding regions, and sums all positions to get the number of bases in the exome.

    :param str build: Reference genome build; must be one of BUILDS.
    :return: Number of bases in the exome.
    :rtype: int
    """
    if build not in BUILDS:
        raise DataException(f"Build must be one of {BUILDS}.")

    logger.info("Reading in gencode gtf")
    if build == "GRCh37":
        ht = grch37.gencode.ht()

    else:
        ht = grch38.gencode.ht()

    logger.info("Filtering gencode gtf to feature type == CDS...")
    ht = ht.filter((ht.feature == "CDS"))

    logger.info("Summing total bases in exome...")
    # GENCODE is 1-based
    ht = ht.annotate(
        cds_len=(ht.interval.end.position - ht.interval.start.position) + 1
    )
    return ht.aggregate(hl.agg.sum(ht.cds_len))


def keep_criteria(ht: hl.Table, exac: bool = False) -> hl.expr.BooleanExpression:
    """
    Returns Boolean expression to filter variants in input Table.

    :param hl.Table ht: Input Table.
    :param bool exac: Whether input Table is ExAC data. Default is False.
    :return: Keep criteria Boolean expression.
    :rtype: hl.expr.BooleanExpression
    """
    # ExAC keep criteria: adjusted AC <= 123 and VQSLOD >= -2.632
    # Also remove variants with median depth < 1
    if exac:
        keep_criteria = (
            (ht.ac <= 123) & (ht.ac > 0) & (ht.vqslod >= -2.632) & (ht.coverage > 1)
        )
    else:
        keep_criteria = (
            (ht.ac > 0) & (ht.af < 0.001) & (ht.pass_filters) & (ht.exome_coverage > 0)
        )

    return keep_criteria


def process_vep(ht: hl.Table, filter_csq: bool = False, csq: str = None) -> hl.Table:
    """
    Filters input Table to canonical transcripts only.

    Option to filter Table to specific variant consequence (csq).

    :param Table ht: Input Table.
    :param bool filter: Whether to filter Table to a specific consequence. Default is False.
    :param str csq: Desired consequence. Default is None. Must be specified if filter is True.
    :return: Table filtered to canonical transcripts with option to filter to specific variant consequence.
    :rtype: hl.Table
    """
    if "was_split" not in ht.row:
        logger.info("Splitting multiallelic variants and filtering to SNPs...")
        ht = hl.split_multi(ht)
        ht = ht.filter(hl.is_snp(ht.alleles[0], ht.alleles[1]))

    logger.info("Filtering to canonical transcripts...")
    ht = fast_filter_vep(
        ht, vep_root="vep", syn=False, canonical=True, filter_empty=True
    )

    logger.info("Annotating HT with most severe consequence...")
    ht = add_most_severe_csq_to_tc_within_ht(ht)
    ht = ht.transmute(transcript_consequences=ht.vep.transcript_consequences)
    ht = ht.explode(ht.transcript_consequences)

    if filter_csq:
        logger.info(f"Filtering to {csq}...")
        ht = ht.filter(ht.transcript_consequences.most_severe_consequence == csq)
    return ht


def filter_to_region_type(ht: hl.Table, region: str) -> hl.Table:
    """
    Filters input Table to autosomes + chrX PAR, chrX non-PAR, or chrY non-PAR. 

    :param hl.Table ht: Input Table to be filtered.
    :param str region: Desired region type. One of 'autosomes', 'chrX', or 'chrY'.
    :return: Table filtered to autosomes/PAR, chrX, or chrY.
    :rtype: hl.Table
    """
    if region == "chrX":
        ht = ht.filter(ht.locus.in_x_nonpar())
    elif region == "chrY":
        ht = ht.filter(ht.locus.in_y_nonpar())
    else:
        ht = ht.filter(ht.locus.in_autosome() | ht.locus.in_x_par())
    return ht


def generate_models(
    coverage_ht: hl.Table,
    coverage_x_ht: hl.Table,
    coverage_y_ht: hl.Table,
    trimers: bool,
    weighted: bool = True,
) -> Tuple[
    Tuple[float, float],
    hl.expr.DictExpression,
    hl.expr.DictExpression,
    hl.expr.DictExpression,
]:
    """
    Calls `build_models` from gnomAD LoF repo to generate models used to adjust expected variants count.

    :param hl.Table coverage_ht: Table with proportion of variants observed by coverage (autosomes/PAR only).
    :param hl.Table coverage_x_ht: Table with proportion of variants observed by coverage (chrX only).
    :param hl.Table coverage_y_ht: Table with proportion of variants observed by coverage (chrY only).
    :param bool trimers: Whether to use trimers instead of heptamers.
    :param bool weighted: Whether to use weighted least squares when building models. Default is True.
    :return: Coverage model, plateau models for autosomes, plateau models for chrX, plateau models for chrY. 
    :rtype: Tuple[Tuple[float, float], hl.expr.DictExpression, hl.expr.DictExpression, hl.expr.DictExpression]
    """
    logger.info("Building autosomes/PAR plateau model and coverage model...")
    coverage_model, plateau_models = build_models(
        coverage_ht, trimers=trimers, weighted=weighted
    )

    logger.info("Building plateau models for chrX and chrY...")
    # TODO: make half_cutoff (for coverage cutoff) True for X/Y?
    # This would also mean saving new coverage model for allosomes
    _, plateau_x_models = build_models(
        coverage_x_ht, trimers=trimers, weighted=weighted
    )
    _, plateau_y_models = build_models(
        coverage_y_ht, trimers=trimers, weighted=weighted
    )
    return (coverage_model, plateau_models, plateau_x_models, plateau_y_models)


def get_coverage_correction_expr(
    coverage: hl.expr.Float64Expression,
    coverage_model: Tuple[float, float],
    high_cov_cutoff: int = 40,
) -> hl.expr.Float64Expression:
    """
    Gets coverage correction for expected variants count.

    .. note:: 
        Default high coverage cutoff taken from gnomAD LoF repo.

    :param hl.expr.Float64Expression ht: Input coverage expression. Should be median coverage at position.
    :param Tuple[float, float] coverage_model: Model to determine coverage correction factor necessary
         for calculating expected variants at low coverage sites.
    :param int high_cov_cutoff: Cutoff for high coverage. Default is 40. 
    :return: Coverage correction expression.
    :rtype: hl.expr.Float64Expression
    """
    return (
        hl.case()
        .when(coverage == 0, 0)
        .when(coverage >= high_cov_cutoff, 1)
        .default(coverage_model[1] * hl.log10(coverage) + coverage_model[0])
    )


def get_plateau_model(
    locus_expr: hl.expr.LocusExpression,
    cpg_expr: hl.expr.BooleanExpression,
    globals_expr: hl.expr.StructExpression,
    include_cpg: bool = False,
) -> hl.expr.Float64Expression:
    """
    Gets model to determine adjustment to mutation rate based on locus type and CpG status.

    .. note::
        This function expects that the context Table has each plateau model (autosome, X, Y) added as global annotations.

    :param hl.expr.LocusExpression locus_expr: Locus expression.
    :param hl.expr.BooleanExpression: Expression describing whether site is a CpG site. Required if include_cpg is False.
    :param hl.expr.StructExpression globals_expr: Expression containing global annotations of context HT. Must contain plateau models as annotations.
    :param bool include_cpg: Whether to return full plateau model dictionary including CpG keys. Default is False.
    :return: Plateau model for locus type.
    :rtype: hl.expr.Float64Expression
    """
    if include_cpg:
        return (
            hl.case()
            .when(locus_expr.in_x_nonpar(), globals_expr.plateau_x_models.total)
            .when(locus_expr.in_y_nonpar(), globals_expr.plateau_y_models.total)
            .default(globals_expr.plateau_models.total)
        )

    return (
        hl.case()
        .when(locus_expr.in_x_nonpar(), globals_expr.plateau_x_models.total[cpg_expr])
        .when(locus_expr.in_y_nonpar(), globals_expr.plateau_y_models.total[cpg_expr])
        .default(globals_expr.plateau_models.total[cpg_expr])
    )
