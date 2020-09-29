import logging
from typing import Dict, Tuple

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad_lof.constraint_utils.constraint_basics import (
    add_most_severe_csq_to_tc_within_ht,
    annotate_constraint_groupings,
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
from rmc.resources.resource_utils import BUILDS


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
    build: str, trimers: bool = True, overwrite: bool = True, n_partitions: int = 1000
) -> None:
    """
    Imports reference fasta (SNPs only, annotated with VEP) as a hail Table.

    Filters to missense variants in canonical protein coding transcripts. 
    Also annotates with probability of mutation for each variant.

    .. note::
        `trimers` needs to be True for gnomAD v2.

    :param str build: Reference genome build; must be one of BUILDS.
    :param bool trimers: Whether to filter to trimers or heptamers. Default is True.
    :param bool overwrite: Whether to overwrite output. Default is True.
    :param int n_partitions: Number of desired partitions for output. Default is 1000.
    :return: None
    :rtype: None
    """
    if build not in BUILDS:
        raise DataException(f"Build must be one of {BUILDS}.")

    logger.info("Reading in SNPs-only, VEP-annotated context ht...")

    if build == "GRCh37":
        ht = grch37.full_context.ht()
        output_path = grch37.processed_context.path
    else:
        ht = grch38.full_context.ht()
        output_path = grch38.processed_context.path

    # prepare_ht annotates HT with: ref, alt, methylation_level, exome_coverage, cpg, transition, variant_type
    ht = prepare_ht(ht, trimers)

    logger.info(
        "Filtering to missense variants in canonical protein coding transcripts..."
    )
    ht = filter_to_missense(ht)

    logger.info("Annotating with mutation rate...")
    # Mutation rate HT is keyed by context, ref, alt, methylation level
    mu_ht = mutation_rate.ht()
    ht, grouping = annotate_constraint_groupings(ht)
    ht = ht.filter(hl.is_defined(ht.exome_coverage))
    ht = ht.select(
        "context", "ref", "alt", "methylation_level", "exome_coverage", *grouping
    )
    ht = ht.annotate(
        mu_snp=mu_ht[ht.context, ht.ref, ht.alt, ht.methylation_level].mu_snp
    )

    logger.info("Writing out context HT...")
    ht = ht.naive_coalesce(n_partitions)
    ht.write(output_path, overwrite=overwrite)


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
    ht = ht.annotate(cds_len=(ht.interval.end - ht.interval.start) + 1)
    return ht.aggregate(hl.agg.sum(ht.cds_len))


def keep_criteria(ht: hl.Table, exac: bool) -> hl.expr.BooleanExpression:
    """
    Returns Boolean expression to filter variants in input Table.

    :param hl.Table ht: Input Table.
    :param bool exac: Whether input Table is ExAC data.
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
        keep_criteria = (ht.ac > 0) & (ht.af < 0.001) & (ht.pass_filters)

    return keep_criteria


def filter_to_missense(ht: hl.Table, n_partitions: int = 5000) -> hl.Table:
    """
    Filters input Table to missense variants in canonical transcripts only.

    :param Table ht: Input Table to be filtered.
    :param int n_partitions: Number of desired partitions for output.
    :return: Table filtered to only missense variants.
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

    logger.info("Filtering to missense variants...")
    ht = ht.filter(
        ht.transcript_consequences.most_severe_consequence == "missense_variant"
    )
    return ht.naive_coalesce(n_partitions)


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
        .when(coverage == 0, 0.0)
        .when(
            (coverage >= 1) & (coverage < high_cov_cutoff),
            coverage_model[1] * hl.log(coverage) + coverage_model[0],
        )
        .default(1.0)
    )


def get_plateau_model(
    locus_expr: hl.expr.LocusExpression,
    cpg_expr: hl.expr.BooleanExpression,
    globals_expr: hl.expr.StructExpression,
) -> hl.expr.Float64Expression:
    """
    Gets model to determine adjustment to mutation rate based on locus type and CpG status.

    .. note::
        This function expects that the context Table has each plateau model (autosome, X, Y) added as global annotations.

    :param hl.expr.LocusExpression locus_expr: Locus expression.
    :param hl.expr.BooleanExpression: Expression describing whether site is a CpG site.
    :param hl.expr.StructExpression globals_expr: Expression containing global annotations of context HT. Must contain plateau models as annotations.
    :return: Plateau model for locus type.
    :rtype: hl.expr.Float64Expression
    """
    return (
        hl.case()
        .when(locus_expr.in_x_nonpar(), globals_expr.plateau_x_models.total[cpg_expr])
        .when(locus_expr.in_y_nonpar(), globals_expr.plateau_y_models.total[cpg_expr])
        .default(globals_expr.plateau_model.total[cpg_expr])
    )
