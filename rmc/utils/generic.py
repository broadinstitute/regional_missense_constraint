import logging
from typing import Dict, Tuple

import hail as hl

from gnomad.resources.grch37.gnomad import coverage, public_release
from gnomad.resources.grch37.reference_data import vep_context
from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists
from gnomad.utils.filtering import filter_to_clinvar_pathogenic

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
    hi_genes,
    mutation_rate,
    MUTATION_RATE_TABLE_PATH,
    TOTAL_EXOME_BASES,
    TOTAL_GNOMAD_MISSENSE,
)
from rmc.resources.grch37.gnomad import constraint_ht, filtered_exomes
import rmc.resources.grch37.reference_data as grch37
import rmc.resources.grch38.reference_data as grch38
from rmc.resources.resource_utils import BUILDS, CURRENT_VERSION, MISSENSE


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_missense_constraint_generic")
logger.setLevel(logging.INFO)


## Resources from Kaitlin
def get_codon_lookup() -> hl.expr.DictExpression:
    """
    Read in codon lookup table and return as dictionary (key: codon, value: amino acid).

    .. note::
        This is only necessary for testing on ExAC and should be replaced with VEP annotations.

    :return: DictExpression of codon translation.
    """
    codon_lookup = {}
    with hl.hadoop_open(CODON_TABLE_PATH) as c:
        c.readline()
        for line in c:
            line = line.strip().split()
            codon_lookup[line[0]] = line[1]
    return hl.literal(codon_lookup)


def get_aa_map() -> Dict[str, str]:
    """
    Create dictionary mapping amino acid 1 letter name to 3 letter name.

    :return: Dictionary mapping amino acid 1 letter name (key) to 3 letter name (value).
    """
    aa_map = {}
    with hl.hadoop_open(ACID_NAMES_PATH) as a:
        a.readline()
        for line in a:
            line = line.strip().split("\t")
            aa_map[line[2]] = line[1]
    return aa_map


def get_mutation_rate() -> Dict[str, Tuple[str, float]]:
    """
    Read in mutation rate table and store as dict.

    :return: Dictionary of mutation rate information (key: context, value: (alt, mu_snp)).
    :rtype: Dict[str, Tuple[str, float]].
    """
    mu = {}
    # from    n_kmer  p_any_snp_given_kmer    mu_kmer to      count_snp       p_snp_given_kmer        mu_snp
    with hl.hadoop_open(MUTATION_RATE_TABLE_PATH) as m:
        for line in m:
            context, _, _, _, new_kmer, _, _, mu_snp = line.strip().split("\t")
            mu[context] = (new_kmer[1], mu_snp)
    return mu


def get_divergence_scores() -> Dict[str, float]:
    """
    Read in divergence score file and store as dict (key: transcript, value: score).

    :return: Divergence score dict.
    :rtype: Dict[str, float]
    """
    div_scores = {}
    with hl.hadoop_open(DIVERGENCE_SCORES_TSV_PATH) as d:
        d.readline()
        for line in d:
            transcript, score = line.strip().split("\t")
            try:
                div_scores[transcript.split(".")[0]] = float(score)
            except ValueError:
                continue
    return div_scores


## Functions to process reference genome related resources
def process_context_ht(
    build: str,
    trimers: bool = True,
    filter_to_missense: bool = True,
    missense_str: str = MISSENSE,
    add_annotations: bool = True,
) -> hl.Table:
    """
    Prepare context HT (SNPs only, annotated with VEP) for regional missense constraint calculations.

    Filter to missense variants in canonical protein coding transcripts.
    Also annotate with probability of mutation for each variant, CpG status, and methylation level.

    .. note::
        `trimers` needs to be True for gnomAD v2.

    :param str build: Reference genome build; must be one of BUILDS.
    :param bool trimers: Whether to filter to trimers (if set to True) or heptamers. Default is True.
    :param bool filter_to_missense: Whether to filter Table to missense variants only. Default is True.
    :param bool add_annotations: Whether to add ref, alt, methylation_level, exome_coverage, cpg, transition,
        and variant_type annotations. Default is True.
    :return: Context HT filtered to canonical transcripts and optionally filtered to missense variants with
        mutation rate, CpG status, and methylation level annotations.
    :rtype: hl.Table
    """
    if build not in BUILDS:
        raise DataException(f"Build must be one of {BUILDS}.")

    logger.info("Reading in SNPs-only, VEP-annotated context ht...")
    if build == "GRCh37":
        ht = vep_context.ht().select_globals()
    else:
        ht = grch38.full_context.ht()

    logger.info(
        "Filtering to canonical transcripts and annotating with most severe consequence...",
    )
    if filter_to_missense:
        ht = process_vep(ht, filter_csq=True, csq=missense_str)
    else:
        ht = process_vep(ht)

    if add_annotations:
        # `prepare_ht` annotates HT with: ref, alt, methylation_level, exome_coverage, cpg, transition, variant_type
        ht = prepare_ht(ht, trimers)

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
    return ht


def filter_context_using_gnomad(
    context_ht: hl.Table,
    gnomad_data_type: str = "exomes",
    adj_freq_index: int = 0,
    filter_context_using_cov: bool = True,
    cov_threshold: int = 0,
) -> hl.Table:
    """
    Filter VEP context Table to sites that aren't seen in gnomAD or are rare in gnomAD.

    :param hl.Table context_ht: VEP context Table.
    :param str gnomad_data_type: gnomAD data type. Used to retrieve public release and coverage resources.
        Must be one of "exomes" or "genomes" (check is done within `public_release`).
        Default is "exomes".
    :param adj_freq_index: Index of frequency array that contains global population filtered calculated on
        high quality (adj) genotypes. Default is 0.
    :param bool filter_context_using_cov: Whether to also filter sites in context Table using gnomAD coverage.
        Default is True.
    :param int cov_threshold: Coverage threshold used to filter context Table if `filter_context_using_cov` is True.
        Default is 0.
    :return: Filtered VEP context Table.
    """
    gnomad = public_release(gnomad_data_type).ht().select_globals()
    gnomad_cov = coverage(gnomad_data_type).ht()
    gnomad = gnomad.select(
        "filters",
        ac=gnomad.freq[adj_freq_index].AC,
        af=gnomad.freq[adj_freq_index].AF,
        gnomad_coverage=gnomad_cov[gnomad.locus].median,
    )

    # Filter to sites not seen in gnomAD or to rare sites in gnomAD
    gnomad_join = gnomad[context_ht.key]
    context_ht = context_ht.filter(
        hl.is_missing(gnomad_join)
        | keep_criteria(
            gnomad_join.ac,
            gnomad_join.af,
            gnomad_join.filters,
            gnomad_join.gnomad_coverage,
            cov_threshold=cov_threshold,
        )
    )

    # Optionally also filter context HT using gnomAD coverage
    if filter_context_using_cov:
        context_ht = context_ht.annotate(
            gnomad_coverage=gnomad_cov[context_ht.locus].median
        )
        context_ht = context_ht.filter(context_ht.gnomad_coverage > cov_threshold)
        # Drop coverage annotation here for consistency
        # (This annotation is only added if `filter_context_using_cov` is True)
        context_ht = context_ht.drop("gnomad_coverage")

    return context_ht


## Functions for obs/exp related resources
def get_exome_bases(build: str) -> int:
    """
    Get number of bases in the exome.

    Read in context HT (containing all coding bases in the genome), remove outlier transcripts, and filter to median coverage >= 5.

    :param str build: Reference genome build; must be one of BUILDS.
    :return: Number of bases in the exome.
    :rtype: int
    """
    if build not in BUILDS:
        raise DataException(f"Build must be one of {BUILDS}.")

    logger.info("Reading in SNPs-only, VEP-annotated context ht...")
    if build == "GRCh37":
        ht = vep_context.ht()
    else:
        ht = grch38.full_context.ht()

    logger.info("Filtering to canonical transcripts...")
    ht = fast_filter_vep(
        ht, vep_root="vep", syn=False, canonical=True, filter_empty=True
    )

    logger.info("Removing outlier transcripts...")
    outlier_transcripts = get_constraint_transcripts(outlier=True)
    ht = ht.transmute(transcript_consequences=ht.vep.transcript_consequences)
    ht = ht.explode(ht.transcript_consequences)
    ht = ht.filter(
        ~outlier_transcripts.contains(ht.transcript_consequences.transcript_id)
    )

    logger.info(
        "Collecting context HT by key (to make sure each locus only gets counted once)..."
    )
    ht = ht.key_by("locus").collect_by_key()

    logger.info("Removing positions with median coverage < 5...")
    # Taking just the first value since the coverage should be the same for all entries at a locus
    ht = ht.filter(ht.values[0].coverage.exomes.median >= 5)
    ht = ht.select()
    return ht.count()


def get_avg_bases_between_mis(
    build: str,
    get_total_exome_bases: bool = False,
    get_total_gnomad_missense: bool = False,
) -> int:
    """
    Return average number of bases between observed missense variation.

    For example, if the total number of bases is 30, and the total number of missense variants is 10,
    this function will return 3.

    This function is used to determine the minimum size window to check for significant missense depletion
    when searching for two simultaneous breaks.

    :param str build: Reference genome build; must be one of BUILDS.
    :param bool get_total_exome_bases: Boolean for whether to recalculate total number of bases in exome.
        If False, will use value from `TOTAL_EXOME_BASES`. Default is False.
    :param bool get_total_gnomad_missense: Boolean for whether to recount total number of missense variants in gnomAD.
        If False, will use value from `TOTAL_GNOMAD_MISSENSE`. Default is False.
    :return: Average number of bases between observed missense variants, rounded to the nearest integer,
    :rtype: int
    """
    total_variants = TOTAL_GNOMAD_MISSENSE[CURRENT_VERSION]
    total_bases = TOTAL_EXOME_BASES[build]

    if get_total_exome_bases:
        if build not in BUILDS:
            raise DataException(f"Build must be one of {BUILDS}.")

        logger.info(
            "Getting total number of bases in the exome from full context HT..."
        )
        total_bases = get_exome_bases(build=build)

    if get_total_gnomad_missense:
        ht = filtered_exomes.ht()
        logger.info("Getting total number of missense variants in gnomAD...")
        total_variants = ht.count()

    logger.info("Total number of bases in the exome: %i", total_bases)
    logger.info(
        "Total number of missense variants in gnomAD exomes: %i", total_variants
    )
    logger.info("Getting average bases between missense variants and returning...")
    return round(total_bases / total_variants)


def keep_criteria(
    ac_expr: hl.expr.Int32Expression,
    af_expr: hl.expr.Float64Expression,
    filters_expr: hl.expr.SetExpression,
    cov_expr: hl.expr.Int32Expression,
    af_threshold: float = 0.001,
    cov_threshold: int = 0,
) -> hl.expr.BooleanExpression:
    """
    Return Boolean expression to filter variants in input Table.

    Default values will filter to rare variants (AC > 0, AF < 0.001) that pass filters and have median coverage greater than 0.

    :param hl.expr.Int32Expression ac_expr: Allele count Int32Expression.
    :param hl.expr.Float64Expression af_expr: Allele frequency (AF) Float64Expression.
    :param hl.expr.SetExpression filters_expr: Filters SetExpression.
    :param hl.expr.Int32Expression cov_expr: gnomAD median coverage Int32Expression.
    :param float af_threshold: Remove rows above this AF threshold. Default is 0.001.
    :param int cov_threshold: Remove rows below this median coverage threshold. Default is 0.
    :return: Boolean expression used to filter variants.
    """
    return (
        (ac_expr > 0)
        & (af_expr < af_threshold)
        & (hl.len(filters_expr) == 0)
        & (cov_expr > cov_threshold)
    )


def process_vep(ht: hl.Table, filter_csq: bool = False, csq: str = None) -> hl.Table:
    """
    Filter input Table to canonical transcripts only.

    Option to filter Table to specific variant consequence (csq).

    :param Table ht: Input Table.
    :param bool filter_csq: Whether to filter Table to a specific consequence. Default is False.
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
        if not csq:
            raise DataException("Need to specify consequence if filter_csq is True!")
        logger.info("Filtering to %s...", csq)
        ht = ht.filter(ht.transcript_consequences.most_severe_consequence == csq)
    return ht


def filter_to_region_type(ht: hl.Table, region: str) -> hl.Table:
    """
    Filter input Table to autosomes + chrX PAR, chrX non-PAR, or chrY non-PAR.

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
    Call `build_models` from gnomAD LoF repo to generate models used to adjust expected variants count.

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
    Get coverage correction for expected variants count.

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
    Get model to determine adjustment to mutation rate based on locus type and CpG status.

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


## Outlier transcript util
def get_constraint_transcripts(outlier: bool = True) -> hl.expr.SetExpression:
    """
    Read in LoF constraint HT results to get set of transcripts.

    Return either set of transcripts to keep (transcripts that passed transcript QC)
    or outlier transcripts.

    Transcripts are removed for the reasons detailed here:
    https://gnomad.broadinstitute.org/faq#why-are-constraint-metrics-missing-for-this-gene-or-annotated-with-a-note

    :param bool outlier: Whether to filter LoF constraint HT to outlier transcripts (if True),
        or QC-pass transcripts (if False). Default is True.
    :return: Set of outlier transcripts or transcript QC pass transcripts.
    :rtype: hl.expr.SetExpression
    """
    logger.warning(
        "Assumes LoF constraint has been separately calculated and that constraint HT exists..."
    )
    if not file_exists(constraint_ht.path):
        raise DataException("Constraint HT not found!")

    constraint_transcript_ht = constraint_ht.ht().key_by("transcript")
    constraint_transcript_ht = constraint_transcript_ht.filter(
        constraint_transcript_ht.canonical
    ).select("constraint_flag")
    if outlier:
        constraint_transcript_ht = constraint_transcript_ht.filter(
            hl.len(constraint_transcript_ht.constraint_flag) > 0
        )
    else:
        constraint_transcript_ht = constraint_transcript_ht.filter(
            hl.len(constraint_transcript_ht.constraint_flag) == 0
        )
    return hl.literal(
        constraint_transcript_ht.aggregate(
            hl.agg.collect_as_set(constraint_transcript_ht.transcript)
        )
    )


## Assessment utils
def import_clinvar_hi_variants(build: str, overwrite: bool) -> None:
    """
    Import ClinVar HT and filter to pathogenic/likely pathogenic missense variants in haploinsufficient genes.

    .. note::
        This function currently only works for build GRCh37.

    :param str build: Reference genome build; must be one of BUILDS.
    :param bool overwrite: Whether to overwrite ClinVar HT.
    :return: None; writes HT to resource path.
    """
    if build not in BUILDS:
        raise DataException(f"Build must be one of {BUILDS}.")
    if build == "GRCh37":
        from gnomad.resources.grch37.reference_data import clinvar

        clinvar_ht_path = grch37.clinvar_path_mis.path
    else:
        from gnomad.resources.grch38.reference_data import clinvar

        raise DataException("ClinVar files currently only exist for GRCh37!")

    if not file_exists(clinvar_ht_path) or overwrite:
        logger.info("Reading in ClinVar HT...")
        clinvar_ht = clinvar.ht()
        clinvar_ht = filter_to_clinvar_pathogenic(clinvar_ht)

        logger.info("Filtering to missense variants...")
        clinvar_ht = clinvar_ht.annotate(mc=clinvar_ht.info.MC)
        clinvar_ht = clinvar_ht.filter(
            clinvar_ht.mc.any(lambda x: x.contains("missense_variant"))
        )
        logger.info(
            "Number of variants after filtering to missense: %i", clinvar_ht.count()
        )

        logger.info("Filtering to variants in haploinsufficient genes...")
        # File header is '#gene'
        hi_ht = hl.import_table(hi_genes)
        hi_gene_set = hi_ht.aggregate(
            hl.agg.collect_as_set(hi_ht["#gene"]), _localize=False
        )

        logger.info("Getting gene information from ClinVar HT...")
        clinvar_ht = clinvar_ht.annotate(gene=clinvar_ht.info.GENEINFO.split(":")[0])
        clinvar_ht = clinvar_ht.filter(hi_gene_set.contains(clinvar_ht.gene))
        clinvar_ht = clinvar_ht.checkpoint(clinvar_ht_path, overwrite=overwrite)
        logger.info(
            "Number of variants after filtering to HI genes: %i", clinvar_ht.count()
        )


def import_de_novo_variants(
    build: str,
    overwrite: bool,
    csq_str: str = "consequence",
    missense_str: str = MISSENSE,
) -> None:
    """
    Import de novo missense variants.

    .. note::
        These files currently only exist for build GRCh37.

    :param str build: Reference genome build; must be one of BUILDS.
    :param bool overwrite: Whether to overwrite de novo Table.
    :param str csq_str: Name of field in de novo TSV that contains variant consequence.
        Default is 'consequence'.
    :param str missense_str: String representing missense variant effect. Default is MISSENSE.
    :return: None; writes HT to resource path.
    """
    if build not in BUILDS:
        raise DataException(f"Build must be one of {BUILDS}.")
    if build == "GRCh37":
        import grch37.de_novo_tsv as tsv_path
        import grch37.de_novo.path as ht_path
    else:
        raise DataException("De novo TSV does not exist for GRCh38!")

    dn_ht = hl.import_table(tsv_path, impute=True)
    dn_ht = dn_ht.filter(dn_ht[csq_str] == missense_str)
    dn_ht = dn_ht.transmute(
        locus=hl.locus(dn_ht.chrom, dn_ht.pos),
        alleles=[dn_ht.ref, dn_ht.alt],
    )
    # Key by locus and alleles for easy MPC annotation
    # (MPC annotation requires input Table to be keyed by locus and alleles)
    dn_ht = dn_ht.key_by("locus", "alleles")
    dn_ht.write(ht_path, overwrite=overwrite)
