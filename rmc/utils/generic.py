import logging
from typing import Dict, List, Tuple

import hail as hl
from gnomad.resources.grch37.gnomad import coverage, public_release
from gnomad.resources.grch37.reference_data import vep_context
from gnomad.resources.resource_utils import DataException
from gnomad.utils.constraint import (
    annotate_exploded_vep_for_constraint_groupings,
    annotate_with_mu,
    build_models,
)
from gnomad.utils.file_utils import file_exists
from gnomad.utils.filtering import filter_to_clinvar_pathogenic
from gnomad.utils.vep import (
    CSQ_NON_CODING,
    add_most_severe_csq_to_tc_within_vep_root,
    filter_vep_to_canonical_transcripts,
)
from gnomad_constraint.utils.constraint import prepare_ht_for_constraint_calculations

from rmc.resources.basics import (
    ACID_NAMES_PATH,
    CODON_TABLE_PATH,
    TEMP_PATH_WITH_FAST_DEL,
)
from rmc.resources.gnomad import constraint_ht, mutation_rate
from rmc.resources.reference_data import (
    autism_de_novo_2022_tsv_path,
    clinvar,
    clinvar_plp_mis_haplo,
    clinvar_plp_mis_triplo,
    dosage_ht,
    dosage_tsv_path,
    haplo_genes_path,
    ndd_de_novo,
    ndd_de_novo_2020_tsv_path,
    triplo_genes_path,
)
from rmc.resources.resource_utils import MISSENSE
from rmc.resources.rmc import DIVERGENCE_SCORES_TSV_PATH, MUTATION_RATE_TABLE_PATH

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_missense_constraint_generic")
logger.setLevel(logging.INFO)


## ExAC mutational model-related utils
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


## Codon and amino acid utils
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


def annotate_and_filter_codons(ht: hl.Table) -> hl.Table:
    """
    Remove non-coding loci and keep informative codons only.

    Split codon annotation to annotate reference and alternate amino acids, then
    remove rows with unknown amino acids.

    Additionally remove rows that are annotated as 'coding_sequence_variant',
    as these variants have either undefined or uninformative codon annotations
    (NA or codon with Ns, e.g. nnG/nnT).

    'coding_sequence_variant' defined as: 'A sequence variant that changes the coding sequence'
    https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html

    :param hl.Table ht: Input Table.
    :return: Table with informative codons only.
    """
    logger.info("Removing non-coding loci from HT...")
    non_coding_csq = hl.literal(CSQ_NON_CODING)
    ht = ht.filter(~non_coding_csq.contains(ht.most_severe_consequence))

    logger.info("Filtering to lines with expected codon annotations...")
    # Codons are in this format: NNN/NNN, so expected length is 7
    ht = ht.filter((hl.is_defined(ht.codons)) & (hl.len(ht.codons) == 7))
    codon_map = get_codon_lookup()
    ht = ht.annotate(
        ref=ht.codons.split("/")[0].upper(),
        alt=ht.codons.split("/")[1].upper(),
    )
    ht = ht.annotate(
        ref=codon_map.get(ht.ref, "Unk"),
        alt=codon_map.get(ht.alt, "Unk"),
    )
    # Remove any lines with "Unk" (unknown amino acids)
    return ht.filter((ht.ref != "Unk") & (ht.alt != "Unk"))


## Functions to process reference genome related resources
def process_context_ht(
    filter_to_missense: bool = True,
    missense_str: str = MISSENSE,
    add_annotations: bool = True,
) -> hl.Table:
    """
    Prepare context HT (SNPs only, annotated with VEP) for regional missense constraint calculations.

    Filter to missense variants in canonical protein coding transcripts.
    Also annotate with probability of mutation for each variant, CpG status, and methylation level.

    :param bool filter_to_missense: Whether to filter Table to missense variants only. Default is True.
    :param bool add_annotations: Whether to add ref, alt, methylation_level, exome_coverage, cpg, transition,
        and mutation_type annotations. Default is True.
    :return: Context HT filtered to canonical transcripts and optionally filtered to missense variants with
        mutation rate, CpG status, and methylation level annotations.
    :rtype: hl.Table
    """
    logger.info("Reading in SNPs-only, VEP-annotated context ht...")
    ht = vep_context.ht().select_globals()

    logger.info(
        "Filtering to canonical transcripts and annotating with most severe consequence...",
    )
    if filter_to_missense:
        ht = process_vep(ht, filter_csq=True, csq=missense_str)
    else:
        ht = process_vep(ht)

    if add_annotations:
        # `prepare_ht_for_constraint_calculations` annotates HT with:
        # ref, alt, methylation_level, exome_coverage, cpg, transition, mutation_type
        # NOTE: `mutation_type` was initially named `variant_type`, but this
        # field was renamed because `variant_type` refers to a different piece of
        # information in the gnomad_methods repo
        # See docstring for `annotate_mutation_type` for more details
        ht = prepare_ht_for_constraint_calculations(ht)

        logger.info("Annotating with mutation rate...")
        # Mutation rate HT is keyed by context, ref, alt, methylation level
        mu_ht = mutation_rate.ht().select("mu_snp")
        ht, grouping = annotate_exploded_vep_for_constraint_groupings(ht)
        ht = ht.select(
            "context",
            "ref",
            "alt",
            "methylation_level",
            "exome_coverage",
            "cpg",
            "transition",
            "mutation_type",
            *grouping,
        )
        return annotate_with_mu(ht, mu_ht)
    return ht


def get_gnomad_public_release(
    gnomad_data_type: str = "exomes", adj_freq_index: int = 0
) -> hl.Table:
    """
    Return gnomAD public sites Table annotated with coverage.

    Also filters HT to fields required by `keep_criteria`:
        - ac
        - af
        - filters
        - coverage

    :param gnomad_data_type: gnomAD data type. Used to retrieve public release and coverage resources.
        Must be one of "exomes" or "genomes" (check is done within `public_release`).
        Default is "exomes".
    :param adj_freq_index: Index of frequency array that contains global population filtered calculated on
        high quality (adj) genotypes. Default is 0.
    :return: gnomAD public sites HT annotated with coverage and filtered to select fields.
    """
    gnomad = public_release(gnomad_data_type).ht().select_globals()
    gnomad_cov = coverage(gnomad_data_type).ht()
    return gnomad.select(
        "filters",
        ac=gnomad.freq[adj_freq_index].AC,
        af=gnomad.freq[adj_freq_index].AF,
        gnomad_coverage=gnomad_cov[gnomad.locus].median,
    )


def filter_context_using_gnomad(
    context_ht: hl.Table,
    gnomad_data_type: str = "exomes",
    adj_freq_index: int = 0,
    filter_context_using_cov: bool = True,
    cov_threshold: int = 0,
) -> hl.Table:
    """
    Filter VEP context Table to sites that aren't seen in gnomAD or are rare in gnomAD.

    :param context_ht: VEP context Table.
    :param gnomad_data_type: gnomAD data type. Used to retrieve public release and coverage resources.
        Must be one of "exomes" or "genomes" (check is done within `public_release`).
        Default is "exomes".
    :param adj_freq_index: Index of frequency array that contains global population filtered calculated on
        high quality (adj) genotypes. Default is 0.
    :param filter_context_using_cov: Whether to also filter sites in context Table using gnomAD coverage.
        Default is True.
    :param cov_threshold: Coverage threshold used to filter context Table if `filter_context_using_cov` is True.
        Default is 0.
    :return: Filtered VEP context Table.
    """
    gnomad = get_gnomad_public_release(gnomad_data_type, adj_freq_index)

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
        gnomad_cov = coverage(gnomad_data_type).ht()
        context_ht = context_ht.annotate(
            gnomad_coverage=gnomad_cov[context_ht.locus].median
        )
        context_ht = context_ht.filter(context_ht.gnomad_coverage > cov_threshold)
        # Drop coverage annotation here for consistency
        # (This annotation is only added if `filter_context_using_cov` is True)
        context_ht = context_ht.drop("gnomad_coverage")

    return context_ht


def get_annotations_from_context_ht_vep(
    ht: hl.Table, annotations: List[str] = ["polyphen", "sift"]
) -> hl.Table:
    """
    Get desired annotations from context Table's VEP annotation and annotate/filter to informative amino acids.

    Function will get 'codons', 'most_severe_consequence', and 'transcript' annotations by default.
    Additional annotations should be specified using `annotations`.

    :param hl.Table ht: Input VEP context Table.
    :param List[str] annotations: Additional annotations to extract from context HT's VEP field.
        Default is ["polyphen", "sift"].
    :return: VEP context Table with rows filtered to remove loci that are non-coding and retain only those with informative amino acids
        and columns filtered to keep only specified annotations.
    """
    annot_expr = {
        "codons": ht.transcript_consequences.codons,
        "most_severe_consequence": ht.transcript_consequences.most_severe_consequence,
        "transcript": ht.transcript_consequences.transcript_id,
    }
    if "polyphen" in annotations:
        annot_expr["polyphen"] = hl.struct(
            prediction=ht.transcript_consequences.polyphen_prediction,
            score=ht.transcript_consequences.polyphen_score,
        )
    if "sift" in annotations:
        annot_expr["sift"] = hl.struct(
            prediction=ht.transcript_consequences.sift_prediction,
            score=ht.transcript_consequences.sift_score,
        )
    ht = ht.annotate(**annot_expr)
    annotations.extend(["codons", "most_severe_consequence", "transcript"])
    ht = ht.select(*annotations)
    return annotate_and_filter_codons(ht)


## Functions for obs/exp related resources
def get_exome_bases() -> int:
    """
    Get number of bases in the exome.

    Read in context HT (containing all coding bases in the genome), remove outlier transcripts, and filter to median coverage >= 5.

    :return: Number of bases in the exome.
    :rtype: int
    """
    logger.info("Reading in SNPs-only, VEP-annotated context ht...")
    ht = vep_context.ht()

    logger.info("Filtering to canonical transcripts...")
    ht = filter_vep_to_canonical_transcripts(ht, vep_root="vep", filter_empty_csq=True)

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


def keep_criteria(
    ac_expr: hl.expr.Int32Expression,
    af_expr: hl.expr.Float64Expression,
    filters_expr: hl.expr.SetExpression,
    cov_expr: hl.expr.Int32Expression,
    af_threshold: float = 0.001,
    cov_threshold: int = 0,
    filter_to_rare: bool = True,
) -> hl.expr.BooleanExpression:
    """
    Return Boolean expression to filter variants in input Table.

    Default values will filter to rare variants (AC > 0, AF < 0.001) that pass filters and have median coverage greater than 0.

    :param ac_expr: Allele count Int32Expression.
    :param af_expr: Allele frequency (AF) Float64Expression.
    :param filters_expr: Filters SetExpression.
    :param cov_expr: gnomAD median coverage Int32Expression.
    :param af_threshold: Remove rows above this AF threshold. Default is 0.001.
    :param cov_threshold: Remove rows below this median coverage threshold. Default is 0.
    :param filter_to_rare: Whether to filter to keep rare variants only. Default is True.
    :return: Boolean expression used to filter variants.
    """
    af_filter_expr = (
        (af_expr < af_threshold) if filter_to_rare else (af_expr > af_threshold)
    )
    return (
        (ac_expr > 0)
        & (af_filter_expr)
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
    ht = filter_vep_to_canonical_transcripts(ht, vep_root="vep", filter_empty_csq=True)

    logger.info("Annotating HT with most severe consequence...")
    ht = add_most_severe_csq_to_tc_within_vep_root(ht)
    ht = ht.transmute(transcript_consequences=ht.vep.transcript_consequences)
    ht = ht.explode(ht.transcript_consequences)

    logger.info("Filtering to non-outlier transcripts...")
    # Keep transcripts used in LoF constraint only (remove all other outlier transcripts)
    constraint_transcripts = get_constraint_transcripts(outlier=False)
    ht = ht.filter(
        constraint_transcripts.contains(ht.transcript_consequences.transcript_id)
    )

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
    :param bool weighted: Whether to use weighted least squares when building models. Default is True.
    :return: Coverage model, plateau models for autosomes, plateau models for chrX, plateau models for chrY.
    :rtype: Tuple[Tuple[float, float], hl.expr.DictExpression, hl.expr.DictExpression, hl.expr.DictExpression]
    """
    logger.info("Building autosomes/PAR plateau model and coverage model...")
    coverage_model, plateau_models = build_models(
        coverage_ht,
        weighted=weighted,
    )

    logger.info("Building plateau models for chrX and chrY...")
    # TODO: make half_cutoff (for coverage cutoff) True for X/Y?
    # This would also mean saving new coverage model for allosomes
    _, plateau_x_models = build_models(coverage_x_ht, weighted=weighted)
    _, plateau_y_models = build_models(coverage_y_ht, weighted=weighted)
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
        # Use log10 here per
        # https://github.com/broadinstitute/gnomad_lof/blob/master/constraint_utils/constraint_basics.py#L555
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
def import_dosage(
    overwrite: bool,
    haplo_threshold: float = 0.86,
    triplo_threshold: float = 0.94,
) -> None:
    """
    Import gene dosage sensitivity information.

    Also create HailExpressions with haploinsufficient
    and triplosensitive genes.

    :param overwrite: Whether to overwrite output data.
    :param haplo_threshold: pHaplo score threshold for determining whether a gene is predicted haploinsufficient.
        Default is 0.86 (from Collins et al. paper).
    :param triplo_threshold: pTriplo score threshold for determining whether a gene is predicted triplosensitive.
        Default is 0.94 (from Collins et al. paper).
    :return: None; function writes data to resource paths.
    """
    ht = hl.import_table(dosage_tsv_path, impute=True, force=True)
    ht = ht.transmute(gene=ht["#gene"])
    ht = ht.key_by("gene")
    ht = ht.annotate_globals(
        haplo_cutoff=haplo_threshold,
        triplo_cutoff=triplo_threshold,
    )
    ht = ht.checkpoint(
        dosage_ht.path, _read_if_exists=not overwrite, overwrite=overwrite
    )

    haplo_genes = ht.filter(ht.pHaplo >= haplo_threshold)
    haplo_genes = haplo_genes.aggregate(hl.agg.collect_as_set(haplo_genes.gene))
    hl.experimental.write_expression(
        haplo_genes,
        haplo_genes_path,
        overwrite=overwrite,
    )
    triplo_genes = ht.filter(ht.pTriplo >= triplo_threshold)
    triplo_genes = triplo_genes.aggregate(hl.agg.collect_as_set(triplo_genes.gene))
    hl.experimental.write_expression(
        triplo_genes,
        triplo_genes_path,
        overwrite=overwrite,
    )


def import_clinvar(overwrite: bool, missense_str: str = MISSENSE) -> None:
    """
    Import ClinVar HT and pathogenic/likely pathogenic missense variants.

    Also filter P/LP missense variants to variants in haploinsufficient
    and triplosensitive genes.

    .. note::
        This function currently only works for build GRCh37.

    :param bool overwrite: Whether to overwrite output data.
    :param missense_str: String that corresponds to missense variant consequence.
        Default is MISSENSE.
    :return: None; writes HTs and HEs to resource paths.
    """
    if (
        not file_exists(clinvar_plp_mis_haplo.path)
        or not file_exists(clinvar_plp_mis_triplo.path)
        or overwrite
    ):
        logger.info("Reading in ClinVar HT...")
        ht = clinvar.ht()
        logger.info("Filtering to P/LP missense variants...")
        ht = filter_to_clinvar_pathogenic(ht)
        ht = ht.annotate(mc=ht.info.MC)
        ht = ht.filter(ht.mc.any(lambda x: x.contains(missense_str)))

        logger.info("Getting gene information from ClinVar HT...")
        ht = ht.annotate(gene=ht.info.GENEINFO.split(":")[0])
        ht = ht.checkpoint(f"{TEMP_PATH_WITH_FAST_DEL}/clinvar.ht", overwrite=True)
        logger.info(
            "Number of variants after filtering to P/LP missense: %i", ht.count()
        )

        logger.info("Filtering to variants in haploinsufficient genes...")
        if not file_exists(haplo_genes_path):
            import_dosage(overwrite)
        hi_genes = hl.experimental.read_expression(haplo_genes_path)
        haplo_ht = ht.filter(hi_genes.contains(ht.gene))
        haplo_ht = haplo_ht.checkpoint(
            clinvar_plp_mis_haplo.path,
            _read_if_exists=not overwrite,
            overwrite=overwrite,
        )
        logger.info(
            "Number of variants after filtering to HI genes: %i", haplo_ht.count()
        )

        logger.info("Filtering to variants in triplosensitive genes...")
        triplo_genes = hl.experimental.read_expression(triplo_genes_path)
        triplo_ht = ht.filter(triplo_genes.contains(ht.gene))
        triplo_ht = triplo_ht.checkpoint(
            clinvar_plp_mis_triplo.path,
            _read_if_exists=not overwrite,
            overwrite=overwrite,
        )
        logger.info(
            "Number of variants after filtering to TS genes: %i",
            triplo_ht.count(),
        )


def import_fu_data(overwrite: bool) -> None:
    """
    Import de novo variants from Fu et al. (2022) paper.

    Function imports variants from TSV into HT, removes malformed rows,
    and lifts data from GRCh38 to GRCh37.

    :param overwrite: Whether to overwrite Table if it exists.
    :return: None; Function writes Table to temporary path.
    """
    fu_ht = hl.import_table(
        autism_de_novo_2022_tsv_path,
        impute=True,
        # Skip blank lines at the bottom of this TSV
        missing="",
        skip_blank_lines=True,
    )
    # Remove lines from bottom of TSV that are parsed incorrectly upon import
    # These lines contain metadata about the TSV, e.g.:
    # "Supplementary Table 20. The de novo SNV/indel variants used in TADA
    # association analyses from assembled ASD cohorts"
    fu_ht = fu_ht.filter(~hl.is_missing(fu_ht.Role))

    # Prepare to lift dataset back to GRCh37
    rg37 = hl.get_reference("GRCh37")
    rg38 = hl.get_reference("GRCh38")
    rg38.add_liftover(
        "gs://hail-common/references/grch38_to_grch37.over.chain.gz", rg37
    )

    fu_ht = fu_ht.annotate(
        locus=hl.parse_locus(
            hl.format(
                "chr%s:%s",
                fu_ht.Variant.split(":")[0],
                fu_ht.Variant.split(":")[1],
            ),
            reference_genome="GRCh38",
        ),
        alleles=[fu_ht.Variant.split(":")[2], fu_ht.Variant.split(":")[3]],
    )
    fu_ht = fu_ht.annotate(
        new_locus=hl.liftover(fu_ht.locus, "GRCh37", include_strand=True),
        old_locus=fu_ht.locus,
    )
    liftover_stats = fu_ht.aggregate(
        hl.struct(
            liftover_fail=hl.agg.count_where(hl.is_missing(fu_ht.new_locus)),
            flip_strand=hl.agg.count_where(fu_ht.new_locus.is_negative_strand),
        )
    )
    logger.info(
        "%i variants failed to liftover, and %i variants flipped strands",
        liftover_stats.liftover_fail,
        liftover_stats.flip_strand,
    )
    fu_ht = fu_ht.filter(
        hl.is_defined(fu_ht.new_locus) & ~fu_ht.new_locus.is_negative_strand
    )
    fu_ht = fu_ht.key_by(locus=fu_ht.new_locus.result, alleles=fu_ht.alleles)

    # Rename 'Proband' > 'ASD' and 'Sibling' > 'control'
    fu_ht = fu_ht.transmute(role=hl.if_else(fu_ht.Role == "Proband", "ASD", "control"))
    fu_ht = fu_ht.group_by("locus", "alleles").aggregate(
        role=hl.agg.collect(fu_ht.role),
    )
    fu_ht.write(f"{TEMP_PATH_WITH_FAST_DEL}/fu_dn.ht", overwrite=overwrite)


def import_kaplanis_data(overwrite: bool) -> None:
    """
    Import de novo variants from Kaplanis et al. (2020) paper.

    Function imports variants from TSV into HT, and filters to cases ascertained for
    developmental delay/intellectual disability only.
    Input TSV also contains autistic individuals, but these individuals overlap with
    data from Fu et al. paper and are therefore not retained.

    :param overwrite: Whether to overwrite Table if it exists.
    :return: None; Function writes Table to temporary path.
    """
    kap_ht = hl.import_table(ndd_de_novo_2020_tsv_path, impute=True)
    kap_ht = kap_ht.transmute(
        locus=hl.locus(kap_ht.chrom, kap_ht.pos),
        alleles=[kap_ht.ref, kap_ht.alt],
    )
    kap_ht = kap_ht.filter(hl.is_snp(kap_ht.alleles[0], kap_ht.alleles[1]))
    kap_ht = kap_ht.filter(kap_ht.case_control == "DD")
    kap_ht = kap_ht.group_by("locus", "alleles").aggregate(
        case_control=hl.agg.collect(kap_ht.case_control)
    )
    kap_ht.write(f"{TEMP_PATH_WITH_FAST_DEL}/kaplanis_dn.ht", overwrite=overwrite)


def import_de_novo_variants(overwrite: bool, n_partitions: int = 5000) -> None:
    """
    Import de novo missense variants.

    .. note::
        These files currently only exist for build GRCh37.

    :param bool overwrite: Whether to overwrite de novo Table.
    :paran n_partitions: Number of partitions for input Tables.
        Used to repartition Tables on read.
        Will also help determine number of partitions in final Table.
    :return: None; writes HT to resource path.
    """
    fu_ht_path = f"{TEMP_PATH_WITH_FAST_DEL}/fu_dn.ht"
    kaplanis_ht_path = f"{TEMP_PATH_WITH_FAST_DEL}/kaplanis_dn.ht"
    if not file_exists(fu_ht_path) or overwrite:
        import_fu_data(overwrite=overwrite)
    if not file_exists(kaplanis_ht_path) or overwrite:
        import_kaplanis_data(overwrite=overwrite)

    fu_ht = hl.read_table(fu_ht_path, _n_partitions=n_partitions)
    kap_ht = hl.read_table(kaplanis_ht_path, _n_partitions=n_partitions)
    ht = kap_ht.join(fu_ht, how="outer")

    # Union sample types (DD, ASD, control)
    ht = ht.transmute(
        sample_set=hl.case()
        .when(
            hl.is_defined(ht.case_control) & hl.is_defined(ht.role),
            ht.case_control.extend(ht.role),
        )
        .when(hl.is_defined(ht.case_control), ht.case_control)
        .when(hl.is_defined(ht.role), ht.role)
        .or_missing()
    )
    ht.write(ndd_de_novo.path, overwrite=overwrite)
