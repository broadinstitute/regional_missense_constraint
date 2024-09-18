import logging
from typing import Dict, List, Set, Tuple

import hail as hl
from gnomad.resources.grch38.gnomad import public_release
from gnomad.resources.grch38.reference_data import vep_context
from gnomad.resources.resource_utils import DataException
from gnomad.utils.constraint import annotate_exploded_vep_for_constraint_groupings
from gnomad.utils.file_utils import file_exists
from gnomad.utils.vep import (
    CSQ_NON_CODING,
    explode_by_vep_annotation,
    filter_vep_transcript_csqs,
    process_consequences,
)
from gnomad_constraint.resources.resource_utils import get_preprocessed_ht

from rmc.resources.basics import (
    ACID_NAMES_PATH,
    CODON_TABLE_PATH,
    TEMP_PATH_WITH_FAST_DEL,
)
from rmc.resources.gnomad import constraint_ht
from rmc.resources.reference_data import VEP_VERSION

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_missense_constraint_generic")
logger.setLevel(logging.INFO)


####################################################################################
## Codon and amino acid utils
####################################################################################
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


####################################################################################
## Reference genome processing-related utils
####################################################################################
def process_context_ht(
    filter_to_canonical: bool = False,
    filter_csq: Set[str] = None,
) -> hl.Table:
    """
    Get context HT for SNPs annotated with VEP in canonical protein-coding transcripts.

    This function offers options to filter to specific variant consequences and add annotations
    to prepare for regional missense constraint calculations.

    :param filter_to_canonical: Whether to filter to canonical transcripts only. Default is False.
    :param filter_csq: Specific consequences to keep. Default is None.
    :return: VEP context HT filtered to canonical transcripts and optionally filtered to variants
        in non-outlier transcripts with specific consequences and annotated with mutation rate etc.
    :rtype: hl.Table
    """
    logger.info("Reading in gene constraint preprocessed context ht...")
    # Created using
    # https://github.com/broadinstitute/gnomad-constraint/blob/0acd2815e59c04d642bb705e6d1ca166f5d79e5f/gnomad_constraint/utils/constraint.py#L78
    # NOTE: The v4.1 constraint table currently only contains autosomes
    # TODO: Add allosomes and PAR
    # Used this table for testing
    # ht = hl.read_table(
    #    "gs://gnomad/v4.1/constraint_an/preprocessed_data/gnomad.v4.1.context.preprocessed.autosome_par.ht"
    # ).select_globals()
    ht = get_preprocessed_ht("context").ht().select_globals()

    if filter_to_canonical:
        logger.info("Filtering to canonical transcripts only...")

    # Filter to protein-coding ENST transcripts
    # Also optionally filter to canonical transcripts and specific consequences
    ht = filter_vep_transcript_csqs(
        t=ht,
        vep_root="vep",
        synonymous=False,
        canonical=filter_to_canonical,
        protein_coding=True,
        ensembl_only=True,
        filter_empty_csq=True,
        csqs=filter_csq,
    )

    logger.info(
        "Annotating context HT with annotations needed for constraint groupings..."
    )
    # See https://github.com/broadinstitute/gnomad-constraint/blob/0acd2815e59c04d642bb705e6d1ca166f5d79e5f/gnomad_constraint/utils/constraint.py#L78
    # and
    # https://github.com/broadinstitute/gnomad_methods/blob/2387430a79068225114c620952275f2f805cb24a/gnomad/utils/constraint.py#L948
    # for documentation on annotations added
    ht, _ = annotate_exploded_vep_for_constraint_groupings(
        ht=ht,
        coverage_expr=ht.exomes_AN,
        vep_annotation="transcript_consequences",
        include_canonical_group=True,
        # NOTE: all canonical transcripts are also the MANE select transcript in
        # GENCODE v39/VEP v105 context HT
        include_mane_select_group=False,
    )

    # Drop other unnecessary annotations
    return ht.select(
        "context",
        "ref",
        "alt",
        "methylation_level",
        "cpg",
        "mutation_type",
        "annotation",
        "modifier",
        "coverage",
        "transcript",
        "exomes_AN",
        "exomes_AN_percent",
    )


def get_aa_from_context(
    overwrite_temp: bool,
    keep_transcripts: Set[str] = None,
    n_partitions: int = 10000,
    filter_to_canonical: bool = False,
    vep_version: str = VEP_VERSION,
) -> hl.Table:
    """
    Extract amino acid information from VEP context HT.

    :param overwrite_temp: Whether to overwrite temporary data.
        If False, will read existing temp data rather than overwriting.
        If True, will overwrite temp data.
    :param keep_transcripts: Desired set of transcripts to keep from the context HT.
        If set, function will filter to keep these trancripts only.
        Default is None.
    :param n_partitions: Desired number of partitions for context HT after filtering.
        Default is 10,000.
    :param filter_to_canonical: Whether to filter to canonical transcripts only. Default is False.
    :param vep_version: VEP version to use. Default is `VEP_VERSION`.
    :return: VEP context HT filtered to keep only transcript ID, protein number, and amino acid information.
    """
    logger.info(
        "Reading in VEP context HT and filtering to coding effects in canonical"
        " transcripts..."
    )
    # TODO: Add option to filter to non-outliers if still desired
    # Drop globals and select only VEP transcript consequences field
    ht = (
        hl.read_table(vep_context.versions[vep_version].path)
        .select_globals()
        .select("vep", "was_split")
    )
    ht = ht.naive_coalesce(n_partitions)
    ht = filter_vep_transcript_csqs(
        t=ht,
        vep_root="vep",
        synonymous=False,
        canonical=filter_to_canonical,
        protein_coding=True,
        ensembl_only=True,
        filter_empty_csq=True,
    )
    ht = process_consequences(ht)
    ht = explode_by_vep_annotation(ht, "transcript_consequences")
    ht = ht.select("transcript_consequences")
    ht = ht.filter(
        ~hl.literal(CSQ_NON_CODING).contains(
            ht.transcript_consequences.most_severe_consequence
        )
    )
    # Unnest annotations from context HT
    ht = ht.transmute(
        transcript=ht.transcript_consequences.transcript_id,
        aa_start_num=ht.transcript_consequences.protein_start,
        aa_end_num=ht.transcript_consequences.protein_end,
        amino_acids=ht.transcript_consequences.amino_acids,
    )

    if keep_transcripts:
        logger.info("Filtering to desired transcripts only...")
        ht = ht.filter(hl.literal(keep_transcripts).contains(ht.transcript))

    ht = ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/vep_amino_acids.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )
    return ht


def get_ref_aa(
    ht: hl.Table,
    overwrite_temp: bool,
    extra_aa_map: Dict[str, str] = {"*": "Ter", "U": "Sec"},
    aa_to_remove: Set[str] = {"X"},
) -> hl.Table:
    """
    Filter input HT to keep only reference amino acid information (identity and number).

    :param ht: Input Table. Should be HT output from `get_aa_from_context`.
    :param overwrite_temp: Whether to overwrite temporary data.
        If False, will read existing temp data rather than overwriting.
        If True, will overwrite temp data.
    :param extra_aa_map: Dictionary mapping any amino acid one letter codes to three letter codes.
        Designed to capture any amino acids not present in `ACID_NAMES_PATH`.
        Default is {"*": "Ter", "U": "Sec"}.
    :param aa_to_remove: Any amino acid one letter codes to remove. Default is {"X"}.
    :return: HT with reference AA identity and number annotated for each locus-transcript combination.
    """
    logger.info("Mapping amino acid one letter codes to three letter codes...")
    aa_map = get_aa_map()
    # Add any extra mappings into amino acid mapping dict
    if extra_aa_map:
        aa_map.update(extra_aa_map)
    aa_map = hl.literal(aa_map)
    ht = ht.annotate(ref_aa_1_letter=ht.amino_acids.split("/")[0])
    ht = ht.annotate(ref_aa=aa_map.get(ht.ref_aa_1_letter, ht.ref_aa_1_letter))
    if aa_to_remove:
        ht = ht.filter(~hl.literal(aa_to_remove).contains(ht.ref_aa))

    # Select fields and checkpoint
    ht = ht.select("ref_aa", "aa_start_num", "aa_end_num", "transcript")
    ht = ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/ref_amino_acids.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    # Check if there are any ref amino acids in HT that aren't in `aa_map`
    ref_aa_check_he_path = f"{TEMP_PATH_WITH_FAST_DEL}/ref_aa_3_letter_check.he"
    overwrite_ref_aa_he = (
        not file_exists(ref_aa_check_he_path) if not overwrite_temp else overwrite_temp
    )
    if overwrite_ref_aa_he:
        ref_aa_check = ht.aggregate(hl.agg.collect_as_set(ht.ref_aa))
        hl.experimental.write_expression(
            ref_aa_check, ref_aa_check_he_path, overwrite=True
        )
    ref_aa_check = hl.eval(hl.experimental.read_expression(ref_aa_check_he_path))
    ref_aa_check = ref_aa_check.difference(set(hl.eval(aa_map).values()))
    if len(ref_aa_check) != 0:
        logger.warning(
            "The following reference amino acids were not mapped to three letter"
            " codes: %s",
            ref_aa_check,
        )

    # Double check that protein start always equals protein end
    protein_num_check_he_path = f"{TEMP_PATH_WITH_FAST_DEL}/protein_num_count.he"
    overwrite_protein_num_he = (
        not file_exists(protein_num_check_he_path)
        if not overwrite_temp
        else overwrite_temp
    )
    if overwrite_protein_num_he:
        protein_num_check = ht.aggregate(
            hl.agg.count_where(ht.aa_start_num != ht.aa_end_num)
        )
        hl.experimental.write_expression(
            protein_num_check, protein_num_check_he_path, overwrite=True
        )
    # Assume file already exists otherwise
    protein_num_check = hl.eval(
        hl.experimental.read_expression(protein_num_check_he_path)
    )
    if protein_num_check != 0:
        raise DataException(
            f"{protein_num_check} sites had different amino acid numbers at start and"
            " end -- please double check!"
        )

    # Reformat reference AA to have both the 3 letter code and number
    ht = ht.annotate(
        ref_aa=hl.or_missing(
            hl.is_defined(ht.ref_aa) & hl.is_defined(ht.aa_start_num),
            hl.format("%s%s", ht.ref_aa, ht.aa_start_num),
        )
    )

    # Collect by key to collapse AA per locus
    ht = ht.key_by("locus", "transcript").drop("alleles", "aa_start_num", "aa_end_num")
    ht = ht.collect_by_key(name="aa_info")
    ht = ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/ref_aa_collected.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    # Check to see if AA info is defined for all alleles associated with a locus/transcript
    missing_aa_check_he_path = f"{TEMP_PATH_WITH_FAST_DEL}/missing_aa_check.he"
    overwrite_he = (
        not file_exists(missing_aa_check_he_path)
        if not overwrite_temp
        else overwrite_temp
    )
    if overwrite_he:
        ht = ht.annotate(
            any_aa_missing=hl.any(hl.map(lambda x: hl.is_missing(x.ref_aa), ht.aa_info))
        )
        missing_aa_check = ht.aggregate(hl.agg.count_where(ht.any_aa_missing))
        missing_aa_check = hl.experimental.write_expression(
            missing_aa_check, missing_aa_check_he_path, overwrite=True
        )

    missing_aa_check = hl.eval(
        hl.experimental.read_expression(missing_aa_check_he_path)
    )
    if missing_aa_check != 0:
        logger.warning(
            "%i locus-transcript combinations had missing AA info for at least 1"
            " allele!",
            missing_aa_check,
        )
    return ht.transmute(ref_aa=ht.aa_info[0].ref_aa)


####################################################################################
## gnomAD and context HT processing-related utils
####################################################################################
def get_gnomad_public_release(
    gnomad_data_type: str = "exomes", adj_freq_index: int = 0
) -> hl.Table:
    """
    Return gnomAD public sites Table annotated with coverage.

    Also filters HT to fields required by `keep_criteria`:
        - ac
        - af
        - filters

    :param gnomad_data_type: gnomAD data type. Used to retrieve public release resource.
        Must be one of "exomes" or "genomes" (check is done within `public_release`).
        Default is "exomes".
    :param adj_freq_index: Index of array that contains allele frequency information calculated on
        high quality (adj) genotypes across genetic ancestry groups. Default is 0.
    :return: gnomAD public sites HT annotated with coverage and filtered to select fields.
    """
    gnomad = public_release(gnomad_data_type).ht().select_globals()
    return gnomad.select(
        "filters",
        ac=gnomad.freq[adj_freq_index].AC,
        af=gnomad.freq[adj_freq_index].AF,
    )


def filter_context_using_gnomad(
    context_ht: hl.Table,
    gnomad_data_type: str = "exomes",
    adj_freq_index: int = 0,
    an_threshold: int = 0,
) -> hl.Table:
    """
    Filter VEP context Table to sites that aren't seen in gnomAD or are rare in gnomAD.

    Also filter sites with zero coverage in gnomAD.

    :param context_ht: VEP context Table.
    :param gnomad_data_type: gnomAD data type. Used to retrieve public release resource.
        Must be one of "exomes" or "genomes" (check is done within `public_release`).
        Default is "exomes".
    :param adj_freq_index: Index of array that contains allele frequency information calculated on
        high quality (adj) genotypes across genetic ancestry groups. Default is 0.
    :param an_threshold: Remove variants at or below this AN threshold (in gnomAD exomes). Default is 0.
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
        )
    )
    return context_ht.filter(
        hl.is_defined(context_ht.exomes_AN) & (context_ht.exomes_AN > an_threshold)
    )


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


def keep_criteria(
    ac_expr: hl.expr.Int32Expression,
    af_expr: hl.expr.Float64Expression,
    filters_expr: hl.expr.SetExpression,
    af_threshold: float = 0.001,
    filter_to_rare: bool = True,
) -> hl.expr.BooleanExpression:
    """
    Return Boolean expression to filter variants in input Table.

    Default values will filter to rare variants (AC > 0, AF < 0.001) that pass filters.

    :param ac_expr: Allele count (AC) Int32Expression.
    :param af_expr: Allele frequency (AF) Float64Expression.
    :param filters_expr: Filters SetExpression.
    :param af_threshold: AF threshold used for filtering variants in combination with `filter_to_rare`. Default is 0.001.
    :param filter_to_rare: Whether to filter to keep rare variants only.
        If True, only variants with AF < `af_threshold` will be kept.
        If False, only variants with AF >= `af_threshold` will be kept.
        Default is True.
    :return: Boolean expression used to filter variants.
    """
    af_filter_expr = (
        (af_expr < af_threshold) if filter_to_rare else (af_expr >= af_threshold)
    )
    return (ac_expr > 0) & (af_filter_expr) & (hl.len(filters_expr) == 0)


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


def get_coverage_correction_expr(
    an_expr: hl.expr.Int32Expression,
    coverage_model: Tuple[float, float],
    high_AN_cutoff: int = 90,
) -> hl.expr.Float64Expression:
    """
    Get 'coverage' correction for expected variants count.

    .. note::
        - As of gnomAD v4, we use allele number (AN) as a proxy for coverage,
            because coverage in v4 was not computed from CRAMs.
        - Default high coverage cutoff taken from gnomAD LoF repo.

    :param ht: Input AN expression. Should be percent of exome AN defined at each locus.
    :param coverage_model: Model to determine coverage correction factor necessary
         for calculating expected variants at low AN sites.
    :param high_AN_cutoff: Cutoff for high AN value. Default is 90.
        NOTE that default should be adjusted based on upstream gene constraint pipeline default.
    :return: Coverage correction expression.
    :rtype: hl.expr.Float64Expression
    """
    return (
        hl.case()
        .when(an_expr == 0, 0)
        .when(an_expr >= high_AN_cutoff, 1)
        # Use log10 here per
        # https://github.com/broadinstitute/gnomad_lof/blob/master/constraint_utils/constraint_basics.py#L555
        # .default(coverage_model[1] * hl.log10(an_expr) + coverage_model[0])
        .default(coverage_model[1] * an_expr + coverage_model[0])
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


####################################################################################
## Outlier transcript utils
####################################################################################
def get_constraint_transcripts(
    all_transcripts: bool = False,
    filter_to_canonical: bool = False,
    outlier: bool = True,
) -> hl.expr.SetExpression:
    """
    Read in LoF constraint HT results to get set of transcripts.

    Return either set of transcripts to keep (transcripts that passed transcript QC)
    or outlier transcripts.

    Transcripts are removed for the reasons detailed here:
    https://gnomad.broadinstitute.org/faq#why-are-constraint-metrics-missing-for-this-gene-or-annotated-with-a-note

    .. note::
        - Function assumes that LoF constraint HT has been filtered to include only
            protein-coding transcripts.

    :param all_transcripts: Whether to filter to all transcripts. Will only keep
        all transcripts if `filter_to_canonical` is False, otherwise toggles
        between removing or keeping non-outlier transcripts. Default is False.
    :param filter_to_canonical: Whether to filter to canonical transcripts only. Default is False.
    :param outlier: Whether to filter LoF constraint HT to outlier transcripts (if True),
        or QC-pass transcripts (if False). Applies only if `all_transcripts` is False.
        Default is True.
    :return: Set of outlier transcripts or transcript QC pass transcripts.
    :rtype: hl.expr.SetExpression
    """
    logger.warning(
        "Assumes LoF constraint has been separately calculated and that constraint HT"
        " exists..."
    )
    if not file_exists(constraint_ht.path):
        raise DataException("Constraint HT not found!")

    constraint_transcript_ht = constraint_ht.ht().key_by("transcript")
    # NOTE: all protein-coding transcripts are ENST transcripts in constraint HT
    constraint_transcript_ht = constraint_transcript_ht.filter(
        constraint_transcript_ht.transcript_type == "protein_coding"
    )
    if filter_to_canonical:
        constraint_transcript_ht = constraint_transcript_ht.filter(
            constraint_transcript_ht.canonical
        )

    if not all_transcripts:
        constraint_transcript_ht = constraint_transcript_ht.select("constraint_flags")
        if outlier:
            constraint_transcript_ht = constraint_transcript_ht.filter(
                hl.len(constraint_transcript_ht.constraint_flags) > 0
            )
        else:
            constraint_transcript_ht = constraint_transcript_ht.filter(
                hl.len(constraint_transcript_ht.constraint_flags) == 0
            )

    return hl.literal(
        constraint_transcript_ht.aggregate(
            hl.agg.collect_as_set(constraint_transcript_ht.transcript)
        )
    )
