import logging
from typing import Dict

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.vep import filter_vep_to_canonical_transcripts
from gnomad_lof.constraint.constraint_basics import (
    add_most_severe_csq_to_tc_within_ht,
    prepare_ht,
)
from rmc.resources.basics import (
    ACID_NAMES_PATH,
    CODON_TABLE_PATH,
    divergence_scores,
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
    build: str, trimers: bool = True, overwrite: bool = True, n_partitions: int = 1000
) -> None:
    """
    Imports reference fasta (SNPs only, VEP'd) as ht
    Filters to canonical protein coding transcripts

    :param str build: Reference genome build; must be one of BUILDS.
    :param bool trimers: Whether to filter to trimers or heptamers.
    :param bool overwrite: Whether to overwrite output.
    :param int n_partitions: Number of desired partitions for output. Default is 1000.
    :return: None
    :rtype: None
    """
    if build not in BUILDS:
        raise DataException(f"Build must be one of {BUILDS}.")

    logger.info("Reading in SNPs-only, VEP-annotated context ht...")

    if build == "GRCh37":
        full_context_ht = grch37.full_context.ht()
        output_path = grch37.processed_context.path
    else:
        full_context_ht = grch38.full_context.ht()
        output_path = grch38.processed_context.path

    full_context_ht = prepare_ht(full_context_ht, trimers)

    logger.info(
        "Filtering to missense variants in canonical protein coding transcripts..."
    )
    context_ht = filter_vep_to_canonical_transcripts(full_context_ht)

    # NOTE: Had issues with filtering using most_severe_consequence in nb, so switched code here
    # Also discovered that transcript_consequences always has a length of 1 in notebook (at least in chr22)
    context_ht = context_ht.filter(
        context_ht.vep.transcript_consequences[0].consequence_terms.contains(
            "missense_variant"
        )
    )

    logger.info(
        "Importing codon translation table and amino acid names and annotating as globals..."
    )
    context_ht = context_ht.annotate_globals(
        codon_translation=get_codon_lookup(), acid_names=get_acid_names()
    )

    # NOTE: trimers must be True for this join to work correctly (for ExAC mu ht)
    logger.info("Importing mutation rates and joining to context HT...")
    mu_ht = mutation_rate.ht()
    context_ht = context_ht.annotate(
        mu_snp=mu_ht[context_ht.context, context_ht.alleles].mu_snp
    )

    logger.info("Importing divergence scores and annotating context HT...")
    div_ht = divergence_scores.ht()
    context_ht = context_ht.annotate(
        div=hl.if_else(
            hl.is_defined(
                div_ht[context_ht.vep.transcript_consequences[0].transcript_id]
            ),
            div_ht[context_ht.vep.transcript_consequences[0].transcript_id].divergence,
            0.0564635,
        )
    )

    logger.info("Moving transcript ID and exon number to top level annotation...")
    context_ht = context_ht.annotate(
        transcript=context_ht.vep.transcript_consequences[0].transcript_id,
        exon=context_ht.vep.transcript_consequences[0].exon.split("\/")[0],
    )

    logger.info("Writing out context HT...")
    context_ht = context_ht.naive_coalesce(n_partitions)
    context_ht.write(output_path, overwrite=overwrite)
    context_ht.describe()


## Functions to process exon/transcript related resourcces
def process_gencode_ht(build: str) -> None:
    """
    Imports gencode gtf as ht,filters to protein coding transcripts, and writes out as ht

    :param str build: Reference genome build; must be one of BUILDS.
    :return: None
    :rtype: None
    """
    if build not in BUILDS:
        raise DataException(f"Build must be one of {BUILDS}.")

    logger.info("Reading in gencode gtf")
    if build == "GRCh37":
        ht = grch37.gencode.ht()

    else:
        ht = grch38.gencode.ht()

    logger.info("Filtering gencode gtf to exons in protein coding genes...")
    ht = ht.filter((ht.feature == "exon") & (ht.gene_type == "protein_coding"))

    logger.info("Keying by transcript and exon number...")
    # Stripping decimal from transcript to match transcript in exome data
    ht = ht.transmute(transcript=ht.transcript_id.split("\.")[0])
    return ht.key_by("transcript", "exon_number").select()


## Functions for obs/exp related resources
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
        # TODO: check about impose_high_af_cutoff upfront
        keep_criteria = (ht.ac > 0) & (ht.pass_filters)

    return keep_criteria


def filter_to_missense(ht: hl.Table, n_partitions: int = 5000) -> hl.Table:
    """
    Filters input Table to missense variants in canonical transcripts only.

    :param Table ht: Input Table to be filtered.
    :param int n_partitions: Number of desired partitions for output.
    :return: Table filtered to only missense variants.
    :rtype: hl.Table
    """
    logger.info("Filtering to canonical transcripts...")
    ht = filter_vep_to_canonical_transcripts(ht)

    logger.info("Annotating HT with most severe consequence...")
    ht = add_most_severe_csq_to_tc_within_ht(ht)
    logger.info(
        f"Consequence count: {ht.aggregate(hl.agg.counter(ht.vep.most_severe_consequence))}"
    )

    logger.info("Filtering to missense variants...")
    ht = ht.filter(hl.literal(MISSENSE).contains(ht.vep.most_severe_consequence))

    return ht.naive_coalesce(n_partitions)


def filter_alt_decoy(ht: hl.Table) -> hl.Table:
    """
    Filters input Table to autosomes, X/X PAR, and Y (not mito or alt/decoy contigs). 

    Also annotates each locus with region type.

    :param Table ht: Input Table to be filtered/annotated.
    :return: Table filtered to autosomes/PAR and annotated with PAR status.
    :rtype: hl.Table
    """
    logger.info(f"HT count before filtration: {ht.count()}")
    ht = ht.annotate(
        region_type=hl.case()
        .when((ht.locus.in_autosome() | ht.locus.in_x_par()), "autosome_xpar")
        .when(ht.locus.in_x_nonpar(), "x_nonpar")
        .when(ht.locus.in_y_nonpar(), "y")
        .or_missing()
    )
    logger.info("Filtering to autosomes + X/Y..")
    ht = ht.filter(hl.is_defined(ht.region_type))
    logger.info(f"HT count after filtration: {ht.count()}")
    return ht
