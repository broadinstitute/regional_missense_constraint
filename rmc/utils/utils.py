import logging

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad_lof.constraint.constraint_basics import (
    add_most_severe_csq_to_tc_within_ht,
    prepare_ht,
)
from rmc.resources.basics import (
    acid_names_path,
    codon_table_path,
    divergence_ht,
    divergence_scores_path,
    get_full_context_ht_path,
    get_gencode_gtf_path,
    get_gencode_ht_path,
    get_processed_context_ht_path,
    get_processed_gencode_ht_path,
    mutation_rate_ht,
    mutation_rate_table_path,
)
from rmc.resources.resource_utils import BUILDS, MISSENSE


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_missense_constraint_generic")
logger.setLevel(logging.INFO)


## Resources from Kaitlin
def get_codon_lookup() -> dict:
    """
    Reads in codon lookup table and returns as dictionary (key: codon, value: amino acid)

    :return: Dictionary of codon translation
    :rtype: dict
    """
    codon_lookup = {}
    with hl.hadoop_open(codon_table_path) as c:
        c.readline()
        for line in c:
            line = line.strip().split(" ")
            codon_lookup[line[0]] = line[1]
    return codon_lookup


def get_acid_names() -> dict:
    """
    Reads in amino acid table and stores as dict (key: 3 letter name, value: (long name, one letter name)

    :return: Dictionary of amino acid names
    :rtype: dict
    """
    acid_map = {}
    with hl.hadoop_open(acid_names_path) as a:
        a.readline()
        for line in a:
            line = line.strip().split("\t")
            acid_map[line[1]] = (line[0], line[2])
    return acid_map


def get_mutation_rate() -> dict:
    """
    Reads in mutation rate table and stores as dict

    :return: Dictionary of mutation rate information (key: context, value: (alt, mu_snp))
    :rtype: dict
    """
    mu = {}
    # from    n_kmer  p_any_snp_given_kmer    mu_kmer to      count_snp       p_snp_given_kmer        mu_snp
    with hl.hadoop_open(mutation_rate_table_path) as m:
        for line in m:
            (
                context,
                n_kmer,
                p_any_snp_given_kmer,
                mu_kmer,
                new_kmer,
                count_snp,
                p_snp_given_kmer,
                mu_snp,
            ) = line.strip().split("\t")
            mu[context] = (new_kmer[1], mu_snp)
    return mu


def get_mutation_rate_ht(overwrite: bool = True) -> hl.Table:
    """
    Reads in mutation rate table and stores as ht

    :param bool overwrite: Whether to overwrite data
    :return: Mutation rate information in ht
    :rtype: Table
    """
    # from	n_kmer	p_any_snp_given_kmer	mu_kmer	to	count_snp	p_snp_given_kmer	mu_snp
    ht = hl.import_table(mutation_rate_table_path, impute=True)
    ht = ht.transmute(context=ht["from"], ref=ht["from"][1], alt=ht.to[1])
    ht = ht.transmute(alleles=[ht.ref, ht.alt])
    ht = ht.key_by("context", "alleles")
    ht = ht.select(ht.mu_snp)
    ht = ht.checkpoint(mutation_rate_ht, overwrite=overwrite)
    return ht


def get_divergence_scores() -> dict:
    """
    Reads in divergence score file and stores as dict (key: transcript, value: score)

    :return: Divergence score dict
    :rtype: dict
    """
    div_scores = {}
    with hl.hadoop_open(divergence_scores_path) as d:
        d.readline()
        for line in d:
            transcript, score = line.strip().split("\t")
            try:
                div_scores[transcript.split(".")[0]] = float(score)
            except:
                continue
    return div_scores


def get_divergence_score_ht(overwrite: bool = True) -> hl.Table:
    """
    Reads in divergence score file and writes out to ht

    :param bool overwrite: Whether to overwrite
    :return: Divergence score ht
    :rtype: Table
    """
    ht = hl.import_table(divergence_scores_path, impute=True)
    ht = ht.transmute(transcript=ht.transcript.split(".")[0])
    ht = ht.key_by("transcript")
    ht.write(divergence_ht, overwrite=overwrite)
    return ht


## Functions to process reference genome related resources
def process_context_ht(build: str, trimers: bool) -> None:
    """
    Imports reference fasta (SNPs only, VEP'd) as ht
    Filters to canonical protein coding transcripts

    :param str build: Reference genome build; one of BUILDS
    :param bool trimers: Whether to filter to trimers or heptamers
    :return: None
    :rtype: None
    """
    if build not in BUILDS:
        raise DataException(f"Build must be one of {BUILDS}.")

    logger.info("Reading in SNPs-only, VEP-annotated context ht")
    full_context_ht = prepare_ht(
        hl.read_table(get_full_context_ht_path(build)), trimers
    )  # from constraint_basics
    logger.info(f"Full ht count: {full_context_ht.count()}")

    logger.info("Filtering to canonical protein coding transcripts")
    full_context_ht = full_context_ht.explode(
        full_context_ht.vep.transcript_consequences
    )
    context_ht = full_context_ht.filter(
        (full_context_ht.vep.transcript_consequences.biotype == "protein_coding")
        & (full_context_ht.vep.transcript_consequences.canonical == 1)
    )
    logger.info(f"Count after filtration: {context_ht.count()}")

    logger.info(
        "Importing codon translation table and amino acid names and annotating as globals"
    )
    full_context_ht = full_context_ht.annotate_globals(
        codon_translation=get_codon_lookup(), acid_names=get_acid_names()
    )

    logger.info("Importing mutation rates and joining to context ht")
    mu_ht = get_mutation_rate_ht()
    context_ht = context_ht.key_by("context", "alleles").join(mu_ht, how="left")

    logger.info("Importing divergence scores and annotating context ht")
    div_ht = get_divergence_score_ht()
    context_ht = context_ht.annotate(
        div=hl.cond(
            hl.is_defined(div_ht[context_ht.vep.transcript_consequences.transcript_id]),
            div_ht[context_ht.vep.transcript_consequences.transcript_id].divergence,
            0.0564635,
        )
    )

    logger.info("Re-keying context ht by locus and alleles")
    # context_ht = context_ht.key_by(context_ht.locus, context_ht.context, context_ht.ref, context_ht.alt)
    context_ht = context_ht.key_by("locus", "alleles")

    logger.info("Writing out context ht")
    context_ht = context_ht.naive_coalesce(10000)
    context_ht.write(get_processed_context_ht_path(build), overwrite=True)
    context_ht.describe()


## Functions to process exon/transcript related resourcces
def process_gencode_ht(build: str) -> None:
    """
    Imports gencode gtf as ht,filters to protein coding transcripts, and writes out as ht

    :param str build: Reference genome build; one of BUILDS
    :return: None
    :rtype: None
    """
    if build not in BUILDS:
        raise DataException(f"Build must be one of {BUILDS}.")

    logger.info("Reading in gencode gtf")
    ht = hl.experimental.import_gtf(
        get_gencode_gtf_path(build),
        reference_genome=f"GRCh{build}",
        skip_invalid_contigs=True,
    )

    logger.info(
        "Filtering gencode gtf to exons in protein coding genes with support levels 1 and 2"
    )
    ht = ht.filter(
        (ht.feature == "exon") & (ht.gene_type == "protein_coding") & (ht.level != "3")
    )
    ht = ht.checkpoint(get_gencode_ht_path(build), overwrite=True)

    logger.info("Grouping gencode ht by transcript ID")
    ht = ht.key_by(ht.transcript_id)
    ht = ht.select(ht.frame, ht.exon_number, ht.gene_name, ht.interval)
    ht = ht.collect_by_key()
    ht.write(get_processed_gencode_ht_path(build), overwrite=True)


## Functions for obs/exp related resources
def filter_to_missense(ht: hl.Table) -> hl.Table:
    """
    Filters input table to missense variants only

    :param Table ht: Input ht to be filtered
    :return: Table filtered to only missense variants
    :rtype: hl.Table
    """
    logger.info(f"ht count before filtration: {ht.count()}")  # this printed 17209972
    logger.info("Annotating ht with most severe consequence")
    ht = add_most_severe_csq_to_tc_within_ht(ht)  # from constraint_basics
    logger.info(
        f"Consequence count: {ht.aggregate(hl.agg.counter(ht.vep.most_severe_consequence))}"
    )

    # vep consequences from https://github.com/macarthur-lab/gnomad_hail/blob/master/utils/constants.py
    # missense definition from seqr searches
    logger.info("Filtering to missense variants")
    ht = ht.filter(hl.literal(MISSENSE).contains(ht.vep.most_severe_consequence))
    logger.info(f"ht count after filtration: {ht.count()}")  # this printed 6818793

    ht = ht.naive_coalesce(5000)
    return ht


def filter_alt_decoy(ht: hl.Table) -> hl.Table:
    """
    Filters input table to autosomes, X/X PAR, and Y (not mito or alt/decoy contigs). 

    Also annotates each locus with region type.

    :param Table ht: Input ht to be filtered/annotated
    :return: Table filtered to autosomes/PAR and annotated with PAR status
    :rtype: hl.Table
    """
    logger.info(f"ht count before filtration: {ht.count()}")
    ht = ht.annotate(
        region_type=hl.case()
        .when((ht.locus.in_autosome() | ht.locus.in_x_par()), "autosome_xpar")
        .when(ht.locus.in_x_nonpar(), "x_nonpar")
        .when(ht.locus.in_y_nonpar(), "y")
        .default("remove")
    )
    logger.info("Filtering to autosomes + X/Y")
    ht = ht.filter(ht.region_type == "remove", keep=False)
    logger.info(f"ht count after filtration: {ht.count()}")
    return ht
