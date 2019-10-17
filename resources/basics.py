import argparse
import hail as hl
import logging
from gnomad_hail.resources import *
from gnomad_hail.utils import *
from gnomad_hail.utils.generic import *
from constraint_utils.constraint_basics import *
#from gnomad_lof.constraint_utils.generic import * 


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("regional_missense_constraint_basics")
logger.setLevel(logging.INFO)


RESOURCE_PREFIX = 'gs://regional_missense_constraint/resources'
BUILDS = ['GRCh37', 'GRCh38']

# original regional missense constraint resource files
codon_table_path = f'{RESOURCE_PREFIX}/codons_lookup.tsv'
acid_names_path = f'{RESOURCE_PREFIX}/acid_names.tsv'
mutation_rate_table_path = f'{RESOURCE_PREFIX}/mutation_rate_table.tsv'
divergence_scores_path = f'{RESOURCE_PREFIX}/divsites_gencodev19_all_transcripts.tsv'

# constraint resource files
FLAGSHIP_LOF = 'gs://gnomad-public/papers/2019-flagship-lof/v1.0/'
MODEL_PREFIX = f'{FLAGSHIP_LOF}/model/'
processed_exomes_ht_path = f'{MODEL_PREFIX}/exomes_processed.ht'
filtered_exomes_ht_path = f'{RESOURCE_PREFIX}/ht/exomes_missense_only.ht'
processed_genomes_ht_path = f'{MODEL_PREFIX}/genomes_processed.ht'

# processed constraint resource files
mutation_rate_ht = f'{RESOURCE_PREFIX}/ht/mutation_rate.ht'


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
            line = line.strip().split(' ')
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
            line = line.strip().split('\t')
            acid_map[line[1]] = (line[0], line[2])
    return acid_map


def get_mutation_rate() -> dict:
    """
    Reads in mutation rate table and stores as ht

    :return: Dictionary of mutation rate information (key: context, value: (alt, mu_snp))
    :rtype: dict
    """
    mu = {}
    # from    n_kmer  p_any_snp_given_kmer    mu_kmer to      count_snp       p_snp_given_kmer        mu_snp
    with hl.hadoop_open(mutation_rate_table_path) as m:
        for line in m:
            context, n_kmer, p_any_snp_given_kmer, mu_kmer, new_kmer, count_snp, p_snp_given_kmer, mu_snp = line.strip().split('\t')
            mu[context] = (new_kmer[1], mu_snp)
    return mu


def get_mutation_rate_ht() -> hl.Table:
    """
    Reads in mutation rate table and stores as ht

    :return: Mutation rate information in ht
    :rtype: Table 
    """
    ht = hl.import_table(f'{mutation_rate_table_path}', impute=True)
    ht = ht.transmute(context=ht['from'], ref=ht['from'][1],
                                 alt=ht.to[1])
    ht = ht.key_by('context', 'ref', 'alt')
    ht = ht.select(ht.mu_snp)
    ht = ht.checkpoint(mutation_rate_ht, overwrite=True)
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
            transcript, score = line.strip().split('\t')
            try:
                div_scores[transcript.split('.')[0]] = float(score)
            except:
                continue
    return div_scores


def get_reference_path(build: str) -> str:
    """
    Returns path to reference fasta files 

    :param str build: Reference genome build; one of BUILDS
    :return: Path to reference fasta + fasta index
    :rtype: str
    """
    if build not in BUILDS:
        raise DataException(f'Build must be one of {BUILDS}.')
        
    if build == 'GRCh37':
        return 'gs://hail-common/references/human_g1k_v37.fasta'
    else:
        return 'gs://hail-common/references/Homo_sapiens_assembly38.fasta' 


def get_full_context_ht_path(build: str) -> str:
    """
    Returns path to reference fasta in ht form, filtered to SNPs, and annotated with VEP

    :param str build: Reference genome build; one of BUILDS
    :return: Path to full SNP reference fasta
    :rtype: str
    """
    # TODO: Add support for b38 
    if build not in BUILDS:
        raise DataException(f'Build must be one of {BUILDS}.')
        
    if build == 'GRCh37':
        return f'{FLAGSHIP_LOF}/context/Homo_sapiens_assembly19.fasta.snps_only.vep_20181129.ht'
    else:
        raise DataException('Sorry, no reference ht for b38 yet')


def get_processed_context_ht_path(build: str) -> str:
    """
    Returns path to reference fasta in ht form, filtered to SNPs, and annotated with VEP

    :param str build: Reference genome build; one of BUILDS
    :return: Path to full SNP reference fasta
    :rtype: str
    """
    # TODO: Add support for b38 
    if build not in BUILDS:
        raise DataException(f'Build must be one of {BUILDS}.')
        
    if build == 'GRCh37':
        return f'{RESOURCE_PREFIX}/ht/context/context_fasta_snps_only_vep_20190430.ht'
    else:
        raise DataException('Sorry, no reference ht for b38 yet')


def get_gencode_gtf_path(build: str) -> str:
    """
    Gets path to gencode gtf

    :param str build: Reference genome build; one of BUILDS
    :return: Full path to gencode gtf
    :rtype: str
    """
    if build not in BUILDS:
        raise DataException(f'Build must be one of {BUILDS}.')

    if build == 'GRCh37':
        return f'{RESOURCE_PREFIX}/gencode.v30lift37.basic.annotation.gtf'
    else:
        return f'{RESOURCE_PREFIX}/gencode.v30.basic.annotation.gtf'


def get_gencode_ht_path(build: str) -> str:
    """
    Gets path to gencode ht

    :param str build: Reference genome build; one of BUILDS
    :return: Full path to gencode ht
    :rtype: str
    """
    if build not in BUILDS:
        raise DataException(f'Build must be one of {BUILDS}.')

    if build == 'GRCh37':
        return f'{RESOURCE_PREFIX}/ht/context/gencode.v30lift37.basic.annotation.ht'
    else:
        return f'{RESOURCE_PREFIX}/ht/context/gencode.v30.basic.annotation.ht'


def get_processed_gencode_ht_path(build: str) -> str:
    """
    Gets path to gencode ht

    :param str build: Reference genome build; one of BUILDS
    :return: Full path to gencode ht
    :rtype: str
    """
    if build not in BUILDS:
        raise DataException(f'Build must be one of {BUILDS}.')

    if build == 'GRCh37':
        return f'{RESOURCE_PREFIX}/ht/context/gencode.v30lift37.exons.ht'
    else:
        return f'{RESOURCE_PREFIX}/ht/context/gencode.v30.basic.exons.ht'


def process_gencode_ht(build: str) -> None:
    """
    Imports gencode gtf as ht,filters to protein coding transcripts, and writes out as ht

    :param str build: Reference genome build; one of BUILDS
    :return: None
    :rtype: None
    """
    if build not in BUILDS:
        raise DataException(f'Build must be one of {BUILDS}.')

    logger.info('Reading in gencode gtf')
    ht = hl.experimental.import_gtf(get_gencode_gtf_path(build), reference_genome=f'GRCh{build}',
                                    skip_invalid_contigs=True)

    logger.info('Filtering gencode gtf to exons in protein coding genes with support levels 1 and 2')
    ht = ht.filter((ht.feature == 'exon') & (ht.gene_type == 'protein_coding') & (ht.level != '3'))
    ht = ht.checkpoint(get_gencode_ht_path(build), overwrite=True)

    logger.info('Grouping gencode ht by transcript ID')
    ht = ht.key_by(ht.transcript_id)
    ht = ht.select(ht.frame, ht.exon_number, ht.gene_name, ht.interval)
    ht = ht.collect_by_key()
    ht.write(get_processed_gencode_ht_path(build), overwrite=True)


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
        raise DataException(f'Build must be one of {BUILDS}.')

    
    logger.info('Reading in SNPs-only, VEP-annotated context ht')
    full_context_ht = prepare_ht(hl.read_table(get_full_context_ht_path(build)), trimers) # from constraint_basics

    logger.info(
                'Importing codon translation table, amino acid names,'
                'mutation rate, divergence scores and annotating as globals')
    full_context_ht = full_context_ht.annotate_globals(
                                                    codon_translation=get_codon_lookup(), acid_names=get_acid_names(),
                                                    mutation_rate=get_mutation_rate(), div_scores=get_divergence_scores())

    logger.info('Filtering to canonical protein coding transcripts')
    full_context_ht = full_context_ht.explode(full_context_ht.vep.transcript_consequences)
    context_ht = full_context_ht.filter(
                                        (full_context_ht.vep.transcript_consequences.biotype == 'protein_coding')
                                        & (full_context_ht.vep.transcript_consequences.canonical == 1))

    #logger.info('Importing mutation rate table and annotating as global')
    #mu_ht = get_mutation_rate_ht()

    context_ht.vep.transcript_consequences.transcript_id.show()

    logger.info('Re-keying context ht')
    context_ht = context_ht.key_by(context_ht.locus, context_ht.context, context_ht.ref, context_ht.alt)

    logger.info('Writing out context ht')
    context_ht.write(get_processed_context_ht_path(build), overwrite=True)
    context_ht.describe()


def filter_to_missense(ht: hl.Table) -> hl.Table:
    """
    Filters input table to missense variants only

    :param Table ht: Input ht to be filtered
    :return: Table filtered to only missense variants
    :rtype: hl.Table
    """
    logger.info(f'ht count before filtration: {ht.count()}') # this printed 17209972
    logger.info('Annotating ht with most severe consequence')
    ht = add_most_severe_csq_to_tc_within_ht(ht) # from constraint_basics

    # vep consequences from https://github.com/macarthur-lab/gnomad_hail/blob/master/utils/constants.py
    # missense definition from seqr searches
    missense = ['stop_lost', 'initiator_codon_variant', 'start_lost', 'protein_altering_variant', 'missense_variant']
    logger.info('Filtering to missense variants')
    ht = ht.filter(hl.literal(missense).contains(ht.vep.most_severe_consequence))
    logger.info(f'ht count after filtration: {ht.count()}') # this printed 6818793

    ht = ht.naive_coalesce(5000)
    return ht


def filter_to_autosome_and_par(ht: hl.Table) -> hl.Table:
    """
    Filters input table to autosomes and PAR regions. Also annotates each locus with whether it's in a PAR

    :param Table ht: Input ht to be filtered/annotated
    :return: Table filtered to autosomes/PAR and annotated with PAR status
    :rtype: hl.Table
    """
    logger.info(f'ht count before filtration: {ht.count()}') 
    logger.info('Filtering to autosomes and PAR')

    ht = ht.filter(ht.locus.in_autosome_or_par())
    ht = ht.annotate(in_par=(ht.locus_in_x_par() | ht.locus_in_y_par()))
    logger.info(f'ht count after filtration: {ht.count()}')
    return ht

    
class DataException(Exception):
    pass
