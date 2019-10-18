import argparse
import hail as hl
import logging
import pickle
from gnomad_hail.resources import *
from gnomad_hail.utils import *
from gnomad_hail.utils.generic import *
from constraint_utils.constraint_basics import *


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("regional_missense_constraint_basics")
logger.setLevel(logging.INFO)


RESOURCE_PREFIX = 'gs://regional_missense_constraint/resources'
BUILDS = ['GRCh37', 'GRCh38']

# missense variant VEP annotations
MISSENSE = ['stop_lost', 'initiator_codon_variant', 'start_lost', 'protein_altering_variant', 'missense_variant']

# original regional missense constraint resource files
codon_table_path = f'{RESOURCE_PREFIX}/codons_lookup.tsv'
acid_names_path = f'{RESOURCE_PREFIX}/acid_names.tsv'
mutation_rate_table_path = f'{RESOURCE_PREFIX}/mutation_rate_table.tsv'
divergence_scores_path = f'{RESOURCE_PREFIX}/divsites_gencodev19_all_transcripts.tsv'

# LoF constraint resource files
FLAGSHIP_LOF = 'gs://gnomad-public/papers/2019-flagship-lof/v1.0/'
MODEL_PREFIX = f'{FLAGSHIP_LOF}/model/'
processed_exomes_ht_path = f'{MODEL_PREFIX}/exomes_processed.ht'
filtered_exomes_ht_path = f'{RESOURCE_PREFIX}/ht/exomes_missense_only.ht'
processed_genomes_ht_path = f'{MODEL_PREFIX}/genomes_processed.ht'

# processed constraint resource files
mutation_rate_ht = f'{RESOURCE_PREFIX}/ht/mutation_rate.ht'

# ExAC files (for direct comparison with Kaitlin's code)
EXAC_PREFIX = f'{RESOURCE_PREFIX}/ExAC'
exac_vcf = f'{EXAC_PREFIX}/ExAC.r0.3.sites.vep.vcf.gz'
exac_ht = f'{EXAC_PREFIX}/ExAC.r0.3.sites.vep.ht'
filtered_exac_ht = f'{EXAC_PREFIX}/ExAC.r0.3.missense_only.ht'
filt_exac_cov_ht = f'{EXAC_PREFIX}/ExAC.r0.3.missense_only_cov.ht' # NOTE: this is only annotated with coverage for chr22


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
    Reads in mutation rate table and stores as dict

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
    ht = ht.transmute(alleles=[ht.ref, ht.alt])
    ht = ht.key_by('context', 'alleles')
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


## Reference genome related resources
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


## Exon/transcript related resourcces
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


## obs/exp related resources
# expected variants resource files
MODEL_PREFIX = 'gs://regional_missense_constraint/model'
EXP_PREFIX = f'{MODEL_PREFIX}/exp/'
exp_var_pickle = f'{EXP_PREFIX}/pickle/expected_variants.pckl'
exac_exp_var_pickle = f'{EXP_PREFIX}/pickle/exac_expected_variants.pckl' 
cov_ht = coverage_ht_path('exomes') # https://github.com/macarthur-lab/gnomad_hail/blob/master/resources/basics.py#L366
exac_tsv_path = 'gs://gnomad-public/legacy/exacv1_downloads/release0.1/coverage'
exac_cov_path = f'{EXAC_PREFIX}/coverage' 


def exp_ht_path(ExAC: bool=False) -> str:
    """
    Returns path of context ht annotated with expected variant counts

    :param bool ExAC: Whether to return table calculated on ExAC
    :return: Path to expected variants ht
    :rtype: str
    """
    return f'{EXP_PREFIX}/ExAC_exp_var.ht' if ExAC else f'{EXP_PREFIX}/exp_var.ht'


def load_exp_var(ExAC: bool=False) -> Dict[hl.Struct, int]: # NOTE: I think this is used for the mu calculations; don't need this yet (as of 10/18)
    """
    Loads saved expected variant count from pickle

    :param bool ExAC: Whether to load expected variants from ExAC
    :return: Dictionary of variant counts
    :rtype: dict
    """
    fname = exac_exp_var_pickle if filtered else exp_var_pickle
    with hl.hadoop_open(fname, 'rb') as f:
        return pickle.load(f)

    
class DataException(Exception):
    pass
