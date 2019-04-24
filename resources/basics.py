import argparse
import logging
from gnomad_hail.resources import *
from gnomad_hail.utils.generic import file_exists,get_reference_genome
from constraint_utils.constraint_basics import *
#from gnomad_lof.constraint_utils.generic import * 


RESOURCE_PREFIX = 'gs://regional_missense_constraint/resources'
MODEL_PREFIX = f'{RESOURCE_PREFIX}/model/'
FLAGSHIP_LOF = 'gs://gnomad-public/papers/2019-flagship-lof/v1.0/'
BUILDS = [37, 38]


def get_codon_table_path() -> str:
    """
    Returns path to codon lookup table (codon > amino acid)

    :return: Path to codon table
    :rtype: str
    """
    return f'{RESOURCE_PREFIX}/codons_lookup.tsv'


def get_acid_names_path() -> str:
    """
    Returns path to amino acid table (full name, 3 letter name, 1 letter name)

    :return: Path to amino acid names table
    :rtype: str
    """
    return f'{RESOURCE_PREFIX}/acid_names.tsv'


def get_reference_path(build) -> str:
    """
    Returns path to reference fasta files 

    :param str build: Reference genome build; one of BUILDS
    :return: Path to reference fasta + fasta index
    :rtype: str
    """
    if build not in BUILDS:
        raise DataException(f'Build must be one of {BUILDS}.')
        
    if build == 37:
        return 'gs://hail-common/references/human_g1k_v37.fasta'
    else:
        return 'gs://hail-common/references/Homo_sapiens_assembly38.fasta' 


def get_raw_context_ht_path() -> str:
    """
    Returns path to raw context Table (reference fasta in Table form)

    :return: Path to raw context Table
    :rtype: str
    """
    return f'{RESOURCE_PREFIX}/raw_fasta.snps_only.unsplit.ht'


def get_vep_context_ht_path() -> str:
    """
    Returns path to VEP'd raw context Table

    :return: Path to raw context Table with VEP annotations
    :rtype: str
    """
    return f'{RESOURCE_PREFIX}/raw_fasta.snps_only.unsplit.vep.ht'


def get_context_ht_path() -> str:
    """
    Returns path to VEP'd raw context Table with multiallelics split

    :return: Path to raw context Table with VEP annotations with multiallelics split
    :rtype: str
    """
    return f'{RESOURCE_PREFIX}/raw_fasta.snps_only.vep.ht'


def get_processed_exomes_ht_path() -> str:
    """
    Returns path to gnomAD exomes Table annotated with context information

    :return: Path to processed gnoMAD exomes Table
    :rtype: str
    """
    return f'{MODEL_PREFIX}/exomes_processed.ht'


def get_processed_genomes_ht_path() -> str:
    """
    Returns path to VEP'd raw context Table with multiallelics split

    :return: Path to raw context Table with VEP annotations with multiallelics split
    :rtype: str
    """
    return f'{MODEL_PREFIX}/genomes_processed.ht'


def pre_process_data(build, overwrite) -> None:
    """
    Preprocesses data for constraint; code stolen from Konrad

    :param str build: Reference genome build; one of BUILDS
    :param bool overwrite: Whether to overwrite files
    :return: None; updates resource paths with processed data
    :rtype: None
    """
    raw_context_txt_path = f'{get_reference_path(build)}.gz'

    import_fasta(raw_context_txt_path, get_raw_context_ht_path, overwrite)
    vep_context_ht(raw_context_ht_path, get_vep_context_ht_path, overwrite)
    split_context_mt(get_vep_context_ht_path, {'exomes': coverage_ht_path('exomes'), 'genomes': coverage_ht_path('genomes')},
                    methylation_sites_mt_path(), get_context_ht_path, overwrite)
    pre_process_data(get_gnomad_public_data('genomes'), get_context_ht_path, get_processed_genomes_ht_path, overwrite)
    pre_process_data(get_gnomad_public_data('exomes'), get_context_ht_path, get_processed_exomes_ht_path, overwrite)


class DataException(Exception):
    pass
