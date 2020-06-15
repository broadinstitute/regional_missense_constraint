import logging
import pickle
from typing import Dict

import hail as hl

from gnomad.resources.resource_utils import DataException
from .resource_utils import BUILDS, RESOURCE_PREFIX


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_missense_constraint_basics")
logger.setLevel(logging.INFO)


## Kaitlin's resources
# Original regional missense constraint resource files
codon_table_path = f"{RESOURCE_PREFIX}/codons_lookup.tsv"
acid_names_path = f"{RESOURCE_PREFIX}/acid_names.tsv"
mutation_rate_table_path = f"{RESOURCE_PREFIX}/mutation_rate_table.tsv"
divergence_scores_path = f"{RESOURCE_PREFIX}/divsites_gencodev19_all_transcripts.tsv"
divergence_ht = f"{RESOURCE_PREFIX}/div_scores.ht"

## Konrad's resources
# LoF constraint resource files
FLAGSHIP_LOF = "gs://gnomad-public/papers/2019-flagship-lof/v1.0/"
MODEL_PREFIX = f"{FLAGSHIP_LOF}/model/"
processed_exomes_ht_path = f"{MODEL_PREFIX}/exomes_processed.ht"
filtered_exomes_ht_path = f"{RESOURCE_PREFIX}/ht/exomes_missense_only.ht"
processed_genomes_ht_path = f"{MODEL_PREFIX}/genomes_processed.ht"

# Processed constraint resource files
mutation_rate_ht = f"{RESOURCE_PREFIX}/ht/mutation_rate.ht"

## ExAC resources
# ExAC files (for direct comparison with Kaitlin's code)
EXAC_PREFIX = f"{RESOURCE_PREFIX}/ExAC"
exac_vcf = f"{EXAC_PREFIX}/ExAC.r0.3.sites.vep.vcf.gz"
exac_ht = f"{EXAC_PREFIX}/ExAC.r0.3.sites.vep.ht"
filtered_exac_ht = f"{EXAC_PREFIX}/ExAC.r0.3.missense_only.ht"
filt_exac_cov_ht = f"{EXAC_PREFIX}/ExAC.r0.3.missense_only_cov.ht"  # NOTE: this is only annotated with coverage for chr22


## Reference genome related resources
def get_reference_path(build: str) -> str:
    """
    Returns path to reference fasta files 

    :param str build: Reference genome build; one of BUILDS
    :return: Path to reference fasta + fasta index
    :rtype: str
    """
    if build not in BUILDS:
        raise DataException(f"Build must be one of {BUILDS}.")

    if build == "GRCh37":
        return "gs://hail-common/references/human_g1k_v37.fasta"
 
    return "gs://hail-common/references/Homo_sapiens_assembly38.fasta"


def get_full_context_ht_path(build: str) -> str:
    """
    Returns path to reference fasta in ht form, filtered to SNPs, and annotated with VEP

    :param str build: Reference genome build; one of BUILDS
    :return: Path to full SNP reference fasta
    :rtype: str
    """
    # TODO: Add support for b38
    if build not in BUILDS:
        raise DataException(f"Build must be one of {BUILDS}.")

    if build == "GRCh37":
        return f"{FLAGSHIP_LOF}/context/Homo_sapiens_assembly19.fasta.snps_only.vep_20181129.ht"

    raise DataException("Sorry, no reference ht for b38 yet")


def get_processed_context_ht_path(build: str) -> str:
    """
    Returns path to reference fasta in ht form, filtered to SNPs, and annotated with VEP

    :param str build: Reference genome build; one of BUILDS
    :return: Path to full SNP reference fasta
    :rtype: str
    """
    # TODO: Add support for b38
    if build not in BUILDS:
        raise DataException(f"Build must be one of {BUILDS}.")

    if build == "GRCh37":
        return f"{RESOURCE_PREFIX}/ht/context/context_fasta_snps_only_vep_20190430.ht"
    raise DataException("Sorry, no reference ht for b38 yet")


## Exon/transcript related resourcces
def get_gencode_gtf_path(build: str) -> str:
    """
    Gets path to gencode gtf

    :param str build: Reference genome build; one of BUILDS
    :return: Full path to gencode gtf
    :rtype: str
    """
    if build not in BUILDS:
        raise DataException(f"Build must be one of {BUILDS}.")

    if build == "GRCh37":
        return f"{RESOURCE_PREFIX}/gencode.v30lift37.basic.annotation.gtf"
    return f"{RESOURCE_PREFIX}/gencode.v30.basic.annotation.gtf"


def get_gencode_ht_path(build: str) -> str:
    """
    Gets path to gencode ht

    :param str build: Reference genome build; one of BUILDS
    :return: Full path to gencode ht
    :rtype: str
    """
    if build not in BUILDS:
        raise DataException(f"Build must be one of {BUILDS}.")

    if build == "GRCh37":
        return f"{RESOURCE_PREFIX}/ht/context/gencode.v30lift37.basic.annotation.ht"
    return f"{RESOURCE_PREFIX}/ht/context/gencode.v30.basic.annotation.ht"


def get_processed_gencode_ht_path(build: str) -> str:
    """
    Gets path to gencode ht

    :param str build: Reference genome build; one of BUILDS
    :return: Full path to gencode ht
    :rtype: str
    """
    if build not in BUILDS:
        raise DataException(f"Build must be one of {BUILDS}.")

    if build == "GRCh37":
        return f"{RESOURCE_PREFIX}/ht/context/gencode.v30lift37.exons.ht"
    return f"{RESOURCE_PREFIX}/ht/context/gencode.v30.basic.exons.ht"


## Observed/expected count related resources
# Expected variants resource files
MODEL_PREFIX = "gs://regional_missense_constraint/model"
EXP_PREFIX = f"{MODEL_PREFIX}/exp/"
exp_var_pickle = f"{EXP_PREFIX}/pickle/expected_variants.pckl"
exac_exp_var_pickle = f"{EXP_PREFIX}/pickle/exac_expected_variants.pckl"
exac_tsv_path = "gs://gnomad-public/legacy/exacv1_downloads/release0.1/coverage"
exac_cov_path = f"{EXAC_PREFIX}/coverage"


def exp_ht_path(exac: bool = False) -> str:
    """
    Returns path of context ht annotated with expected variant counts

    :param bool exac: Whether to return table calculated on ExAC
    :return: Path to expected variants ht
    :rtype: str
    """
    return f"{EXP_PREFIX}/ExAC_exp_var.ht" if exac else f"{EXP_PREFIX}/exp_var.ht"


def load_exp_var(
    exac: bool = False,
) -> Dict[
    hl.Struct, int
]:  # NOTE: I think this is used for the mu calculations; don't need this yet (as of 10/18)
    """
    Loads saved expected variant count from pickle

    :param bool ExAC: Whether to load expected variants from ExAC
    :return: Dictionary of variant counts
    :rtype: dict
    """
    fname = exac_exp_var_pickle if exac else exp_var_pickle
    with hl.hadoop_open(fname, "rb") as f:
        return pickle.load(f)
