from gnomad.resources.resource_utils import (
    TableResource,
    VersionedTableResource,
)
from rmc.resources.basics import MPC_PREFIX
from rmc.resources.resource_utils import (
    CURRENT_VERSION,
    RESOURCE_PREFIX,
)


## Reference genome related resources
processed_context = VersionedTableResource(
    default_version="v1",
    versions={
        "v1": TableResource(
            path=f"{RESOURCE_PREFIX}/GRCh37/reference_data/ht/context_fasta_snps_only_vep_v1.ht",
        )
    },
)
"""
Context Table filtered to missense variants in canonical protein coding transcripts and annotated with
probability of mutation for each variant, CpG status, gnomAD exome coverage, and methylation level.

Used to calculate the cumulative observed and expected missense values per locus and
generate regional missense constraint results.
"""

gene_model = TableResource(path=f"{RESOURCE_PREFIX}/GRCh37/browser/b37_transcripts.ht")
"""
Table containing transcript start and stop positions displayed in the browser.
"""


## Assessment related resources
clinvar_path_mis = TableResource(
    path=f"{RESOURCE_PREFIX}/GRCh37/reference_data/ht/clinvar_pathogenic_missense.ht",
)
"""
ClinVar pathogenic/likely pathogenic missense variants in haploinsufficient genes that cause severe disease.
"""

de_novo_tsv = f"{RESOURCE_PREFIX}/GRCh37/reference_data/fordist_KES_combined_asc_dd_dnms_2020_04_21_annotated.txt"
"""
De novo missense variants from 31,058 cases with developmental disorders, 6,430 cases with autism spectrum disorders, and 2,179 controls.
Controls are the siblings of the autism cases.
Samples are from:
Kaplanis et al. (Evidence for 28 genetic disorders discovered by combining healthcare and research data.)
Satterstrom et al. (Large-Scale Exome Sequencing Study Implicates Both Developmental and Functional Changes in the Neurobiology of Autism.)
"""

de_novo = TableResource(
    path=f"{RESOURCE_PREFIX}/GRCh37/reference_data/ht/dd_de_novo.ht",
)
"""
De novo missense variants from 37,488 cases and 2,179 controls (same as above).
"""

DE_NOVO_COUNTS = {
    "dd_only": 31058,
    "asd_only": 6430,
    "controls": 2179,
}
"""
Dictionary with counts of neurodevelopmental disorders (NDD) cases and controls in `de_novo` resources.
"""

## MPC related resources
cadd = TableResource(
    path="gs://seqr-reference-data/GRCh37/CADD/CADD_snvs_and_indels.v1.6.ht"
)
"""
Table with CADD (v1.6) raw and phredd scores.
"""


def get_mpc_case_control_ht_path(asd_only: bool, dd_only: bool) -> str:
    """
    Return path to case control Table.

    Table is used to create stacked histogram of MPC scores in NDD cases vs controls
    and rate ratio Table.

    Function will return path to one of three possible Tables:
        - Table with developmental disorders (DD) cases + controls only (no Autism Spectrum Disorders cases)
        - Table with Autism Spectrum Disorders (ASD) cases + controls only (no DD cases)
        - Table with all NDD cases + controls

    :param bool asd_only: Whether to return path to Table with ASD cases and controls only (no DD cases).
    :param bool dd_only: Whether to return path to Table with developmental disorders (DD) cases and controls only (no ASD cases).
    :return: Path to desired case control Table
    """
    if asd_only:
        return f"{MPC_PREFIX}/{CURRENT_VERSION}/asd_mpc_hist.ht"
    elif dd_only:
        return f"{MPC_PREFIX}/{CURRENT_VERSION}/dd_mpc_hist.ht"
    return f"{MPC_PREFIX}/{CURRENT_VERSION}/all_ndd_mpc_hist.ht"
