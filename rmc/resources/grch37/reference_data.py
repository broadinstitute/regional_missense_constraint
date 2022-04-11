from gnomad.resources.resource_utils import (
    TableResource,
    VersionedTableResource,
)

from rmc.resources.resource_utils import RESOURCE_PREFIX


## Reference genome related resources
VEP_VERSION = 85
"""
VEP version used to annotate full context Table.
"""

processed_context = VersionedTableResource(
    default_version="v1",
    versions={
        "v1": TableResource(
            path=f"{RESOURCE_PREFIX}/GRCh37/reference_data/ht/context_fasta_snps_only_vep_v1.ht",
        )
    },
)

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
