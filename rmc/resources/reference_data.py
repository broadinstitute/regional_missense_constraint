"""Script containing reference resources."""
from gnomad.resources.resource_utils import (
    TableResource,
    VersionedTableResource,
)

from rmc.resources.basics import RESOURCE_BUILD_PREFIX
from rmc.resources.resource_utils import CURRENT_BUILD


####################################################################################
## Reference genome related resources
####################################################################################
filtered_context = VersionedTableResource(
    default_version="v1",
    versions={
        "v1": TableResource(
            path=f"{RESOURCE_BUILD_PREFIX}/reference_data/ht/context_fasta_snps_only_vep_v1.ht",
        )
    },
)
"""
Context Table filtered to missense variants in canonical protein coding transcripts and annotated with
probability of mutation for each variant, CpG status, gnomAD exome coverage, and methylation level.

Used to calculate the cumulative observed and expected missense values per locus and
generate regional missense constraint results.

NOTE: This resource is created with `process_vep`, which now filters to non-outlier transcripts by default.
However, for v2, this resource contains *all* canonical transcripts (including outliers).
"""

gene_model = TableResource(path=f"{RESOURCE_BUILD_PREFIX}/browser/b37_transcripts.ht")
"""
Table containing transcript start and stop positions displayed in the browser.
"""


####################################################################################
## Assessment related resources
####################################################################################
clinvar_path_mis = TableResource(
    path=f"{RESOURCE_BUILD_PREFIX}/reference_data/ht/clinvar_pathogenic_missense.ht",
)
"""
ClinVar pathogenic/likely pathogenic missense variants in haploinsufficient genes that cause severe disease.
"""

de_novo_tsv = f"{RESOURCE_BUILD_PREFIX}/reference_data/fordist_KES_combined_asc_dd_dnms_2020_04_21_annotated.txt"
"""
De novo missense variants from 31,058 cases with developmental disorders, 6,430 cases with autism spectrum disorders, and 2,179 controls.
Controls are the siblings of the autism cases.
Samples are from:
Kaplanis et al. (Evidence for 28 genetic disorders discovered by combining healthcare and research data.)
Satterstrom et al. (Large-Scale Exome Sequencing Study Implicates Both Developmental and Functional Changes in the Neurobiology of Autism.)
"""

de_novo = TableResource(
    path=f"{RESOURCE_BUILD_PREFIX}/reference_data/ht/dd_de_novo.ht",
)
"""
De novo missense variants from 37,488 cases and 2,179 controls (same as above).
"""

####################################################################################
## MPC related resources
####################################################################################
cadd = TableResource(
    path=f"gs://seqr-reference-data/{CURRENT_BUILD}/CADD/CADD_snvs_and_indels.v1.6.ht"
)
"""
Table with CADD (v1.6) raw and phredd scores.
"""
