"""Script containing reference resources."""
from gnomad.resources.resource_utils import (
    TableResource,
    VersionedTableResource,
)

from rmc.resources.basics import (
    AMINO_ACIDS_PREFIX,
    RESOURCE_PREFIX,
    RESOURCE_BUILD_PREFIX,
)
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


######################################################################
## Gene/transcript resources
######################################################################
hi_genes = f"{RESOURCE_PREFIX}/HI_genes.rCNV.txt"
"""
Path to haploinsufficient genes that cause severe disease.

List is from Ryan Collins.
"""


####################################################################################
## Assessment related resources
####################################################################################
REF_DATA_PREFIX = f"{RESOURCE_BUILD_PREFIX}/reference_data"
"""
Path to bucket containing reference data resources.
"""

clinvar_path_mis = TableResource(
    path=f"{REF_DATA_PREFIX}/ht/clinvar_pathogenic_missense.ht",
)
"""
ClinVar pathogenic/likely pathogenic missense variants in haploinsufficient genes that cause severe disease.
"""

ddd_autism_de_novo_tsv = (
    f"{REF_DATA_PREFIX}/fordist_KES_combined_asc_dd_dnms_2020_04_21_annotated.txt"
)
"""
De novo variants from 31,058 cases with developmental disorders, 6,430 cases with autism spectrum disorders, and 2,179 controls.

Controls are the siblings of the autism cases.
Samples are from:
Kaplanis et al. (Evidence for 28 genetic disorders discovered by combining healthcare and research data.)
Satterstrom et al. (Large-Scale Exome Sequencing Study Implicates Both Developmental and Functional Changes in the Neurobiology of Autism.)
"""

autism_spark_de_novo_tsv = f"{REF_DATA_PREFIX}/fu_2022_supp20.txt"
"""
De novo variants from 20,528 samples.

Sample count* breakdown:
- 15,036 probands
- 5,492 siblings
(28,522 parents)

.. note::
    Fu et al. note that:
    "one family is in both SPARK and unpublished ASC data, with different probands;
    one mother in the unpublished ASC data is also a proband in a different trio
    in the same dataset"

Samples are from the following cohorts:
- Autism Sequencing Consortium (ASC)
- Simons Foundation Autism Research Initiative (SFARI) Simons Simplex Collection (SSC)
- Simons Foundation Powering Autism Research for Knowledge (SPARK initiative)

Samples are from:
Fu et al. (Rare coding variation provides insight into the genetic architecture and phenotypic context of autism)
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

blosum_txt_path = f"{AMINO_ACIDS_PREFIX}/blosum62.txt"
"""
Text file containing matrix of BLOSUM scores for each amino acid pair.
"""

blosum = TableResource(
    path=f"{AMINO_ACIDS_PREFIX}/ht/blosum.ht",
)
"""
Table containing BLOSUM scores for each amino acid pair.

Hail Table representation of scores in `blosum_txt_path`.
"""

grantham_txt_path = f"{AMINO_ACIDS_PREFIX}/grantham.matrix.txt"
"""
Text file containing matrix of Grantham scores for each amino acid pair.
"""

grantham = TableResource(
    path=f"{AMINO_ACIDS_PREFIX}/ht/grantham.ht",
)
"""
Table containing Grantham scores for each amino acid pair.

Hail Table representation of scores in `grantham_txt_path`.
"""
