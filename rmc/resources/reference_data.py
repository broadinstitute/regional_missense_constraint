"""Script containing reference resources."""
from gnomad.resources.resource_utils import (
    TableResource,
    VersionedTableResource,
)

from rmc.resources.basics import (
    AMINO_ACIDS_PREFIX,
    RESOURCE_BUILD_PREFIX,
)
from rmc.resources.resource_utils import CURRENT_BUILD


REF_DATA_PREFIX = f"{RESOURCE_BUILD_PREFIX}/reference_data"
"""
Path to bucket containing reference data resources.
"""


####################################################################################
## Reference genome related resources
####################################################################################
filtered_context = VersionedTableResource(
    default_version="v1",
    versions={
        "v1": TableResource(
            path=f"{REF_DATA_PREFIX}/ht/context_fasta_snps_only_vep_v1.ht",
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
dosage_tsv = f"{REF_DATA_PREFIX}/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz"
"""
Path to TSV of genes and dosage sensitivity scores.

Scores are pHaplo (haplosensitivity score) and pTriplo
(triplosensitivity score).

TSV was downloaded from: https://zenodo.org/record/6347673.

Data are from Collins et al. A cross-disorder dosage sensitivity map of the human genome.
(2022)
"""

dosage_ht = TableResource(path=f"{REF_DATA_PREFIX}/ht/dosage_sensitivity.ht")
"""
HT of genes and genes and dosage sensitivity scores.

Imported from `dosage_tsv`.
"""

haplo_genes_he = f"{REF_DATA_PREFIX}/ht/phaplo_genes.he"
"""
HailExpression of haploinsufficient genes.

List of HI genes was determined by filtering to genes with pHaplo >= 0.86.
"""

triplo_genes_he = f"{REF_DATA_PREFIX}/ht/ptriplo_genes.he"
"""
HailExpression of triplosensitive genes.

List of triplosensitive genes was determined by filtering to genes with pTriplo >= 0.94.
"""


####################################################################################
## Assessment related resources
####################################################################################
clinvar = TableResource(
    path="gs://seqr-reference-data/GRCh37/clinvar/clinvar.GRCh37.ht",
)
"""

Table of ClinVar variants maintained by the seqr team.

Last version of this HT accessed by RMC team corresponds to 20230121 ClinVar release.
"""

clinvar_plp_mis_haplo = TableResource(
    path=f"{REF_DATA_PREFIX}/ht/clinvar_pathogenic_missense_haplo.ht",
)
"""
ClinVar pathogenic/likely pathogenic missense variants in haploinsufficient (HI)
genes.
"""

clinvar_plp_mis_triplo = TableResource(
    path=f"{REF_DATA_PREFIX}/ht/clinvar_pathogenic_missense_triplo.ht",
)
"""
ClinVar pathogenic/likely pathogenic missense variants in triplosensitive genes.
"""

ddd_autism_de_novo_tsv = (
    f"{REF_DATA_PREFIX}/fordist_KES_combined_asc_dd_dnms_2020_04_21_annotated.txt"
)
"""
Path to de novo variants from 39,667 samples.

Sample count breakdown:
- 31,058 individuals with developmental delay or intellectual disability (DD/ID)
- 6,430 autistic individuals
- 2,179 control individuals without NDDs

Controls are the siblings of the autism cases.
Samples are from:
Kaplanis et al. Evidence for 28 genetic disorders discovered by
combining healthcare and research data. (2020)
Satterstrom et al. Large-Scale Exome Sequencing Study Implicates Both Developmental
and Functional Changes in the Neurobiology of Autism. (2020)
"""

asc_ssc_spark_de_novo_tsv = f"{REF_DATA_PREFIX}/fu_2022_supp20.txt"
"""
De novo variants from 20,528 samples.

Sample count* breakdown:
- 15,036 autistic probands
- 5,492 control siblings without NDDs
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
Fu et al. Rare coding variation provides insight into the genetic architecture
and phenotypic context of autism (2022)
"""

de_novo = TableResource(
    path=f"{RESOURCE_BUILD_PREFIX}/reference_data/ht/ddd_autism_de_novo.ht",
)
"""
De novo missense variants from 46,094 neurodevelopmental disorder (NDD) cases and 5,492 controls.

Cases:
- 31,058 individuals with DD/ID (from `ndd_de_novo_tsv`)
- 15,036 autistic individuals (from `autism_de_novo_2022_tsv`)

Controls:
- 5,492 siblings without NDDs (from `autism_de_novo_2022_tsv`)
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
