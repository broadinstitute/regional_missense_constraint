"""Script containing reference resources."""
from gnomad.resources.resource_utils import DataException, TableResource

from rmc.resources.basics import AMINO_ACIDS_PREFIX, RESOURCE_BUILD_PREFIX
from rmc.resources.resource_utils import CURRENT_BUILD

REF_DATA_PREFIX = f"{RESOURCE_BUILD_PREFIX}/reference_data"
"""
Path to bucket containing reference data resources.
"""


####################################################################################
## Reference genome related resources
####################################################################################
gene_model = TableResource(path=f"{RESOURCE_BUILD_PREFIX}/browser/b37_transcripts.ht")
"""
Table containing transcript start and stop positions displayed in the browser.

Contains all transcripts displayed in the browser. (Contains more transcripts than `transcript_ref` below).
"""

transcript_ref = TableResource(
    path=f"{REF_DATA_PREFIX}/ht/canonical_transcripts_genes_coordinates.ht"
)
"""
Table containing canonical transcripts with key reference info:
gene name from gnomAD annotations, gene name from GENCODE v19,
Ensembl gene ID in GENCODE v19, overall and CDS start and end
coordinates from gnomAD annotations, and transcript strand.
"""

transcript_cds = TableResource(path=f"{REF_DATA_PREFIX}/ht/b37_cds_coords.ht")
"""
Table containing coordinates for coding parts of transcripts excluding introns and UTRs.
"""


######################################################################
## Gene/transcript resources
######################################################################
FOLD_K = 5
"""
Number of folds in the training set.
"""


def train_val_test_transcripts_path(
    is_test: bool = False,
    fold: int = None,
    is_val: bool = False,
) -> str:
    """
    Return path to HailExpression of transcripts used for model training or testing.

    By default, all training transcripts are returned.

    .. note::
        - `is_test` cannot be True if `fold` is defined.
        - `is_test` and `is_val` cannot be True simultaneously.
        - If `is_val` is True, `fold` must also be defined.
        - If specified, `fold` value must be between 1 and `FOLD_K`.

    :param is_test: Whether to return test transcripts.
        If False, training transcripts will be returned. If True, test transcripts will be returned.
        Default is False.
    :param fold: Fold number in training set to select transcripts from.
        If None, all training transcripts will be returned.
        If not None, only validation or training transcripts from the specified fold will be returned.
        Default is None.
    :param is_val: Whether to return validation transcripts.
        If False, training transcripts from the specified fold of the training set will be returned.
        If True, validation transcripts from the specified fold of the training set will be returned.
        Default is False.
    :return: Path to Table.
    """
    if is_test and is_val:
        raise DataException("Test and validation sets are mutually exclusive!")
    if is_test and fold is not None:
        raise DataException("Fold number cannot be specified for test set!")
    if is_val and fold is None:
        raise DataException("Fold number must be specified for validation transcripts!")
    if fold is not None and fold not in range(1, FOLD_K + 1):
        raise DataException(
            f"Fold number must be an integer between 1 and {FOLD_K}, inclusive!"
        )

    transcript_type = "test" if is_test else ("val" if is_val else "train")
    fold_name = f"_fold{fold}" if fold is not None else ""
    return f"{RESOURCE_BUILD_PREFIX}/{transcript_type}_transcripts{fold_name}.he"


dosage_tsv_path = (
    f"{REF_DATA_PREFIX}/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz"
)
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

Imported from TSV at `dosage_tsv_path`.
"""

haplo_genes_path = f"{REF_DATA_PREFIX}/ht/phaplo_genes.he"
"""
Path to HailExpression of haploinsufficient genes.

List of HI genes was determined by filtering to genes with pHaplo >= 0.86.
"""

triplo_genes_path = f"{REF_DATA_PREFIX}/ht/ptriplo_genes.he"
"""
Path to HailExpression of triplosensitive genes.

List of triplosensitive genes was determined by filtering to genes with pTriplo >= 0.94.
"""


####################################################################################
## Assessment related resources
####################################################################################
clinvar = TableResource(
    path=f"{REF_DATA_PREFIX}/ht/clinvar.GRCh37.ht",
)
"""

Table of ClinVar variants.

HT corresponds to 20230305 ClinVar release.
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
ClinVar pathogenic/likely pathogenic missense variants in triplosensitive (TS) genes.
"""

ndd_de_novo_2020_tsv_path = (
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

autism_de_novo_2022_tsv_path = f"{REF_DATA_PREFIX}/fu_2022_supp20.txt"
"""
Path to de novo variants from 20,528 samples.

Sample count* breakdown:
- 15,036 autistic probands
- 5,492 control siblings without NDDs

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

ndd_de_novo = TableResource(
    path=f"{REF_DATA_PREFIX}/ht/ndd_de_novo.ht",
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
