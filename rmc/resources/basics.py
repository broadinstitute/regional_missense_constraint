"""Script containing generic resources and Google cloud bucket paths."""


######################################################################
## Google bucket resources
######################################################################
RMC_PREFIX = "gs://regional_missense_constraint"
"""
Path to bucket attached to regional missense constraint (RMC) project.

Contains all RMC and MPC score related data.

MPC stands for Missense badness, Polyphen-2, and Constraint.
"""

RESOURCE_PREFIX = "gs://regional_missense_constraint/resources"
"""
Path to any non-gnomAD or VEP context resource files required for RMC or MPC.
"""

LOGGING_PATH = "gs://regional_missense_constraint/logs"
"""
Path to bucket for hail logs.
"""

TEMP_PATH = "gs://regional_missense_constraint/temp"
"""
Path to bucket for temporary files.

Used when checkpointing intermediate files that can't be deleted immediately.
"""

MODEL_PREFIX = f"{RMC_PREFIX}/model"
"""
Path to bucket containing resources related to building the mutational models.

Data in this bucket is used to set up regional missense constraint (RMC) calculations.
"""

SIMUL_BREAK_TEMP = f"{TEMP_PATH}/simul_breaks"
"""
Bucket to store temporary results for simultaneous results.
"""

CONSTRAINT_PREFIX = f"{RMC_PREFIX}/constraint"
"""
Path to bucket containing RMC output Tables.
"""

MPC_PREFIX = f"{RMC_PREFIX}/MPC"
"""
Path to bucket containing resources related to building MPC score.
"""


######################################################################
## Amino acid resources
######################################################################
CODON_TABLE_PATH = f"{RESOURCE_PREFIX}/amino_acids/codons_lookup.tsv"
"""
TSV file containing two columns: codon and three letter amino acid code.
"""

ACID_NAMES_PATH = f"{RESOURCE_PREFIX}/amino_acids/acid_names.tsv"
"""
TSV file containing amino acid names.

File has three columns: amino acid name, amino acid 3 letter code,
and amino acid 1 letter code.
"""

blosum_txt_path = f"{RESOURCE_PREFIX}/amino_acids/blosum62.txt"
"""
Text file containing matrix of BLOSUM scores for each amino acid pair.
"""

blosum = TableResource(
    path=f"{RESOURCE_PREFIX}/amino_acids/ht/blosum.ht",
)
"""
Table containing BLOSUM scores for each amino acid pair.

Hail Table representation of scores in `blosum_txt_path`.
"""

grantham_txt_path = f"{RESOURCE_PREFIX}/amino_acids/grantham.matrix.txt"
"""
Text file containing matrix of Grantham scores for each amino acid pair.
"""

grantham = TableResource(
    path=f"{RESOURCE_PREFIX}/amino_acids/ht/grantham.ht",
)
"""
Table containing Grantham scores for each amino acid pair.

Hail Table representation of scores in `grantham_txt_path`.
"""
