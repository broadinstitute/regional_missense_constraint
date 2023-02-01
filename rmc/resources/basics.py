"""Script containing generic resources and Google cloud bucket paths."""
from gnomad.resources.resource_utils import TableResource

from rmc.resources.resource_utils import CURRENT_BUILD


######################################################################
## Google bucket resources
######################################################################
RMC_PREFIX = "gs://regional_missense_constraint"
"""
Path to bucket attached to regional missense constraint (RMC) project.

Contains all RMC and MPC score related data.

MPC stands for Missense badness, Polyphen-2, and Constraint.
"""

RESOURCE_PREFIX = f"{RMC_PREFIX}/resources"
"""
Path to any non-gnomAD or VEP context resource files required for RMC or MPC.
"""

RESOURCE_BUILD_PREFIX = f"{RESOURCE_PREFIX}/{CURRENT_BUILD}"
"""
Path to bucket for genome build-specific resource files required for RMC or MPC.
"""

LOGGING_PATH = f"{RMC_PREFIX}/logs"
"""
Path to bucket for hail logs.
"""

TEMP_PATH = f"{RMC_PREFIX}/temp"
"""
Path to bucket for temporary files.

Used when checkpointing intermediate files that can't be deleted immediately.
"""

TEMP_PATH_WITH_FAST_DEL = "gs://gnomad-tmp-4day/rmc"
"""
Path to bucket for temporary files.

Used when checkpointing intermediate files that can be removed immediately
(e.g., files generated by hail).
Used for RMC-relevant temporary files, including files for missense badness and MPC.
"""

TEMP_PATH_WITH_SLOW_DEL = "gs://gnomad-tmp/rmc"
"""
Path to bucket for temporary files.

Used when checkpointing intermediate files to be retained for a short period
and then can be removed (e.g., files to be accessed in a downstream step).
Used for RMC-relevant temporary files, including files for missense badness and MPC.
"""

MODEL_PREFIX = f"{RMC_PREFIX}/model"
"""
Path to bucket containing resources related to building the mutational models.

Data in this bucket is used to set up regional missense constraint (RMC) calculations.
"""

SINGLE_BREAK_TEMP_PATH = f"{TEMP_PATH}/single_breaks"
"""
Path to bucket to store temporary results for single breaks searches.
"""

SIMUL_BREAK_TEMP_PATH = f"{TEMP_PATH}/simul_breaks"
"""
Path to bucket to store temporary results for simultaneous searches.
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
AMINO_ACIDS_PREFIX = f"{RESOURCE_PREFIX}/amino_acids"
"""
Path to any amino-acid related resource files used to build MPC.
"""

CODON_TABLE_PATH = f"{AMINO_ACIDS_PREFIX}/codons_lookup.tsv"
"""
TSV file containing two columns: codon and three letter amino acid code.
"""

ACID_NAMES_PATH = f"{AMINO_ACIDS_PREFIX}/acid_names.tsv"
"""
TSV file containing amino acid names.

File has three columns: amino acid name, amino acid 3 letter code,
and amino acid 1 letter code.
"""
