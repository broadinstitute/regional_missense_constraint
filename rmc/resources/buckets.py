"""Script containing Google cloud bucket paths."""


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

MPC stands for Missense badness, Polyphen-2, and Constraint.
"""
