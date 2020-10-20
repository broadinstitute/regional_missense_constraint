import hail as hl

from gnomad.resources.resource_utils import TableResource
from gnomad_lof.constraint_utils.constraint_basics import get_old_mu_data
from rmc.resources.resource_utils import FLAGSHIP_LOF, RESOURCE_PREFIX


LOGGING_PATH = "gs://regional_missense_constraint/logs"
"""
Path to bucket that stores hail logs.
"""

temp_path = "gs://regional_missense_constraint/temp"
"""
Path to bucket to store temporary files. 
Used when checkpointing intermediate files.
"""


## Kaitlin's resources
# Original regional missense constraint resource files
CODON_TABLE_PATH = f"{RESOURCE_PREFIX}/amino_acids/codons_lookup.tsv"
ACID_NAMES_PATH = f"{RESOURCE_PREFIX}/amino_acids/acid_names.tsv"
MUTATION_RATE_TABLE_PATH = f"{RESOURCE_PREFIX}/GRCh37/exac/mutation_rate_table.tsv"
DIVERGENCE_SCORES_TSV_PATH = (
    f"{RESOURCE_PREFIX}/GRCh37/exac/divsites_gencodev19_all_transcripts.tsv"
)
divergence_scores = TableResource(
    path=f"{RESOURCE_PREFIX}/GRCh37/exac/ht/div_scores.ht",
    import_func=hl.import_table,
    import_args={
        "path": DIVERGENCE_SCORES_TSV_PATH,
        "key": "transcript",
        "min_partitions": 50,
        "impute": True,
    },
)
"""
Table with divergence score between humans and macaques for each canonical transcript in Gencode v19.
"""

old_mutation_rate = TableResource(
    path=f"{RESOURCE_PREFIX}/GRCh37/exac/ht/mutation_rate.ht",
    import_func=get_old_mu_data,
    import_args={
        "path": MUTATION_RATE_TABLE_PATH,
        "min_partitions": 50,
        "impute": True,
    },
)
"""
Table with mutation rate calculated for ExAC constraint.
"""

## gnomAD resources
mutation_rate = TableResource(
    path=f"{FLAGSHIP_LOF}/model/mutation_rate_methylation_bins.ht",
)
"""
Table with mutation rate recalculated for gnomAD constraint.

This was calculated with `calculate_mu_by_downsampling` in 
https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/constraint_basics.py.
"""

## Observed/expected count related resources
# Expected variants resource files
MODEL_PREFIX = "gs://regional_missense_constraint/model"
EXP_PREFIX = f"{MODEL_PREFIX}/exp/"
