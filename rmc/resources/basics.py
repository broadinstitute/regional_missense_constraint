import hail as hl

from gnomad.resources.resource_utils import TableResource, VersionedTableResource
from rmc.resources.resource_utils import (
    FLAGSHIP_LOF,
    GNOMAD_VER,
    RESOURCE_PREFIX,
    RMC_PREFIX,
    UKBB_VER,
)


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

## gnomAD resources
mutation_rate = TableResource(
    path=f"{FLAGSHIP_LOF}/model/mutation_rate_methylation_bins.ht",
)
"""
Table with mutation rate recalculated for gnomAD constraint.

This was calculated with `calculate_mu_by_downsampling` in
https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/constraint_basics.py.
"""

## Observed/expected related resources
MODEL_PREFIX = f"{RMC_PREFIX}/model"

constraint_prep = VersionedTableResource(
    default_version=GNOMAD_VER,
    versions={
        GNOMAD_VER: TableResource(
            path=f"{MODEL_PREFIX}/{GNOMAD_VER}/context_obs_exp_annot.ht"
        ),
        UKBB_VER: TableResource(
            path=f"{MODEL_PREFIX}/ukbb_freeze{UKBB_VER}/context_obs_exp_annot.ht"
        ),
    },
)
"""
Context Table ready for RMC calculations.

HT is annotated with observed and expected variant counts per base.
"""


## Constraint related resources
CONSTRAINT_PREFIX = f"{RMC_PREFIX}/constraint"
one_break = VersionedTableResource(
    default_version=GNOMAD_VER,
    versions={
        GNOMAD_VER: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{GNOMAD_VER}/one_break.ht"
        ),
        UKBB_VER: TableResource(
            path=f"{CONSTRAINT_PREFIX}/ukbb_freeze{UKBB_VER}/one_break.ht"
        ),
    },
)
"""
Table containing transcripts with at least one break.

Found when searching constraint_prep HT for transcripts for a single (first) break.
"""

not_one_break = VersionedTableResource(
    default_version=GNOMAD_VER,
    versions={
        GNOMAD_VER: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{GNOMAD_VER}/not_one_break.ht"
        ),
        UKBB_VER: TableResource(
            path=f"{CONSTRAINT_PREFIX}/ukbb_{UKBB_VER}/not_one_break.ht"
        ),
    },
)
"""
Table containing transcripts without one significant break.

Transcripts in this table will be processed to check for two simultaneous breaks.
"""

not_one_break_grouped = VersionedTableResource(
    default_version=GNOMAD_VER,
    versions={
        GNOMAD_VER: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{GNOMAD_VER}/not_one_break_grouped.ht"
        ),
    },
)
"""
Not one break Table grouped by transcript with observed missense, expected missense, and positions collected into lists.

Input to searching for simultaneous breaks.
"""

multiple_breaks = VersionedTableResource(
    default_version=GNOMAD_VER,
    versions={
        GNOMAD_VER: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{GNOMAD_VER}/multiple_breaks.ht"
        ),
        UKBB_VER: TableResource(
            path=f"{CONSTRAINT_PREFIX}/ukbb_freeze{UKBB_VER}/multiple_breaks.ht"
        ),
    },
)
"""
Table containing transcripts with multiple breaks.
"""

simul_break_under_10k = f"{MODEL_PREFIX}/{GNOMAD_VER}/transcripts_under_10k.he"
"""
SetExpression containing transcripts with < 10,000 possible missense positions.
"""

simul_break_over_10k = f"{MODEL_PREFIX}/{GNOMAD_VER}/transcripts_over_10k.he"
"""
SetExpression containing transcripts with >= 10,000 possible missense positions.
"""

simul_break_temp = f"{temp_path}/simul_breaks/"
"""
Bucket to store temporary results for simultaneous results
"""

simul_break = VersionedTableResource(
    default_version=GNOMAD_VER,
    versions={
        GNOMAD_VER: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{GNOMAD_VER}/simul_break.ht"
        ),
        UKBB_VER: TableResource(
            path=f"{CONSTRAINT_PREFIX}/ukbb_freeze{UKBB_VER}/simul_break.ht"
        ),
    },
)
"""
Table containing transcripts with two simultaneous breaks.
"""

breaks = VersionedTableResource(
    default_version=GNOMAD_VER,
    versions={
        GNOMAD_VER: TableResource(path=f"{CONSTRAINT_PREFIX}/{GNOMAD_VER}/breaks.ht"),
        UKBB_VER: TableResource(
            path=f"{CONSTRAINT_PREFIX}/ukbb_freeze{UKBB_VER}/breaks.ht"
        ),
    },
)
"""
Table containing transcripts with any evidence of RMC (one break, multiple breaks, simultaneous breaks).
"""

no_breaks = VersionedTableResource(
    default_version=GNOMAD_VER,
    versions={
        GNOMAD_VER: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{GNOMAD_VER}/no_breaks.ht"
        ),
        UKBB_VER: TableResource(
            path=f"{CONSTRAINT_PREFIX}/ukbb_freeze{UKBB_VER}/no_breaks.ht"
        ),
    },
)
"""
Table containing transcripts with no significant breaks.
"""


TOTAL_EXOME_BASES = {"GRCh37": 54426835}
"""
Dictionary containing total number of bases in the exome.

Calculated using `get_exome_bases`.
"""


TOTAL_GNOMAD_MISSENSE = {"2.1.1": 5257859}
"""
Dictionary containing total number of missense variants seen in gnomAD.

Calculated by filtering gnomAD release HT to missense variants only and running `ht.count()`.
"""
