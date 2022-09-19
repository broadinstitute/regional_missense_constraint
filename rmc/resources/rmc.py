"""
Script containing RMC and MPC related resources.

RMC: Regional missense constraint
MPC: Missense badness, Polyphen-2, and Constraint score
"""
import hail as hl

from gnomad.resources.resource_utils import TableResource, VersionedTableResource

from rmc.resources.basics import (
    CONSTRAINT_PREFIX,
    MODEL_PREFIX,
    MPC_PREFIX,
    RESOURCE_PREFIX,
    SIMUL_BREAK_TEMP_PATH,
    SINGLE_BREAK_TEMP_PATH,
    TEMP_PATH,
)
from rmc.resources.resource_utils import CURRENT_GNOMAD_VERSION


FREEZES = [1, 2]
"""
RMC/MPC data versions computed with current gnomAD version.
"""

CURRENT_FREEZE = 2
"""
Current RMC/MPC data version.
"""

####################################################################################
## Original regional missense constraint resource files
####################################################################################
EXAC_PREFIX = f"{RESOURCE_PREFIX}/GRCh37/exac"
"""
Path to bucket containing ExAC constraint files.
"""

MUTATION_RATE_TABLE_PATH = f"{EXAC_PREFIX}/mutation_rate_table.tsv"
"""
Path to TSV containing ExAC mutation rates.
"""

DIVERGENCE_SCORES_TSV_PATH = f"{EXAC_PREFIX}/divsites_gencodev19_all_transcripts.txt"
"""
Path to text file with divergence scores per transcript.
"""

divergence_scores = TableResource(
    path=f"{EXAC_PREFIX}/ht/div_scores.ht",
    import_func=hl.import_table,
    import_args={
        "path": DIVERGENCE_SCORES_TSV_PATH,
        "key": "transcript",
        "min_partitions": 50,
        "impute": True,
    },
)
"""
Table with divergence score between humans and macaques
for each canonical transcript in Gencode v19.
"""

####################################################################################
## RMC-related resources
####################################################################################
constraint_prep = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{MODEL_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/context_obs_exp_annot.ht"
        )
        for freeze in FREEZES
    },
)
"""
Context Table ready for RMC calculations.

HT is annotated with observed and expected variant counts per base.
"""


def single_search_ht_path(
    search_num: int,
    is_break_found: bool,
    is_breakpoint_only: bool,
    is_rescue: bool,
) -> str:
    """
    Return path to a Table associated with results from a specified round of single
    break search.

    Function returns path to HT based on search number, break status,
    breakpoint status, and whether HT is associated with "rescue" pathway
    (pathway with lowered chi square significance cutoff).

    Break status refers to whether transcripts/sections in HT have at least one
    single significant breakpoint.

    Breakpoint status refers to whether HT contains breakpoint positions only
    or all positions in transcripts/sections.

    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param is_break_found: Whether to return path to HT with transcript/sections
        that have significant single break results.
    :param is_breakpoint_only: Whether to return path to HT with breakpoint positions
        only.
    :param is_rescue: Whether to return path to HT created in rescue pathway.
    :return: Path to specified single break found or no single break found HT.
    """
    rescue = "rescue_" if is_rescue else ""
    break_status = "break_found" if is_break_found else "no_break_found"
    breakpoint_status = "_breakpoint_only" if is_breakpoint_only else ""
    return f"{SINGLE_BREAK_TEMP_PATH}/{rescue}round{search_num}_single_{break_status}{breakpoint_status}.ht"


def grouped_no_single_break_ht_path(
    search_num: int,
    is_rescue: bool,
) -> str:
    """
    Return path to Table with results from single break search where no break was
    found, grouped by transcript/transcript section. This Table is used in
    preparation for simultaneous break search.

    Function returns path to HT based on search number, break status,
    breakpoint status, and whether HT is associated with "rescue" pathway
    (pathway with lowered chi square significance cutoff).

    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param is_rescue: Whether to return path to HT created in rescue pathway.
    :return: Path to specified single break found or no single break found HT.
    """
    rescue = "rescue_" if is_rescue else ""
    return f"{SIMUL_BREAK_TEMP_PATH}/{rescue}round{search_num}_grouped_single_no_break_found.ht"
    
    
def merged_search_ht_path(
    search_num: int,
    is_break_found: bool,
    is_rescue: bool,
) -> str:
    """
    Return path to Table with merged single and simultaneous breaks search results.

    Function returns path to HT based on search number, break status,
    and whether HT is associated with "rescue" pathway
    (pathway with lowered chi square significance cutoff).

    Break status refers to whether transcripts/sections in HT have at least one
    single significant breakpoint.

    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param is_break_found: Whether to return path to HT with transcript/sections
        that have significant single break results.
    :param is_rescue: Whether to return path to HT created in rescue pathway.
    :return: Path to merged break found HT or no break found HT.
    """
    rescue = "rescue_" if is_rescue else ""
    break_status = "break_found" if is_break_found else "no_break_found"
    return (
        f"{TEMP_PATH}/{rescue}round{search_num}_merged_{break_status}.ht"
    )


one_break = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/one_break.ht"
        )
        for freeze in FREEZES
    },
)
"""
Table containing transcripts with at least one break.

Found when searching constraint_prep HT for transcripts for a single (first) break.
"""

not_one_break = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/not_one_break.ht"
        )
        for freeze in FREEZES
    },
)
"""
Table containing transcripts without one significant break.

Transcripts in this table will be processed to check for two simultaneous breaks.
"""

not_one_break_grouped = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/not_one_break_grouped.ht"
        )
        for freeze in FREEZES
    },
)
"""
Not one break Table grouped by transcript with observed missense, expected missense, and positions collected into lists.

Input to searching for simultaneous breaks.
"""

multiple_breaks = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/multiple_breaks.ht"
        )
        for freeze in FREEZES
    },
)
"""
Table containing transcripts with multiple breaks.
"""

simul_break_under_threshold_path = f"{MODEL_PREFIX}/{CURRENT_GNOMAD_VERSION}/{CURRENT_FREEZE}/transcripts_under_threshold.he"
"""
SetExpression containing transcripts with fewer possible missense positions than cutoff specified in `run_simultaneous_breaks.py`.
"""

simul_break_over_threshold_path = f"{MODEL_PREFIX}/{CURRENT_GNOMAD_VERSION}/{CURRENT_FREEZE}/transcripts_over_threshold.he"
"""
SetExpression containing transcripts with greater than or equal to the cutoff for possible missense positions.
"""

simul_break = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/simul_break.ht"
        )
        for freeze in FREEZES
    },
)
"""
Table containing transcripts with two simultaneous breaks.
"""

no_breaks = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/no_breaks.ht"
        )
        for freeze in FREEZES
    },
)
"""
Table containing transcripts with no significant breaks.
"""

rmc_results = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/all_rmc.ht"
        )
        for freeze in FREEZES
    },
)
"""
Table containing all transcripts with evidence of regional missense constraint.

Contains transcripts with one or additional breaks plus simultaneous breaks results.
"""

context_with_oe = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/context_with_oe.ht"
        )
        for freeze in FREEZES
    },
)
"""
Table containing missense variants in canonical transcripts annotated with missense OE.

Each row is annotated with RMC transcript subsection missense OE if it exists, otherwise
the transcript level missense OE.
"""

context_with_oe_dedup = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/context_with_oe_dedup.ht"
        )
        for freeze in FREEZES
    },
)
"""
Deduplicated version of `context_with_oe`.

Some locus/alleles combinations are present more than once in `context_with_oe` if they have
a most severe consequence of 'missense_variant' in more than one canonical transcript.

This Table contains only one row per each unique locus/alleles combination.
"""

rmc_browser = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/rmc_browser.ht"
        )
        for freeze in FREEZES
    },
)
"""
Table containing all transcripts with evidence of regional missense constraint.

Contains same information as `rmc_results` but has different formatting for gnomAD browser.
"""

####################################################################################
## Missense badness related resources
####################################################################################
amino_acids_oe = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{MPC_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/amino_acid_oe.ht"
        )
        for freeze in FREEZES
    },
)
"""
Table containing all possible amino acid substitutions and their missense OE ratio.

Input to missense badness calculations.
"""

misbad = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{MPC_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/missense_badness.ht"
        )
        for freeze in FREEZES
    },
)
"""
Table containing all possible amino acid substitutions and their missense badness scores.
"""

####################################################################################
## MPC related resources
####################################################################################
joint_clinvar_gnomad = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{MPC_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/joint_clinvar_gnomad.ht"
        )
        for freeze in FREEZES
    },
)
"""
Table containing "population" and "pathogenic" variants.

Table contains common (AF > 0.001) gnomAD variants ("population") and
ClinVar pathogenic/likely pathogenic missense variants in haploinsufficient genes
that cause severe disease ("pathogenic") with defined CADD, BLOSUM, Grantham, missense observed/expected ratios,
missense badness, and PolyPhen-2 scores.

Input to MPC (missense badness, polyphen-2, and constraint) calculations.
"""

mpc_model_pkl_path = (
    f"{MPC_PREFIX}/{CURRENT_GNOMAD_VERSION}/{CURRENT_FREEZE}/mpc_model.pkl"
)
"""
Path to model (stored as pickle) that contains relationship of MPC variables.

Created using logistic regression.
"""

gnomad_fitted_score = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{MPC_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/gnomad_fitted_scores.ht"
        )
        for freeze in FREEZES
    },
)
"""
Table of gnomAD variants and their fitted scores (from MPC model regression).

Input to MPC (missense badness, polyphen-2, and constraint) calculations on other datasets.
"""

gnomad_fitted_score_group = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{MPC_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/gnomad_fitted_scores_group.ht"
        )
        for freeze in FREEZES
    },
)
"""
Table of fitted scores for common (AF > 0.001) variants in gnomAD, grouped by score.

Annotated with the total number of variants with and less than each score.
"""

mpc_release = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{MPC_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/mpc.ht"
        )
        for freeze in FREEZES
    },
)
"""
Table containing missense variants in canonical transcripts annotated with MPC.
"""

mpc_release_dedup = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{MPC_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/mpc_dedup.ht"
        )
        for freeze in FREEZES
    },
)
"""
Table containing missense variants in canonical transcripts annotated with MPC.

This Table contains only one row per each unique locus/alleles combination.
"""

####################################################################################
## Assessment related resources
####################################################################################
oe_bin_counts_tsv = f"{CONSTRAINT_PREFIX}/{CURRENT_FREEZE}/{CURRENT_FREEZE}/oe_bin.tsv"
"""
TSV with RMC regions grouped by obs/exp (OE) bin.

Annotated with proportion coding base pairs, proportion de novo missense (controls),
proportion de novo missense (case), and proportion ClinVar pathogenic/likely pathogenic
severe haploinsufficient missense.
"""
