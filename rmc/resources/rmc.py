"""
Script containing RMC and MPC related resources.

RMC: Regional missense constraint
MPC: Missense badness, Polyphen-2, and Constraint score
"""
from typing import Set

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
    default_version=1,
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

CONSTRAINT_ANNOTATIONS = {
    "mu_snp",
    "observed",
    "coverage",
    "total_exp",
    "total_mu",
    "total_obs",
    "cumulative_obs",
    "cumulative_exp",
    "forward_oe",
    "mu_scan",
    "section_mu",
    "section_exp",
    "section_obs",
    "section_oe",
    "reverse",
    "reverse_obs_exp",
    "total_null",
    "total_alt",
    "chisq",
    "max_chisq",
}
"""
Set of annotations used to calculate constraint and to hold resulting statistics.

Used to drop unnecessary fields when starting rescue break search.

TODO: assess which annotations in this list can be removed
"""

FINAL_ANNOTATIONS = {
    "section_obs",
    "section_exp",
    "section_oe",
}
"""
Set of annotations to keep from individual break search round result HTs when finalizing release HT.
"""


def single_search_bucket_path(
    is_rescue: bool,
    search_num: int = None,
) -> str:
    """
    Return path to bucket associated with single break search inputs and results.

    Function returns path to top level initial or "rescue"
    (search with lowered chi square significance cutoff) bucket,
    or bucket based on search number in either initial or rescue bucket.

    :param is_rescue: Whether to return path corresponding to rescue pathway.
    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
        Default is None.
    :return: Path to single break search round bucket.
    """
    rescue = "rescue" if is_rescue else "initial"
    return (
        f"{SINGLE_BREAK_TEMP_PATH}/{rescue}/round{search_num}"
        if search_num
        else f"{SINGLE_BREAK_TEMP_PATH}/{rescue}"
    )


def single_search_round_ht_path(
    is_rescue: bool,
    search_num: int,
    is_break_found: bool,
    is_breakpoint_only: bool,
) -> str:
    """
    Return path to a Table with results from a specified round of single break search.

    Function returns path to HT based on search number, break status,
    breakpoint status, and whether HT is associated with "rescue" pathway
    (pathway with lowered chi square significance cutoff).

    Break status refers to whether transcripts/sections in HT have at least one
    single significant breakpoint.

    Breakpoint status refers to whether HT contains breakpoint positions only
    or all positions in transcripts/sections.

    :param is_rescue: Whether to return path to HT created in rescue pathway.
    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param is_break_found: Whether to return path to HT with transcript/sections
        that have significant single break results.
    :param is_breakpoint_only: Whether to return path to HT with breakpoint positions
        only.
    :return: Path to specified HT resulting from single break search.
    """
    rescue = "rescue" if is_rescue else "initial"
    break_status = "break_found" if is_break_found else "no_break_found"
    breakpoint_status = "_breakpoint_only" if is_breakpoint_only else ""
    return f"{SINGLE_BREAK_TEMP_PATH}/{rescue}/round{search_num}/{break_status}{breakpoint_status}.ht"


SIMUL_SEARCH_BUCKET_NAMES = {"prep", "raw_results", "final_results", "success_files"}
"""
Names of buckets nested within round bucket of `SIMUL_BREAK_TEMP_PATH`.

Bucket structure:
    `SIMUL_BREAK_TEMP_PATH`
        initial/ or rescue/
            round/
            (anything not specific to round number at this level)
                prep/
                raw_results/
                final_results/
                success_files/
"""

SIMUL_SEARCH_ANNOTATIONS = {"max_chisq", "breakpoints"}
"""
Set of annotations to keep from two simultaneous breaks search.

Used when merging sections found in over and under length threshold search.

`max_chisq`: Chi square value associated with two breaks.
`breakpoints`: Tuple of breakpoints with adjusted inclusiveness/exclusiveness.

Note that this field will also be kept (`section` is a key field):
`section`: Transcript section that was searched.
    Format: <transcript>_<start position>_<end position>.
"""


def simul_search_bucket_path(
    is_rescue: bool,
    search_num: int = None,
) -> str:
    """
    Return path to bucket associated with simultaneous break search inputs and results.

    Function returns path to top level initial or "rescue"
    (search with lowered chi square significance cutoff) bucket,
    or bucket based on search number in either initial or rescue bucket.

    :param is_rescue: Whether to return path corresponding to rescue pathway.
    :param search_num: Search iteration number
        (e.g., second round of searching for simultaneous break would be 2).
        Default is None.
    :return: Path to simultaneous break search round bucket.
    """
    rescue = "rescue" if is_rescue else "initial"
    return (
        f"{SIMUL_BREAK_TEMP_PATH}/{rescue}/round{search_num}"
        if search_num
        else f"{SIMUL_BREAK_TEMP_PATH}/{rescue}"
    )


def simul_search_round_bucket_path(
    is_rescue: bool,
    search_num: int,
    bucket_type: str,
    bucket_names: Set[str] = SIMUL_SEARCH_BUCKET_NAMES,
) -> str:
    """
    Return path to bucket with  Tables resulting from a specific round of simultaneous break search.

    Function returns path to bucket based on search number, whether search is in
    "rescue" pathway (pathway with lowered chi square significance cutoff), and
    bucket type.

    :param is_rescue: Whether to return path corresponding to rescue pathway.
    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param bucket_type: Bucket type.
        Must be in `bucket_names`.
    :param bucket_names: Possible bucket names for simultaneous search bucket type.
        Default is `SIMUL_SEARCH_BUCKET_NAMES`.
    :return: Path to a bucket in the simultaneous break search round bucket.
    """
    assert bucket_type in bucket_names, f"Bucket type must be one of {bucket_names}!"
    return f"{simul_search_bucket_path(is_rescue, search_num)}/{bucket_type}"


def grouped_single_no_break_ht_path(
    is_rescue: bool,
    search_num: int,
) -> str:
    """
    Return path to Table of transcripts/transcript sections without a significant break in a single break search round, grouped by transcript/transcript section.

    Function returns path to Table based on search number and whether search is
    in "rescue" pathway (pathway with lowered chi square significance cutoff).

    :param is_rescue: Whether to return path corresponding to rescue pathway.
    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :return: Path to grouped Table.
    """
    bucket_path = simul_search_round_bucket_path(
        is_rescue=is_rescue,
        search_num=search_num,
        bucket_type="prep",
    )
    return f"{bucket_path}/grouped_single_no_break_found.ht"


def simul_sections_split_by_len_path(
    is_rescue: bool,
    search_num: int,
    is_over_threshold: bool,
) -> str:
    """
    Return path to transcripts/transcript sections entering a specific round of simultaneous break search.

    Function returns path to SetExpression based on search number, whether search is
    in "rescue" pathway (pathway with lowered chi square significance cutoff), and
    whether the transcripts/transcript sections have greater than or equal to the
    cutoff for possible missense positions.

    :param is_rescue: Whether to return path corresponding to rescue pathway.
    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param is_over_threshold: Whether to return path for transcripts/transcript
        sections with more than the cutoff for possible missense positions specified
        in `prepare_transcripts.py`. If True, those with greater than or equal to
        this cutoff will be returned. If False, those with fewer than this cutoff
        will be returned.
    :return: Path to SetExpression containing transcripts/transcript sections.
    """
    bucket_path = simul_search_round_bucket_path(
        is_rescue=is_rescue,
        search_num=search_num,
        bucket_type="prep",
    )
    threshold_relation = "over" if is_over_threshold else "under"
    return f"{bucket_path}/sections_to_simul_{threshold_relation}_threshold.he"


def merged_search_ht_path(
    is_rescue: bool,
    search_num: int,
    is_break_found: bool,
) -> str:
    """
    Return path to Table with merged single and simultaneous breaks search results.

    Function returns path to HT based on search number, break status,
    and whether HT is associated with "rescue" pathway
    (pathway with lowered chi square significance cutoff).

    Break status refers to whether transcripts/sections in HT have at least one
    single significant breakpoint.

    :param is_rescue: Whether to return path to HT created in rescue pathway.
    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param is_break_found: Whether to return path to HT with transcript/sections
        that have significant single break results.
    :return: Path to merged break found HT or no break found HT.
    """
    rescue = "rescue_" if is_rescue else ""
    break_status = "break_found" if is_break_found else "no_break_found"
    return f"{TEMP_PATH}/{rescue}round{search_num}_merged_{break_status}.ht"


no_breaks = (
    f"{CONSTRAINT_PREFIX}/{CURRENT_GNOMAD_VERSION}/{CURRENT_FREEZE}/no_breaks.he"
)
"""
SetExpression containing transcripts with no significant breaks.
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
