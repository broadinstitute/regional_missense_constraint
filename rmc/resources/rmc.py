"""
Script containing RMC and MPC related resources.

RMC: Regional missense constraint
MPC: Missense badness, Polyphen-2, and Constraint score
"""
from typing import Set

import hail as hl
import scipy
from gnomad.resources.resource_utils import (
    DataException,
    TableResource,
    VersionedTableResource,
)

from rmc.resources.basics import (
    CONSTRAINT_PREFIX,
    MODEL_PREFIX,
    MPC_PREFIX,
    RESOURCE_PREFIX,
    SIMUL_BREAK_TEMP_PATH,
    SINGLE_BREAK_TEMP_PATH,
    TEMP_PATH,
)
from rmc.resources.reference_data import FOLD_K
from rmc.resources.resource_utils import CURRENT_GNOMAD_VERSION

FREEZES = [1, 2, 3, 4, 5, 6, 7]
"""
RMC/MPC data versions computed with current gnomAD version.
"""

CURRENT_FREEZE = 7
"""
Current RMC/MPC data version.
"""

P_VALUE = 0.001
"""
Default p-value significance threshold.

Used to determine whether chi square values determining RMC breakpoints
are significant.

Default is 0.001.

Hail reference:
https://hail.is/docs/0.2/functions/stats.html#hail.expr.functions.qchisqtail
Look-up table reference:
https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm
"""

MIN_CHISQ_THRESHOLD = scipy.stats.chi2.ppf(0.975, 2)
"""
Minimum chi square significance.

Used only in two simultaneous breaks search.
Any breakpoint combinations with a chi square value less than this threshold
will not be emitted for storage.
"""

MIN_EXP_MIS = 16.0
"""
Minimum number of expected missense variants within each RMC section.


Sections that have fewer than this number of expected missense variants
will not be computed (chi square will be annotated as a missing value).

Default is 16.

This number was chosen as it corresponds to the minimum number
of expected missense variants needed on each side of a
potential breakpoint to reach exome-wide significance,
in the simple scenario where both sides of a potential breakpoint
have an equal number of expected missense variants,
with one side having an O/E value of 0 and the other having an O/E value of 1.

Exome-wide significance here is 2.7e-6,
calculated based on testing 18,629 transcripts total.

Calculated using a power curve with code at `plotting/expected_vs_significance.R`.

The same minimum expected value is used in both the single breaks and
simultaneous breaks searches.
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
coverage_plateau_models_path = f"{MODEL_PREFIX}/{CURRENT_GNOMAD_VERSION}/{CURRENT_FREEZE}/coverage_plateau_models.he"
"""
Path to HailExpression containing struct for coverage correction model and plateau models.

Schema:
    struct {
        'plateau_models': struct {
            total: dict<bool, array<float64>>
        }
        'plateau_x_models': struct {
            total: dict<bool, array<float64>>
        }
        'plateau_y_models': struct {
            total: dict<bool, array<float64>>
        }
        'coverage_model': tuple (
            float64,
            float64
        )
    }

Used to compute expected variant counts.
"""


# TODO: Delete current version
filtered_context = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{MODEL_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/context_coding_snps_annot.ht"
        )
        for freeze in FREEZES
    },
)
"""
Variant-level VEP context Table filtered to missense, nonsense, read-through and synonymous
variants in all canonical protein-coding transcripts.

Table contains constraint-related annotations, including observed variant counts,
expected variant counts, probability of mutation, CpG status, gnomAD exome coverage,
and methylation level.

Schema:
----------------------------------------
Global fields:
    'plateau_models': struct {
        total: dict<bool, array<float64>>
    }
    'plateau_x_models': struct {
        total: dict<bool, array<float64>>
    }
    'plateau_y_models': struct {
        total: dict<bool, array<float64>>
    }
    'coverage_model': tuple (
        float64,
        float64
    )
----------------------------------------
Row fields:
    'locus': locus<GRCh37>
    'alleles': array<str>
    'context': str
    'ref': str
    'alt': str
    'methylation_level': int32
    'cpg': bool
    'mutation_type': str
    'annotation': str
    'modifier': str
    'coverage': int32
    'transcript': str
    'expected': float64
    'coverage_correction': float64
    'observed': int32
----------------------------------------
Key: ['locus', 'alleles']
----------------------------------------

Used to create the constraint prep Table.
"""


constraint_prep = VersionedTableResource(
    default_version=CURRENT_FREEZE,
    versions={
        freeze: TableResource(
            path=f"{MODEL_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/constraint_prep.ht"
        )
        for freeze in FREEZES
    },
)
"""
Locus-level Table used in first step of regional constraint calculation.

Filtered to only one specific coding variant consequence but contains all canonical, protein-coding transcripts.

Schema:
----------------------------------------
Global fields:
    'plateau_models': struct {
        total: dict<bool, array<float64>>
    }
    'plateau_x_models': struct {
        total: dict<bool, array<float64>>
    }
    'plateau_y_models': struct {
        total: dict<bool, array<float64>>
    }
    'coverage_model': tuple (
        float64,
        float64
    )
----------------------------------------
Row fields:
    'locus': locus<GRCh37>
    'section': str
    'observed': int32
    'expected': float64
----------------------------------------
Key: ['locus', 'section']
----------------------------------------

NOTE: The content of the freeze 1 HT for RMC is slightly different:
    - Expected counts were calculated using the workaround method instead of directly from the plateau model.
    - `mu_snp` contains the coverage-corrected not raw mutation rates.
"""

FINAL_ANNOTATIONS = {
    "section_obs",
    "section_exp",
    "section_oe",
}
"""
Set of annotations to keep from individual break search round result HTs when finalizing release HT.
"""


# NOTE: Removed all references to rescue search pathway,
# but freeze 2 and 3 temp results were written with code that
# differentiated between "initial" and "rescue" search
def single_search_bucket_path(
    search_num: int = None,
    freeze: int = CURRENT_FREEZE,
) -> str:
    """
    Return path to bucket associated with single break search inputs and results.

    Function returns path to top level bucket or bucket based on search number.

    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
        Default is None.
    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: Path to single break search round bucket.
    """
    return (
        f"{SINGLE_BREAK_TEMP_PATH}/{freeze}/round{search_num}"
        if search_num
        else f"{SINGLE_BREAK_TEMP_PATH}/{freeze}"
    )


def single_search_round_ht_path(
    search_num: int,
    is_break_found: bool,
    is_breakpoint_only: bool,
    freeze: int = CURRENT_FREEZE,
) -> str:
    """
    Return path to a Table with results from a specified round of single break search.

    Function returns path to HT based on search number, break status, and breakpoint status.

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
    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: Path to specified HT resulting from single break search.
    """
    break_status = "break_found" if is_break_found else "no_break_found"
    breakpoint_status = "_breakpoint_only" if is_breakpoint_only else ""
    return (
        f"{single_search_bucket_path(search_num, freeze)}/{break_status}{breakpoint_status}.ht"
    )


SIMUL_SEARCH_BUCKET_NAMES = {"prep", "raw_results", "final_results", "success_files"}
"""
Names of buckets nested within round bucket of `SIMUL_BREAK_TEMP_PATH`.

Bucket structure:
    `SIMUL_BREAK_TEMP_PATH`
        freeze/
            round/
            (anything not specific to round number at this level)
                prep/
                raw_results/
                final_results/
                success_files/
"""

SIMUL_SEARCH_ANNOTATIONS = {"chisq", "breakpoints"}
"""
Set of annotations to keep from two simultaneous breaks search.

Used when merging sections found in over and under length threshold search.

`chisq`: Chi square value associated with two breaks.
`breakpoints`: Tuple of breakpoints with adjusted inclusiveness/exclusiveness.

Note that this field will also be kept (`section` is a key field):
`section`: Transcript section that was searched.
    Format: <transcript>_<start position>_<end position>.
"""


def simul_search_bucket_path(
    search_num: int = None,
    freeze: int = CURRENT_FREEZE,
) -> str:
    """
    Return path to bucket associated with simultaneous break search inputs and results.

    Function returns path to top level bucket or bucket based on search number.

    :param search_num: Search iteration number
        (e.g., second round of searching for simultaneous break would be 2).
        Default is None.
    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: Path to simultaneous break search round bucket.
    """
    return (
        f"{SIMUL_BREAK_TEMP_PATH}/{freeze}/round{search_num}"
        if search_num
        else f"{SIMUL_BREAK_TEMP_PATH}/{freeze}"
    )


def simul_search_round_bucket_path(
    search_num: int,
    bucket_type: str,
    bucket_names: Set[str] = SIMUL_SEARCH_BUCKET_NAMES,
    freeze: int = CURRENT_FREEZE,
) -> str:
    """
    Return path to bucket with  Tables resulting from a specific round of simultaneous break search.

    Function returns path to bucket based on search number and bucket type.

    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param bucket_type: Bucket type.
        Must be in `bucket_names`.
    :param bucket_names: Possible bucket names for simultaneous search bucket type.
        Default is `SIMUL_SEARCH_BUCKET_NAMES`.
    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: Path to a bucket in the simultaneous break search round bucket.
    """
    assert bucket_type in bucket_names, f"Bucket type must be one of {bucket_names}!"
    return f"{simul_search_bucket_path(search_num, freeze)}/{bucket_type}"


def grouped_single_no_break_ht_path(
    search_num: int,
    freeze: int = CURRENT_FREEZE,
) -> str:
    """
    Return path to Table of transcripts/transcript sections without a significant break in a single break search round, grouped by transcript/transcript section.

    Function returns path to Table based on search number.

    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: Path to grouped Table.
    """
    bucket_path = simul_search_round_bucket_path(
        search_num=search_num,
        bucket_type="prep",
        freeze=freeze,
    )
    return f"{bucket_path}/grouped_single_no_break_found.ht"


def simul_sections_split_by_len_path(
    search_num: int,
    is_over_threshold: bool,
    freeze: int = CURRENT_FREEZE,
) -> str:
    """
    Return path to transcripts/transcript sections entering a specific round of simultaneous break search.

    Function returns path to SetExpression based on search number and
    whether the transcripts/transcript sections have greater than or equal to the
    cutoff for possible missense positions.

    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param is_over_threshold: Whether to return path for transcripts/transcript
        sections with more than the cutoff for possible missense positions specified
        in `prepare_transcripts.py`. If True, those with greater than or equal to
        this cutoff will be returned. If False, those with fewer than this cutoff
        will be returned.
    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: Path to SetExpression containing transcripts/transcript sections.
    """
    bucket_path = simul_search_round_bucket_path(
        search_num=search_num,
        bucket_type="prep",
        freeze=freeze,
    )
    threshold_relation = "over" if is_over_threshold else "under"
    return f"{bucket_path}/sections_to_simul_{threshold_relation}_threshold.he"


def merged_search_ht_path(
    search_num: int,
    is_break_found: bool = True,
    freeze: int = CURRENT_FREEZE,
) -> str:
    """
    Return path to Table with merged single and simultaneous breaks search results.

    Function returns path to HT for break found sections based on search number.

    Function also has ability to return path to HailExpression containing
    set of sections without breakpoints, though this functionality isn't
    currently being used.

    Break status refers to whether transcripts/sections in HT have at least one
    single significant breakpoint.

    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param is_break_found: Whether to return path to HT with transcript/sections
        that have significant single break results.
        Default is True.
    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: Path to merged break found HT or no break found HailExpression.
    """
    if is_break_found:
        return f"{TEMP_PATH}/freeze{freeze}/round{search_num}_merged_break_found.ht"
    return f"{TEMP_PATH}/freeze{freeze}/round{search_num}_no_break_found.he"


def no_breaks_he_path(
    freeze: int = CURRENT_FREEZE,
) -> str:
    """
    Return path to HailExpression containing transcripts with no evidence of RMC.

    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: Path to no break transcripts HailExpression.
    """
    return f"{CONSTRAINT_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/no_breaks.he"


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
def amino_acids_oe_path(
    fold: int = None,
    freeze: int = CURRENT_FREEZE,
) -> str:
    """
    Return path to Table containing all possible amino acid substitutions and their missense OE ratio.

    Table is input to missense badness calculations.

    :param int fold: Fold number in training set to select training transcripts from.
        If not None, the Table is generated from variants in only training transcripts
        from the specified fold of the overall training set. If None, the Table is generated from
        variants in all training transcripts. Default is None.
    :param int freeze: RMC data freeze number. Default is CURRENT_FREEZE.
    :return: Path to Table.
    """
    if fold is not None and fold not in range(1, FOLD_K + 1):
        raise DataException(
            f"Fold number must be an integer between 1 and {FOLD_K}, inclusive!"
        )
    fold_name = f"_fold{fold}" if fold is not None else ""
    return f"{MPC_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/train{fold_name}/amino_acid_oe.ht"


def misbad_path(
    fold: int = None,
    freeze: int = CURRENT_FREEZE,
) -> str:
    """
    Table containing all possible amino acid substitutions and their missense badness scores.

    :param int fold: Fold number in training set to select training transcripts from.
        If not None, the Table is generated from variants in only training transcripts
        from the specified fold of the overall training set. If None, the Table is generated from
        variants in all training transcripts. Default is None.
    :param int freeze: RMC data freeze number. Default is CURRENT_FREEZE.
    :return: Path to Table.
    """
    if fold is not None and fold not in range(1, FOLD_K + 1):
        raise DataException(
            f"Fold number must be an integer between 1 and {FOLD_K}, inclusive!"
        )
    fold_name = f"_fold{fold}" if fold is not None else ""
    return f"{MPC_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/train{fold_name}/missense_badness.ht"


####################################################################################
## MPC related resources
####################################################################################
def joint_clinvar_gnomad_path(
    fold: int = None,
    freeze: int = CURRENT_FREEZE,
) -> str:
    """
    Return path to Table containing 'population' and 'pathogenic' variants.

    Table contains common (AF > 0.001) gnomAD variants ('population') and
    ClinVar pathogenic/likely pathogenic missense variants in haploinsufficient genes
    that cause severe disease ('pathogenic') with defined CADD, BLOSUM, Grantham, missense observed/expected ratios,
    missense badness, and PolyPhen-2 scores.

    Table is input to MPC (missense badness, polyphen-2, and constraint) calculations.

    :param int fold: Fold number in training set to select training transcripts from.
        If not None, the Table is generated from variants in only  training transcripts
        from the specified fold of the overall training set. If None, the Table is generated from
        variants in all training transcripts. Default is None.
    :param int freeze: RMC data freeze number. Default is CURRENT_FREEZE.
    :return: Path to Table.
    """
    if fold is not None and fold not in range(1, FOLD_K + 1):
        raise DataException(
            f"Fold number must be an integer between 1 and {FOLD_K}, inclusive!"
        )
    fold_name = f"_fold{fold}" if fold is not None else ""
    return f"{MPC_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/train{fold_name}/joint_clinvar_gnomad.ht"


def mpc_model_pkl_path(
    fold: int = None,
    freeze: int = CURRENT_FREEZE,
) -> str:
    """
    Return path to model (stored as pickle) that contains relationship of MPC variables.

    Model created using logistic regression.

    :param int fold: Fold number in training set to select training transcripts from.
        If not None, the model is generated from variants in only training transcripts
        from the specified fold of the overall training set. If None, the model is generated from
        variants in all training transcripts. Default is None.
    :param int freeze: RMC data freeze number. Default is CURRENT_FREEZE.
    :return: Path to model.
    """
    if fold is not None and fold not in range(1, FOLD_K + 1):
        raise DataException(
            f"Fold number must be an integer between 1 and {FOLD_K}, inclusive!"
        )
    fold_name = f"_fold{fold}" if fold is not None else ""
    return (
        f"{MPC_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/train{fold_name}/mpc_model.pkl"
    )


def gnomad_fitted_score_path(
    is_grouped: bool = False,
    fold: int = None,
    freeze: int = CURRENT_FREEZE,
) -> str:
    """
    Return path to fitted scores (from MPC model regression) of common (AF > 0.001) gnomAD variants.

    Table is input to MPC (missense badness, polyphen-2, and constraint) calculations on other datasets.

    :param bool is_grouped: Whether the Table is grouped by score. Default is False.
    :param int fold: Fold number in training set to select training transcripts from.
        If not None, the Table is generated from variants in only training transcripts
        from the specified fold of the overall training set. If None, the Table is generated from
        variants in all training transcripts. Default is None.
    :param int freeze: RMC data freeze number. Default is CURRENT_FREEZE.
    :return: Path to Table.
    """
    if fold is not None and fold not in range(1, FOLD_K + 1):
        raise DataException(
            f"Fold number must be an integer between 1 and {FOLD_K}, inclusive!"
        )
    fold_name = f"_fold{fold}" if fold is not None else ""
    group = "_group" if is_grouped else ""
    return f"{MPC_PREFIX}/{CURRENT_GNOMAD_VERSION}/{freeze}/train{fold_name}/gnomad_fitted_scores{group}.ht"


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

This Table contains only one row per each unique locus/alleles combination.
"""

####################################################################################
## Assessment related resources
####################################################################################
oe_bin_counts_tsv = (
    f"{CONSTRAINT_PREFIX}/{CURRENT_GNOMAD_VERSION}/{CURRENT_FREEZE}/oe_bin.tsv"
)
"""
TSV with RMC regions grouped by obs/exp (OE) bin.

Annotated with proportion coding base pairs, proportion de novo missense (controls),
proportion de novo missense (case), and proportion ClinVar pathogenic/likely pathogenic
severe haploinsufficient missense.
"""
