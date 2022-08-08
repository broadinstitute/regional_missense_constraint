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
)
from rmc.resources.resource_utils import CURRENT_VERSION, GNOMAD_VERSIONS


# Original regional missense constraint resource files
MUTATION_RATE_TABLE_PATH = f"{RESOURCE_PREFIX}/GRCh37/exac/mutation_rate_table.tsv"
"""
Path to TSV containing ExAC mutation rates.
"""

DIVERGENCE_SCORES_TSV_PATH = (
    f"{RESOURCE_PREFIX}/GRCh37/exac/divsites_gencodev19_all_transcripts.txt"
)
"""
Path to text file with divergence scores per transcript.
"""

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


# RMC-related resources
constraint_prep = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(
            path=f"{MODEL_PREFIX}/{version}/context_obs_exp_annot.ht"
        )
        for version in GNOMAD_VERSIONS
    },
)
"""
Context Table ready for RMC calculations.

HT is annotated with observed and expected variant counts per base.
"""

one_break = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(path=f"{CONSTRAINT_PREFIX}/{version}/one_break.ht")
        for version in GNOMAD_VERSIONS
    },
)
"""
Table containing transcripts with at least one break.

Found when searching constraint_prep HT for transcripts for a single (first) break.
"""

not_one_break = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(path=f"{CONSTRAINT_PREFIX}/{version}/not_one_break.ht")
        for version in GNOMAD_VERSIONS
    },
)
"""
Table containing transcripts without one significant break.

Transcripts in this table will be processed to check for two simultaneous breaks.
"""

not_one_break_grouped = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        CURRENT_VERSION: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{CURRENT_VERSION}/not_one_break_grouped.ht"
        ),
    },
)
"""
Not one break Table grouped by transcript with observed missense, expected missense, and positions collected into lists.

Input to searching for simultaneous breaks.
"""

multiple_breaks = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(path=f"{CONSTRAINT_PREFIX}/{version}/multiple_breaks.ht")
        for version in GNOMAD_VERSIONS
    },
)
"""
Table containing transcripts with multiple breaks.
"""

simul_break_under_threshold = (
    f"{MODEL_PREFIX}/{CURRENT_VERSION}/transcripts_under_5k.he"
)
"""
SetExpression containing transcripts with fewer possible missense positions than cutoff specified in `run_simultaneous_breaks.py`.
"""

simul_break_over_threshold = f"{MODEL_PREFIX}/{CURRENT_VERSION}/transcripts_over_5k.he"
"""
SetExpression containing transcripts with greater than or equal to the cutoff for possible missense positions.
"""

simul_break = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(path=f"{CONSTRAINT_PREFIX}/{version}/simul_break.ht")
        for version in GNOMAD_VERSIONS
    },
)
"""
Table containing transcripts with two simultaneous breaks.
"""

no_breaks = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(path=f"{CONSTRAINT_PREFIX}/{version}/no_breaks.ht")
        for version in GNOMAD_VERSIONS
    },
)
"""
Table containing transcripts with no significant breaks.
"""

rmc_results = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(path=f"{CONSTRAINT_PREFIX}/{version}/all_rmc.ht")
        for version in GNOMAD_VERSIONS
    },
)
"""
Table containing all transcripts with evidence of regional missense constraint.

Contains transcripts with one or additional breaks plus simultaneous breaks results.
"""

context_with_oe = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(path=f"{CONSTRAINT_PREFIX}/{version}/context_with_oe.ht")
        for version in GNOMAD_VERSIONS
    },
)
"""
Table containing missense variants in canonical transcripts annotated with missense OE.

Each row is annotated with RMC transcript subsection missense OE if it exists, otherwise
the transcript level missense OE.
"""

context_with_oe_dedup = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{version}/context_with_oe_dedup.ht"
        )
        for version in GNOMAD_VERSIONS
    },
)
"""
Deduplicated version of `context_with_oe`.

Some locus/alleles combinations are present more than once in `context_with_oe` if they have
a most severe consequence of 'missense_variant' in more than one canonical transcript.

This Table contains only one row per each unique locus/alleles combination.
"""

rmc_browser = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(
            path=f"{CONSTRAINT_PREFIX}/{CURRENT_VERSION}/rmc_browser.ht"
        )
        for version in GNOMAD_VERSIONS
    },
)
"""
Table containing all transcripts with evidence of regional missense constraint.

Contains same information as `rmc_results` but has different formatting for gnomAD browser.
"""

# Missense badness related resources
amino_acids_oe = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(path=f"{MPC_PREFIX}/{version}/amino_acid_oe.ht")
        for version in GNOMAD_VERSIONS
    },
)
"""
Table containing all possible amino acid substitutions and their missense OE ratio.

Input to missense badness calculations.
"""

misbad = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(path=f"{MPC_PREFIX}/{version}/missense_badness.ht")
        for version in GNOMAD_VERSIONS
    },
)
"""
Table containing all possible amino acid substitutions and their missense badness scores.
"""

joint_clinvar_gnomad = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(path=f"{MPC_PREFIX}/{version}/joint_clinvar_gnomad.ht")
        for version in GNOMAD_VERSIONS
    },
)
"""
Table containing "population" and "pathogenic" variants.

Table contains common (AF > 0.01) gnomAD variants ("population") and
ClinVar pathogenic/likely pathogenic missense variants in haploinsufficient genes
that cause severe disease ("pathogenic") with defined CADD, BLOSUM, Grantham, missense observed/expected ratios,
missense badness, and PolyPhen-2 scores.

Input to MPC (missense badness, polyphen-2, and constraint) calculations.
"""

# MPC related resources
mpc_model_pkl_path = f"{MPC_PREFIX}/{CURRENT_VERSION}/mpc_model.pkl"
"""
Path to model (stored as pickle) that contains relationship of MPC variables.

Created using logistic regression.
"""

gnomad_fitted_score = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(path=f"{MPC_PREFIX}/{version}/gnomad_fitted_scores.ht")
        for version in GNOMAD_VERSIONS
    },
)
"""
Table of gnomAD variants and their fitted scores (from MPC model regression).

Input to MPC (missense badness, polyphen-2, and constraint) calculations on other datasets.
"""

gnomad_fitted_score_group = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(
            path=f"{MPC_PREFIX}/{version}/gnomad_fitted_scores_group.ht"
        )
        for version in GNOMAD_VERSIONS
    },
)
"""
Table of fitted scores for common (AF > 0.01) variants in gnomAD, grouped by score.

Annotated with the total number of variants with and less than each score.
"""

mpc_release = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(path=f"{MPC_PREFIX}/{version}/mpc.ht")
        for version in GNOMAD_VERSIONS
    },
)
"""
Table containing missense variants in canonical transcripts annotated with MPC.
"""

mpc_release_dedup = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(path=f"{MPC_PREFIX}/{version}/mpc_dedup.ht")
        for version in GNOMAD_VERSIONS
    },
)
"""
Table containing missense variants in canonical transcripts annotated with MPC.

This Table contains only one row per each unique locus/alleles combination.
"""

# Assessment related resources
oe_bin_counts_tsv = f"{CONSTRAINT_PREFIX}/{CURRENT_VERSION}/oe_bin.tsv"
"""
TSV with RMC regions grouped by obs/exp (OE) bin.

Annotated with proportion coding base pairs, proportion de novo missense (controls),
proportion de novo missense (case), and proportion ClinVar pathogenic/likely pathogenic
severe haploinsufficient missense.
"""

# Reference related resources
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
