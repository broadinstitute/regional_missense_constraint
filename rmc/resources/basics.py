import hail as hl

from gnomad.resources.resource_utils import TableResource, VersionedTableResource

from rmc.resources.resource_utils import (
    CURRENT_VERSION,
    FLAGSHIP_LOF,
    GNOMAD_VERSIONS,
    RESOURCE_PREFIX,
    RMC_PREFIX,
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
MUTATION_RATE_TABLE_PATH = f"{RESOURCE_PREFIX}/GRCh37/exac/mutation_rate_table.tsv"
DIVERGENCE_SCORES_TSV_PATH = (
    f"{RESOURCE_PREFIX}/GRCh37/exac/divsites_gencodev19_all_transcripts.txt"
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

## Amino acid resources
CODON_TABLE_PATH = f"{RESOURCE_PREFIX}/amino_acids/codons_lookup.tsv"
"""
TSV file containing two columns: codon and three letter amino acid code.
"""

ACID_NAMES_PATH = f"{RESOURCE_PREFIX}/amino_acids/acid_names.tsv"
"""
TSV file containing three columns: amino acid name, amino acid 3 letter code, and amino acid 1 letter code.
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

## gnomAD resources
mutation_rate = TableResource(
    path=f"{FLAGSHIP_LOF}/model/mutation_rate_methylation_bins.ht",
)
"""
Table with mutation rate recalculated for gnomAD constraint.

This was calculated with `calculate_mu_by_downsampling` in
https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/constraint_basics.py.
"""

## Gene/Transcript related resources
MODEL_PREFIX = f"{RMC_PREFIX}/model"
"""
Path to bucket containing resources related to building the mutational models.

Bucket also contains transcript-related resources.
"""

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

hi_genes = f"{RESOURCE_PREFIX}/HI_genes.rCNV.txt"
"""
Path to haploinsufficient genes that cause severe disease.

List is from Ryan Collins.
"""

## Constraint related resources
CONSTRAINT_PREFIX = f"{RMC_PREFIX}/constraint"
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

simul_break_temp = f"{temp_path}/simul_breaks"
"""
Bucket to store temporary results for simultaneous results.
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

amino_acids_oe = VersionedTableResource(
    default_version=CURRENT_VERSION,
    versions={
        version: TableResource(path=f"{MODEL_PREFIX}/{version}/amino_acid_oe.ht")
        for version in GNOMAD_VERSIONS
    },
)
"""
Table containing all possible amino acid substitutions and their missense OE ratio.

Input to missense badness calculations.
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


oe_bin_counts_tsv = f"{CONSTRAINT_PREFIX}/{CURRENT_VERSION}/oe_bin.tsv"
"""
TSV with RMC regions grouped by obs/exp (OE) bin.

Annotated with proportion coding base pairs, proportion de novo missense (controls),
proportion de novo missense (case), and proportion ClinVar pathogenic/likely pathogenic
severe haploinsufficient missense.
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
