"""Script containing resources for the current gnomAD version."""
from gnomad.resources.resource_utils import TableResource, VersionedTableResource

from rmc.resources.basics import RESOURCE_BUILD_PREFIX
from rmc.resources.resource_utils import CURRENT_GNOMAD_VERSION

FLAGSHIP_LOF = "gs://gnomad-public-requester-pays/papers/2019-flagship-lof/v1.0"
"""
Path to bucket with gnomAD v2 loss-of-function (LoF) constraint results.
"""

FLAGSHIP_LOF_MODEL_PREFIX = f"{FLAGSHIP_LOF}/model"

processed_exomes = TableResource(
    path=f"{FLAGSHIP_LOF_MODEL_PREFIX}/exomes_processed.ht"
)
"""
Processed gnomAD exomes Table.

Dropped colocated variants in vep annotation and removed all non-pass variants.
Also annotated with context Table (sequence context, transcript information, most severe consequence).
"""
filtered_exomes = TableResource(
    path=f"{RESOURCE_BUILD_PREFIX}/gnomad/ht/exomes_missense_only.ht"
)
"""
Processed gnomAD exomes Table filtered to rare (AF < 0.001) missense variants in canonical transcripts only.

NOTE: This resource is created with `process_vep`, which now filters to non-outlier transcripts by default.
However, for v2, this resource contains *all* canonical transcripts (including outliers).
"""

processed_genomes = TableResource(
    path=f"{FLAGSHIP_LOF_MODEL_PREFIX}/genomes_processed.ht"
)
"""
Processed gnomAD genomes Table.

Dropped colocated variants in vep annotation and removed all non-pass variants.
Also annotated with context Table (sequence context, transcript information, most severe consequence).
"""

prop_obs_coverage = TableResource(
    path=f"{FLAGSHIP_LOF_MODEL_PREFIX}/prop_observed_by_coverage_no_common_pass_filtered_bins.ht"
)
"""
Table with proportion of variants observed by coverage.

Calculated using `get_proportion_observed_by_coverage`.
Input for `build_models`.
"""

possible_variants_ht = TableResource(
    path=f"{FLAGSHIP_LOF_MODEL_PREFIX}/possible_data/possible_transcript_pop_standard.ht"
)
"""
Table with all observed SNPs in hg19 fasta (context) Table.

Calculated using `get_proportion_observed`.

Contains multiple mutation rate annotations:
	- mu_snp: Raw mutation rate calculated using gnomAD v2 genomes.
	- mu_agg: `mu_snp` multiplied by the number of times the variant was seen in the context HT (`possible_variants`).
	- adjusted_mutation_rate: `mu_agg` corrected with plateau model.
	- mu: `mu_agg` multipled by coverage correction.
"""

mutation_rate = TableResource(
    path=f"{FLAGSHIP_LOF_MODEL_PREFIX}/mutation_rate_methylation_bins.ht",
)
"""
Table with mutation rate recalculated for gnomAD constraint.

This was calculated with `calculate_mu_by_downsampling` in
https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/constraint_basics.py.
"""

constraint_ht = VersionedTableResource(
    default_version=CURRENT_GNOMAD_VERSION,
    versions={
        CURRENT_GNOMAD_VERSION: TableResource(
            path="gs://gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.ht"
        )
    },
)
"""
Table with total observed and expected variant counts per transcript.

Calculated as part of gnomAD LoF release for v2.1.1.
Observed variants count is annotated as `obs_mis` and expected variants count is annotated as `exp_mis`.
"""
