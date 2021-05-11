from gnomad.resources.resource_utils import TableResource, VersionedTableResource

from ukbb_qc.resources.basics import release_ht_path

from rmc.resources.resource_utils import DATA_SOURCE, UKBB_VER, RESOURCE_PREFIX


UKBB_DIR = f"{RESOURCE_PREFIX}/GRCh38/ukbb/freeze_{UKBB_VER}"
"""
Path to UKBB resources bucket.
"""

processed_exomes = VersionedTableResource(
    default_version=UKBB_VER,
    versions={
        UKBB_VER: TableResource(path=f"{release_ht_path(DATA_SOURCE, UKBB_VER)}")
    },
)
"""
UKBB release sites Table.

Contains all release annotations, including VEP and frequency.
"""

filtered_exomes = VersionedTableResource(
    default_version=UKBB_VER,
    version={UKBB_VER: TableResource(path=f"{UKBB_DIR}/exomes_missense_only.ht")},
)
"""
UKBB exomes Table filtered to missense variants only.
"""

prop_obs_coverage = VersionedTableResource(
    default_version=UKBB_VER,
    version={
        UKBB_VER: TableResource(
            path=f"{UKBB_DIR}/prop_observed_by_coverage_no_common_pass_filtered_bins.ht"
        )
    },
)
"""
Table with proportion of variants observed by coverage.

Calculated using `get_proportion_observed_by_coverage`.
Input for `build_models`.
"""

possible_variants_ht = VersionedTableResource(
    default_version=UKBB_VER,
    versions={
        UKBB_VER: TableResource(
            path=f"{UKBB_DIR}/possible_data/possible_transcript_pop_standard.ht"
        )
    },
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
