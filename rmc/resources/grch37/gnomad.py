from rmc.resources.resource_utils import RESOURCE_PREFIX, FLAGSHIP_LOF

## Konrad's resources
# LoF constraint resource files
FLAGSHIP_LOF = "gs://gnomad-public/papers/2019-flagship-lof/v1.0/"
FLAGSHIP_MODEL_PREFIX = f"{FLAGSHIP_LOF}/model/"

processed_exomes = TableResource(path=f"{FLAGSHIP_MODEL_PREFIX}/exomes_processed.ht")
"""
Processed gnomAD exomes Table.

Dropped colocated variants in vep annotation and removed all non-pass variants.
Also annotated with context Table (sequence context, transcript information, most severe consequence).
"""
filtered_exomes = TableResource(path=f"{RESOURCE_PREFIX}/ht/exomes_missense_only.ht")
"""
Processed gnomAD exomes Table filtered to missense variants only.
"""
processed_genomes = TableResource(path=f"{FLAGSHIP_MODEL_PREFIX}/genomes_processed.ht")
"""
Processed gnomAD genomes Table.

Dropped colocated variants in vep annotation and removed all non-pass variants.
Also annotated with context Table (sequence context, transcript information, most severe consequence).
"""
