from rmc.resources.resource_utils import RESOURCE_PREFIX, FLAGSHIP_LOF

## Kaitlin's resources
# Original regional missense constraint resource files
codon_table_path = f"{RESOURCE_PREFIX}/codons_lookup.tsv"
acid_names_path = f"{RESOURCE_PREFIX}/acid_names.tsv"
mutation_rate_table_path = f"{RESOURCE_PREFIX}/mutation_rate_table.tsv"
divergence_scores_path = f"{RESOURCE_PREFIX}/divsites_gencodev19_all_transcripts.tsv"
divergence_ht = f"{RESOURCE_PREFIX}/div_scores.ht"

## Konrad's resources
# LoF constraint resource files
FLAGSHIP_LOF = "gs://gnomad-public/papers/2019-flagship-lof/v1.0/"
MODEL_PREFIX = f"{FLAGSHIP_LOF}/model/"
processed_exomes_ht_path = f"{MODEL_PREFIX}/exomes_processed.ht"
filtered_exomes_ht_path = f"{RESOURCE_PREFIX}/ht/exomes_missense_only.ht"
processed_genomes_ht_path = f"{MODEL_PREFIX}/genomes_processed.ht"

# Processed constraint resource files
mutation_rate_ht = f"{RESOURCE_PREFIX}/ht/mutation_rate.ht"