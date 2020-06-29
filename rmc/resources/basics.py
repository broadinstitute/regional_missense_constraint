import hail as hl

from gnomad.resources.resource_utils import TableResource
from rmc.resources.resource_utils import FLAGSHIP_LOF, RESOURCE_PREFIX


LOGGING_PATH = "gs://regional_missense_constraint/logs"
"""
Path to bucket that stores hail logs.
"""


def _import_mu(**kwargs) -> hl.Table:
    """
    Imports mutation rate information from Kaitlin into Table.

    :return: Table with context, alleles, and mutation rate
    :rtype: hl.Table
    """
    # from  n_kmer  p_any_snp_given_kmer    mu_kmer to  count_snp   p_snp_given_kmer    mu_snp
    mutation = hl.import_table(**kwargs)
    mutation = mutation.transmute(
        context=mutation["from"], alleles=[mutation["from"][1], mutation.to[1]]
    )
    return mutation.key_by("context", "alleles").select("mu_snp")


## Kaitlin's resources
# Original regional missense constraint resource files
CODON_TABLE_PATH = f"{RESOURCE_PREFIX}/codons_lookup.tsv"
ACID_NAMES_PATH = f"{RESOURCE_PREFIX}/acid_names.tsv"
MUTATION_RATE_TABLE_PATH = f"{RESOURCE_PREFIX}/mutation_rate_table.tsv"
DIVERGENCE_SCORES_TSV_PATH = (
    f"{RESOURCE_PREFIX}/divsites_gencodev19_all_transcripts.tsv"
)
divergence_scores = TableResource(
    path=f"{RESOURCE_PREFIX}/div_scores.ht",
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
mutation_rate = TableResource(
    path=f"{RESOURCE_PREFIX}/ht/mutation_rate.ht",
    import_func=_import_mu,
    import_args={
        "path": MUTATION_RATE_TABLE_PATH,
        "min_partitions": 50,
        "impute": True,
    },
)


## Konrad's resources
# LoF constraint resource files
FLAGSHIP_LOF = "gs://gnomad-public/papers/2019-flagship-lof/v1.0/"
FLAGSHIP_MODEL_PREFIX = f"{FLAGSHIP_LOF}/model/"

processed_exomes = TableResource(path=f"{FLAGSHIP_MODEL_PREFIX}/exomes_processed.ht",)
"""
Processed gnomAD exomes Table.

Dropped colocated variants in vep annotation and removed all non-pass variants.
Also annotated with context Table (sequence context, transcript information, most severe consequence).
"""
filtered_exomes = TableResource(path=f"{RESOURCE_PREFIX}/ht/exomes_missense_only.ht",)
"""
Processed gnomAD exomes Table filtered to missense variants only.
"""
processed_genomes = TableResource(path=f"{FLAGSHIP_MODEL_PREFIX}/genomes_processed.ht",)
"""
Processed gnomAD genomes Table.

Dropped colocated variants in vep annotation and removed all non-pass variants.
Also annotated with context Table (sequence context, transcript information, most severe consequence).
"""


## Observed/expected count related resources
# Expected variants resource files
MODEL_PREFIX = "gs://regional_missense_constraint/model"
EXP_PREFIX = f"{MODEL_PREFIX}/exp/"
