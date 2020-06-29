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
CODON_TABLE_PATH = f"{RESOURCE_PREFIX}/amino_acids/codons_lookup.tsv"
ACID_NAMES_PATH = f"{RESOURCE_PREFIX}/amino_acids/acid_names.tsv"
MUTATION_RATE_TABLE_PATH = f"{RESOURCE_PREFIX}/GRCh37/exac/mutation_rate_table.tsv"
DIVERGENCE_SCORES_TSV_PATH = (
    f"{RESOURCE_PREFIX}/GRCh37/exac/divsites_gencodev19_all_transcripts.tsv"
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
mutation_rate = TableResource(
    path=f"{RESOURCE_PREFIX}/GRCh37/exac/ht/mutation_rate.ht",
    import_func=_import_mu,
    import_args={
        "path": MUTATION_RATE_TABLE_PATH,
        "min_partitions": 50,
        "impute": True,
    },
)


## Observed/expected count related resources
# Expected variants resource files
MODEL_PREFIX = "gs://regional_missense_constraint/model"
EXP_PREFIX = f"{MODEL_PREFIX}/exp/"
