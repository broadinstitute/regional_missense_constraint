import hail as hl


FLAGSHIP_LOF = "gs://gnomad-public/papers/2019-flagship-lof/v1.0/"
RESOURCE_PREFIX = "gs://regional_missense_constraint/resources"


# Missense variant VEP annotations
MISSENSE = [
    "stop_lost",
    "initiator_codon_variant",
    "start_lost",
    "protein_altering_variant",
    "missense_variant",
]


# Import related resources
def import_gencode(**kwargs) -> hl.Table:
    """
	Converts Gencode GTF to Table.

	:return: Table
	:rtype: hl.Table
	"""
    gencode = hl.experimental.import_gtf(**kwargs)
    gencode = gencode.filter(
        (gencode.feature == "exon")
        & (gencode.gene_type == "protein_coding")
        & (gencode.level != "3")
    )
    return gencode
