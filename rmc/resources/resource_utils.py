import hail as hl

BUILDS = ["GRCh37", "GRCh38"]
FLAGSHIP_LOF = "gs://gnomad-public/papers/2019-flagship-lof/v1.0"
RESOURCE_PREFIX = "gs://regional_missense_constraint/resources"
GNOMAD_VER = "2.1.1"


MISSENSE = "missense_variant"
"""
String representing missense variant VEP annotation.
"""


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
