import hail as hl

BUILDS = ["GRCh37", "GRCh38"]
FLAGSHIP_LOF = "gs://gnomad-public/papers/2019-flagship-lof/v1.0"
RMC_PREFIX = "gs://regional_missense_constraint"
RESOURCE_PREFIX = "gs://regional_missense_constraint/resources"
GNOMAD_VER = "2.1.1"
GNOMAD_VERS = ["2.1.1"]
UKBB_FREEZE = 5
UKBB_FREEZES = [5]


DATA_SOURCE = "broad"
"""
Source of UKBB sequencing data.

Original options were 'regeneron' or 'broad', but the data source will always be 'broad' for 
anything larger than the 100K callset.
"""


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
    gencode = gencode.filter(gencode.gene_type == "protein_coding")
    return gencode
