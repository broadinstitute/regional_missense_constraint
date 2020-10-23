import hail as hl

from gnomad.resources.resource_utils import (
    import_sites_vcf,
    TableResource,
)
from rmc.resources.resource_utils import RESOURCE_PREFIX


def _import_coverage(**kwargs) -> hl.Table:
    """
	Imports ExAC coverage TSV into Table.

	:return: Table with ExAC coverage information
	:rtype: hl.Table
	"""
    coverage = hl.import_table(**kwargs)
    coverage = coverage.rename(
        {
            "1": "over_1",
            "5": "over_5",
            "10": "over_10",
            "15": "over_15",
            "20": "over_20",
            "25": "over_25",
            "30": "over_30",
            "50": "over_50",
            "100": "over_100",
        }
    )
    coverage = coverage.key_by(locus=hl.locus(coverage["#chrom"], coverage.pos))
    coverage = coverage.naive_coalesce(100)
    return coverage


# Files for direct comparison with Kaitlin's code
EXAC_PREFIX = f"{RESOURCE_PREFIX}/GRCh37/exac"

exac = TableResource(
    path=f"{EXAC_PREFIX}/ht/ExAC.r1.sites.vep.ht",
    import_func=import_sites_vcf,
    import_args={
        "path": f"{EXAC_PREFIX}/ExAC.r1.sites.vep.vcf.gz",
        "force_bgz": True,
        "min_partitions": 500,
        "reference_genome": "GRCh37",
    },
)
"""
Resource for full ExAC dataset
"""

coverage = TableResource(
    path=f"{EXAC_PREFIX}/ht/22_coverage.ht",
    import_func=_import_coverage,
    import_args={
        "path": "gs://gnomad-public/legacy/exacv1_downloads/release0.1/coverage/Panel.chr22.coverage.txt.gz",
        "min_partitions": 100,
        "impute": True,
        "force_bgz": True,
    },
)
"""
Resource with ExAC coverage
"""

filtered_exac = TableResource(path=f"{EXAC_PREFIX}/ht/ExAC.r1.missense_only.ht")
"""
ExAC dataset filtered to missense variants only on chromosome 22 (specifically MYH9, PI4KA, MAPK1).
Also annotated with trimer context and coverage information.
"""
