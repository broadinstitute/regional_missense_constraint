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
    coverage = coverage.transmute(
        locus=hl.parse_locus(hl.format("%s:%s", coverage["#chrom"], coverage.pos))
    )
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
    coverage = coverage.key_by("locus")
    coverage = coverage.naive_coalesce(100)
    return coverage


# Files for direct comparison with Kaitlin's code
EXAC_PREFIX = f"{RESOURCE_PREFIX}/ExAC"

exac = TableResource(
    path=f"{EXAC_PREFIX}/ExAC.r0.3.sites.vep.vcf.gz",
    import_func=import_sites_vcf,
    import_args={
        "path": f"{EXAC_PREFIX}/ExAC.r0.3.sites.vep.vcf.gz",
        "force_bgz": True,
        "min_partitions": 500,
        "reference_genome": "GRCh37",
    },
)
"""
Resource for full ExAC dataset
"""

coverage = TableResource(
    path=f"{EXAC_PREFIX}/coverage/22_coverage.ht",
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

filtered_exac = TableResource(path=f"{EXAC_PREFIX}/ExAC.r0.3.missense_only.ht")
"""
Resource for ExAC dataset filtered to missense variants only and annotated with trimer context
"""

filtered_exac_cov = TableResource(path=f"{EXAC_PREFIX}/ExAC.r0.3.missense_only_cov.ht")
"""
ExAC dataset filtered to missense variants only and annotated with trimer context and coverage information.
Contains coverage infromation for ONLY chr22.
"""


# Observed/expected count related resources
MODEL_PREFIX = "gs://regional_missense_constraint/model"
EXP_PREFIX = f"{MODEL_PREFIX}/exp/"

expected = TableResource(path=f"{EXP_PREFIX}/ExAC_exp_var.ht")
"""
Table annotated with expected variant counts
"""
