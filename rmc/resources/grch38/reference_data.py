import hail as hl

from gnomad.resources.grch38.reference_data import vep_context
from gnomad.resources.resource_utils import (
    TableResource,
    VersionedTableResource,
)

from rmc.resources.resource_utils import RESOURCE_PREFIX

## Reference genome related resources
full_context = VersionedTableResource(
    default_version="101",
    versions={"101": vep_context.versions["101"]},
)

filtered_context = VersionedTableResource(
    default_version="101",
    versions={
        "101": TableResource(
            path=f"{RESOURCE_PREFIX}/GRCh38/reference_data/ht/context_fasta_snps_only.v101.ht",
        )
    },
)
