import hail as hl

from gnomad.resources.grch38.reference_data import vep_context
from gnomad.resources.resource_utils import (
    TableResource,
    VersionedTableResource,
)

from rmc.resources.resource_utils import import_gencode, RESOURCE_PREFIX

## Reference genome related resources
full_context = VersionedTableResource(
    default_version="101", versions={"101": vep_context.versions["101"]},
)

processed_context = VersionedTableResource(
    default_version="101",
    versions={
        "101": TableResource(
            path=f"{RESOURCE_PREFIX}/GRCh38/reference_data/ht/context_fasta_snps_only.v101.ht",
        )
    },
)

## Exon/transcript related resources
gencode = VersionedTableResource(
    default_version="v30",
    versions={
        "v30": TableResource(
            path=f"{RESOURCE_PREFIX}/ht/context/gencode.v30.basic.annotation.ht",
            import_func=import_gencode,
            import_args={
                "path": f"{RESOURCE_PREFIX}/gencode.v30.basic.annotation.gtf",
                "reference_genome": "GRCh38",
                "skip_invalid_contigs": True,
                "min_partitions": 1000,
            },
        )
    },
)

processed_gencode = VersionedTableResource(
    default_version="v30",
    versions={
        "v30": TableResource(
            path=f"{RESOURCE_PREFIX}/ht/context/gencode.v30.basic.exons.ht",
        )
    },
)
