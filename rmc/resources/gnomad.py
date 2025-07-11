"""Script containing resources for the current gnomAD version."""
from gnomad.resources.resource_utils import TableResource, VersionedTableResource

constraint_ht = VersionedTableResource(
    default_version="4.1.1",
    versions={
        "2.1.1": TableResource(
            path="gs://gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.ht"
        ),
        "4.1": TableResource(
            path="gs://gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.ht"
        ),
        # TODO: Update this path when the new public version is released.
        "4.1.1": TableResource(
            path="gs://gnomad/v4.1/constraint_coverage_corrected/metrics/transcript_consequences/gnomad.v4.1.constraint_metrics.coverage_corrected.ht"
        ),
    },
)
"""
Public gnomAD gene constraint Table.

v2.1.1: Observed variants count is annotated as `ht.obs_mis` and expected variants count is annotated as `ht.exp_mis`.
v4.1: Observed variants count is annotated as `ht.mis.obs` and expected variants count is annotated as `ht.mis.exp`.
NOTE: default version is manually set to point to most updated resource. When public resource is updated,
please set the default version to CURRENT_GNOMAD_VERSION and add this import
from rmc.resources.resource_utils import CURRENT_GNOMAD_VERSION
"""
