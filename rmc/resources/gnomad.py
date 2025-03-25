"""Script containing resources for the current gnomAD version."""
from gnomad.resources.resource_utils import TableResource, VersionedTableResource

from rmc.resources.resource_utils import CURRENT_GNOMAD_VERSION

constraint_ht = VersionedTableResource(
    default_version=CURRENT_GNOMAD_VERSION,
    versions={
        "2.1.1": TableResource(
            path="gs://gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.ht"
        ),
        "4.1": TableResource(
            # TODO: change this back when new constraint metrics are made public
            # path="gs://gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.ht"
            path="gs://gnomad/v4.1/constraint_an_3_10_25/metrics/transcript_consequences/gnomad.v4.1.constraint_metrics.ht"
        ),
    },
)
"""
Public gnomAD gene constraint Table.

v2.1.1: Observed variants count is annotated as `ht.obs_mis` and expected variants count is annotated as `ht.exp_mis`.
v4.1: Observed variants count is annotated as `ht.mis.obs` and expected variants count is annotated as `ht.mis.exp`.
"""
