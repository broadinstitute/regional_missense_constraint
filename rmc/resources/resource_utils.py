"""Script containing resource utility constants."""
CURRENT_BUILD = "GRCh37"
BUILDS = ["GRCh37"]

CURRENT_GNOMAD_VERSION = "2.1.1"
GNOMAD_VERSIONS = ["2.1.1"]

MISSENSE = "missense_variant"
"""
String representing missense variant VEP annotation.
"""
NONSENSES = {"stop_gained", "splice_donor_variant", "splice_acceptor_variant"}
"""
Set containing loss-of-function (LoF) single-nucleotide variant VEP annotations.
"""
READ_THROUGH = "stop_lost"
"""
String representing read-through (stop lost) variant VEP annotation.
"""
SYNONYMOUS = "synonymous_variant"
"""
String representing synonymous variant VEP annotation.
"""
KEEP_CODING_CSQ = {MISSENSE, READ_THROUGH, SYNONYMOUS}.union(NONSENSES)
"""
Set of variant consequences to keep in filtered context Table.
"""
