import logging

import hail as hl

from gnomad.resources.grch37.gnomad import public_release
from gnomad.utils.vep import CSQ_NON_CODING

from rmc.resources.basics import amino_acids_oe, constraint_prep, rmc_results, temp_path
from rmc.resources.grch37.gnomad import constraint_ht
from rmc.utils.generic import (
    filter_context_using_gnomad,
    get_codon_lookup,
    get_outlier_transcripts,
    keep_criteria,
    process_context_ht,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("calculate_missense_badness")
logger.setLevel(logging.INFO)


def filter_codons(ht: hl.Table) -> hl.Table:
    """
    Remove non-coding loci and keep informative codons only.

    Remove rows with unknown amino acids. This also removes rows that are annotated as
    'coding_sequence_variant', as these variants have either undefined or uninformative codon annotations
    (NA or codon with Ns, e.g. nnG/nnT).

    'coding_sequence_variant' defined as: 'At sequence variant that changes the coding sequence'
    https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html
    :param hl.Table ht: Input Table.
    :return: Table with informative codons only.
    """
    logger.info("Removing non-coding loci from HT...")
    non_coding_csq = hl.literal(CSQ_NON_CODING)
    ht = ht.filter(
        ~non_coding_csq.contains(ht.transcript_consequences.most_severe_consequence)
    )

    logger.info("Filtering to lines with expected codon annotations...")
    # Codons are in this format: NNN/NNN, so expected length is 7
    ht = ht.filter((hl.is_defined(ht.codons)) & (hl.len(ht.codons) == 7))
    codon_map = get_codon_lookup()
    ht = ht.annotate(
        ref=ht.codons.split("/")[0].upper(), alt=ht.codons.split("/")[1].upper(),
    )
    ht = ht.annotate(
        ref=codon_map.get(ht.ref, "Unk"), alt=codon_map.get(ht.alt, "Unk"),
    )
    # Remove any lines with "Unk" (unknown amino acids)
    return ht.filter((ht.ref != "Unk") & (ht.alt != "Unk"))


def get_oe_annotation(ht: hl.Table) -> hl.Table:
    """
    Annotate input Table with observed to expected missense (OE) ratio per transcript.

    Use regional OE value if available, otherwise use transcript OE value.

    :param hl.Table ht: Input Table.
    :return: Table with `oe` annotation.
    """
    overall_oe_ht = (
        constraint_prep.ht().select_globals().select("total_obs", "total_exp")
    )
    group_ht = overall_oe_ht.group_by("transcript").aggregate(
        obs=hl.agg.take(overall_oe_ht.total_obs, 1)[0],
        exp=hl.agg.take(overall_oe_ht.total_exp, 1)[0],
    )
    group_ht = group_ht.annotate(oe=group_ht.obs / group_ht.exp)

    # Read in LoF constraint HT to get OE ratio for five transcripts missing in v2 RMC results
    # # 'ENST00000304270', 'ENST00000344415', 'ENST00000373521', 'ENST00000381708', 'ENST00000596936'
    # All 5 of these transcripts have extremely low coverage in gnomAD
    # Will keep for consistency with v2 LoF results but they look terrible
    lof_ht = constraint_ht.ht().select("oe_mis")
    ht = ht.annotate(
        gnomad_oe=constraint_ht[ht.transcript].oe_mis,
        rmc_oe=group_ht[ht.transcript].oe,
    )
    ht = ht.transmute(overall_oe=hl.coalesce(ht.rmc_oe, ht.gnomad_oe))

    rmc_ht = rmc_results.ht().select("section_oe")
    ht = ht.annotate(rmc_oe=rmc_ht[ht.locus, ht.transcript].section_oe)
    return ht.transmute(oe=hl.coalesce(ht.rmc_oe, ht.overall_oe))


def prepare_amino_acid_ht(gnomad_data_type: str = "exomes") -> None:
    """
    Prepare Table with all possible amino acid substitutions and their missense observed to expected (OE) ratio.

    Steps:
        - Import VEP context Table and filter to keep every possible amino acid substitution
        (every codon > codon change).
        - Filter Table to rows that aren't present in gnomAD or are rare in gnomAD (using `keep_criteria`).
        - Add observed and OE annotation
        - Write to `amino_acids_oe` resource path

    :param str gnomad_data_type: gnomAD data type. Used to retrieve public release and coverage resources.
        Default is "exomes".
    :return: None; writes amino acid Table to resource path.
    """
    logger.info("Importing set of transcripts to keep...")
    transcripts = get_outlier_transcripts(keep=True)

    logger.info("Reading in VEP context HT...")
    context_ht = process_context_ht(
        build="GRCh37", filter_to_missense=False, add_annotations=False
    )

    logger.info(
        "Filtering to transcripts to keep and selecting relevant annotations..."
    )
    context_ht = context_ht.filter(
        transcripts.contains(context_ht.transcript_consequences.transcript_id)
    )
    context_ht = context_ht.select(
        transcript=context_ht.transcript_consequences.transcript_id,
        consequence_terms=context_ht.transcript_consequences.consequence_terms,
        most_severe_consequence=context_ht.transcript_consequences.most_severe_consequence,
        amino_acids=context_ht.transcript_consequences.amino_acids,
        codons=context_ht.transcript_consequences.codons,
    )

    logger.info(
        "Filtering non-coding rows and rows with uninformative/unknown codons..."
    )
    context_ht = filter_codons(context_ht)

    logger.info("Checkpointing HT before joining with gnomAD data...")
    context_ht = context_ht.checkpoint(f"{temp_path}/codons.ht", overwrite=True)

    logger.info("Filtering sites using gnomAD %i...", gnomad_data_type)
    context_ht = filter_context_using_gnomad(context_ht)

    logger.info("Adding observed annotation...")
    gnomad = public_release(gnomad_data_type)
    gnomad = gnomad.filter(keep_criteria(gnomad))
    context_ht = context_ht.annotate(_obs=gnomad.index(context_ht.key))
    context_ht = context_ht.transmute(observed=hl.int(hl.is_defined(context_ht._obs)))

    logger.info("Checkpointing hT after joining with gnomAD data...")
    context_ht = context_ht.checkpoint(f"{temp_path}/codons_filt.ht", overwrite=True)

    logger.info(
        "Getting observed to expected missense ratio, rekeying Table, and writing to output path..."
    )
    context_ht = get_oe_annotation(context_ht)
    context_ht = context_ht.key_by().select(
        "ref", "alt", "observed", "codons", "amino_acids", "oe",
    )
    context_ht.write(amino_acids_oe.path, overwrite=True)
