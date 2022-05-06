import logging

import hail as hl

from gnomad.resources.grch37.gnomad import coverage, public_release
from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists
from gnomad.utils.vep import CSQ_NON_CODING

from rmc.resources.basics import (
    amino_acids_oe,
    constraint_prep,
    misbad,
    rmc_results,
    temp_path,
)
from rmc.resources.grch37.gnomad import constraint_ht
from rmc.utils.generic import (
    filter_context_using_gnomad,
    get_codon_lookup,
    get_constraint_transcripts,
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

    'coding_sequence_variant' defined as: 'A sequence variant that changes the coding sequence'
    https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html
    :param hl.Table ht: Input Table.
    :return: Table with informative codons only.
    """
    logger.info("Removing non-coding loci from HT...")
    non_coding_csq = hl.literal(CSQ_NON_CODING)
    ht = ht.filter(~non_coding_csq.contains(ht.most_severe_consequence))

    logger.info("Filtering to lines with expected codon annotations...")
    # Codons are in this format: NNN/NNN, so expected length is 7
    ht = ht.filter((hl.is_defined(ht.codons)) & (hl.len(ht.codons) == 7))
    codon_map = get_codon_lookup()
    ht = ht.annotate(
        ref=ht.codons.split("/")[0].upper(),
        alt=ht.codons.split("/")[1].upper(),
    )
    ht = ht.annotate(
        ref=codon_map.get(ht.ref, "Unk"),
        alt=codon_map.get(ht.alt, "Unk"),
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
    # Recalculating transcript level OE ratio because previous OE ratio (`overall_oe`)
    # is capped at 1 for regional missense constraint calculation purposes
    group_ht = group_ht.annotate(transcript_oe=group_ht.obs / group_ht.exp)

    # Read in LoF constraint HT to get OE ratio for five transcripts missing in v2 RMC results
    # # 'ENST00000304270', 'ENST00000344415', 'ENST00000373521', 'ENST00000381708', 'ENST00000596936'
    # All 5 of these transcripts have extremely low coverage in gnomAD
    # Will keep for consistency with v2 LoF results but they look terrible
    # NOTE: LoF HT is keyed by gene and transcript, but `_key_by_assert_sorted` doesn't work here for v2 version
    # Throws this error: hail.utils.java.FatalError: IllegalArgumentException
    lof_ht = constraint_ht.ht().select("oe_mis").key_by("transcript")
    ht = ht.annotate(
        gnomad_transcript_oe=lof_ht[ht.transcript].oe_mis,
        transcript_oe=group_ht[ht.transcript].oe,
    )
    ht = ht.transmute(transcript_oe=hl.coalesce(ht.rmc_oe, ht.gnomad_oe))

    rmc_ht = rmc_results.ht().select("section_oe")
    ht = ht.annotate(section_oe=rmc_ht[ht.locus, ht.transcript].section_oe)
    return ht.transmute(oe=hl.coalesce(ht.section_oe, ht.transcript_oe))


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
        Must be one of "exomes" or "genomes" (check is done within `public_release`).
        Default is "exomes".
    :return: None; writes amino acid Table to resource path.
    """
    logger.info("Importing set of transcripts to keep...")
    transcripts = get_constraint_transcripts(outlier=False)

    logger.info("Reading in VEP context HT...")
    # NOTE: Keeping all variant types here because need synonymous and nonsense variants to calculate missense badness
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
    # context_ht = context_ht.checkpoint(f"{temp_path}/codons.ht", overwrite=True)
    context_ht = hl.read_table(f"{temp_path}/codons.ht")

    logger.info("Filtering sites using gnomAD %s...", gnomad_data_type)
    context_ht = filter_context_using_gnomad(context_ht, gnomad_data_type)
    logger.info("Adding observed annotation...")
    gnomad = public_release(gnomad_data_type).ht()
    gnomad_cov = coverage(gnomad_data_type).ht()
    gnomad = gnomad.select(
        "filters",
        ac=gnomad.freq[0].AC,
        af=gnomad.freq[0].AF,
        gnomad_coverage=gnomad_cov[gnomad.locus].median,
    )
    gnomad = gnomad.filter(
        keep_criteria(gnomad.ac, gnomad.af, gnomad.filters, gnomad.gnomad_coverage)
    )
    context_ht = context_ht.annotate(_obs=gnomad.index(context_ht.key))
    context_ht = context_ht.transmute(observed=hl.int(hl.is_defined(context_ht._obs)))

    logger.info("Checkpointing HT after joining with gnomAD data...")
    # context_ht = context_ht.checkpoint(f"{temp_path}/codons_filt.ht", overwrite=True)
    context_ht = hl.read_table(f"{temp_path}/codons_filt.ht")

    logger.info(
        "Getting observed to expected ratio, rekeying Table, and writing to output path..."
    )
    # Note that `get_oe_annotation` is pulling the missense OE ratio
    context_ht = get_oe_annotation(context_ht)
    context_ht = context_ht.key_by().select(
        "ref",
        "alt",
        "observed",
        "codons",
        "amino_acids",
        "oe",
    )
    context_ht.write(amino_acids_oe.path, overwrite=True)


def variant_csq_expr(
    ref_expr: hl.expr.StringExpression, alt_expr: hl.expr.StringExpression
) -> hl.expr.StringExpression:
    """
    Determine variant consequence using reference and alternate amino acid annotations.

    :param hl.expr.StringExpression ref_expr: Reference amino acid StringExpression.
    :param hl.expr.StringExpression alt_expr: Alternate amino acid StringExpression.
    :return: Variant type StringExpression. One of 'syn', 'non', 'mis', 'rdt' (stop lost).
    """
    return (
        hl.case()
        .when(ref_expr == alt_expr, "syn")
        .when(alt_expr == "STOP", "non")
        .when(ref_expr == "STOP", "rdt")
        .default("mis")
    )


def split_ht_by_oe(
    ht: hl.Table, keep_high_oe: bool, oe_threshold: float = 0.6
) -> hl.Table:
    """
    Split Table with all possible amino acid substitutions based on missense observed to expected (OE) ratio cutoff.

    :param hl.Table ht: Input Table with amino acid substitutions.
    :param bool keep_high_oe: Whether to filter to high missense OE values.
        If True, returns "boring" HT.
        If False, gets "bad" (low missense OE) Table.
    :param float oe_threshold: OE Threshold used to split Table.
        Rows with OE less than this threshold will be filtered if `keep_high_oe` is True, and
        rows with OE greater than or equal to this threshold will be kept.
        Default is 0.6.
    :return: Table filtered based on missense OE. Schema:
        ----------------------------------------
        Row fields:
            'ref': str
            'alt': str
            'oe': float64
            'obs': int64
            'possible': int64
            'mut_type': str
        ----------------------------------------
    """
    logger.info("Filtering HT on missense OE values...")
    oe_filter_expr = (ht.oe > oe_threshold) if keep_high_oe else (ht.oe <= oe_threshold)
    ht = ht.filter(oe_filter_expr)

    logger.info("Grouping HT and aggregating observed and possible variant counts...")
    ht = ht.group_by("ref", "alt").aggregate(
        obs=hl.agg.sum(ht.observed), possible=hl.agg.count()
    )

    logger.info("Adding variant consequence (mut_type) annotation and returning...")
    return ht.annotate(mut_type=variant_csq_expr(ht.ref, ht.alt))


def get_total_csq_count(ht: hl.Table, csq: str, count_field: str) -> int:
    """
    Filter input Table using specified variant consequence and aggregate total value of specified field.

    :param hl.Table ht: Input Table (Table with amino acid substitutions filtered to have high or low missense OE).
    :param str csq: Desired variant consequence. One of "syn" or "non".
    :param str count_field: Desired count type. One of "obs" or "possible".
    :return: Int of total value of `count_field` for specified consequence.
    """
    return ht.aggregate(hl.agg.filter(ht.mut_type == csq, hl.agg.sum(ht[count_field])))


def calculate_misbad(use_exac_oe_cutoffs: bool) -> None:
    """
    Calculate missense badness score using Table with all amino acid substitutions and their missense observed/expected (OE) ratio.

    If `use_exac_oe_cutoffs` is set, will remove all rows with 0.6 < OE <= 0.8.

    :param bool use_exac_oe_cutoffs: Whether to use the same missense OE cutoffs as in ExAC missense badness calculation.
    :return: None; writes Table with missense badness score to resource path.
    """
    if not file_exists(amino_acids_oe.path):
        raise DataException(
            "Table with all amino acid substitutions and missense OE doesn't exist!"
        )

    ht = amino_acids_oe.ht()
    if use_exac_oe_cutoffs:
        logger.info("Removing rows with OE greater than 0.6 and less than 0.8...")
        ht = ht.filter((ht.oe <= 0.6) | (ht.oe > 0.8))

    logger.info(
        "Splitting input Table by OE to get synonymous and nonsense rates for high and low OE groups..."
    )
    logger.info("Creating high missense OE (OE > 0.6) HT...")
    high_ht = split_ht_by_oe(ht, keep_high_oe=True)
    high_ht = high_ht.checkpoint(f"{temp_path}/amino_acids_high_oe.ht", overwrite=True)

    logger.info("Creating low missense OE (OE <= 0.6) HT...")
    low_ht = split_ht_by_oe(ht, keep_high_oe=False)
    low_ht = low_ht.checkpoint(f"{temp_path}/amino_acids_low_oe.ht", overwrite=True)

    logger.info("Calculating synonymous rates...")
    syn_obs_high = get_total_csq_count(high_ht, csq="syn", count_field="obs")
    syn_pos_high = get_total_csq_count(high_ht, csq="syn", count_field="possible")
    syn_obs_low = get_total_csq_count(low_ht, csq="syn", count_field="obs")
    syn_pos_low = get_total_csq_count(low_ht, csq="syn", count_field="possible")
    syn_rate = (syn_obs_high / syn_pos_high) / (syn_obs_low / syn_pos_low)
    logger.info("Synonymous rate: %i", syn_rate)

    logger.info("Calculating nonsense rates...")
    non_obs_high = get_total_csq_count(high_ht, csq="non", count_field="obs")
    non_pos_high = get_total_csq_count(high_ht, csq="non", count_field="possible")
    non_obs_low = get_total_csq_count(low_ht, csq="non", count_field="obs")
    non_pos_low = get_total_csq_count(low_ht, csq="non", count_field="possible")
    non_rate = (non_obs_high / non_pos_high) / (non_obs_low / non_pos_low)
    logger.info("Nonsense rate: %i", non_rate)

    logger.info("Re-joining split HTs to calculate missense badness...")
    high_ht = high_ht.transmute(
        high_obs=high_ht.obs,
        high_pos=high_ht.possible,
    )
    low_ht = low_ht.transmute(
        low_obs=low_ht.obs,
        low_pos=low_ht.possible,
    )
    ht = high_ht.join(low_ht, how="outer")
    mb_ht = ht.group_by("ref", "alt").aggregate(
        high_low=(
            (hl.agg.sum(ht.high_obs) / hl.agg.sum(ht.high_pos))
            / (hl.agg.sum(ht.low_obs) / hl.agg.sum(ht.low_pos))
        )
    )
    mb_ht = mb_ht.annotate(mut_type=variant_csq_expr(mb_ht.ref, mb_ht.alt))
    mb_ht = mb_ht.annotate(
        misbad=hl.or_missing(
            mb_ht.mut_type == "mis",
            # Cap missense badness at 1
            hl.min((mb_ht.high_low - syn_rate) / (non_rate - syn_rate), 1),
        ),
    )
    mb_ht.write(misbad.path, overwrite=True)
