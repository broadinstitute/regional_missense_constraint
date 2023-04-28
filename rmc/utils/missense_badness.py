import logging

import hail as hl
from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists

from rmc.resources.basics import TEMP_PATH_WITH_FAST_DEL
from rmc.resources.reference_data import (
    fold_k,
    test_transcripts_path,
    training_transcripts_path,
)
from rmc.resources.rmc import CURRENT_FREEZE, amino_acids_oe_path, misbad
from rmc.utils.constraint import add_obs_annotation, get_oe_annotation
from rmc.utils.generic import (
    annotate_and_filter_codons,
    filter_context_to_transcript_cds,
    filter_context_using_gnomad,
    process_context_ht,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("calculate_missense_badness")
logger.setLevel(logging.INFO)


def prepare_amino_acid_ht(
    overwrite_temp: bool,
    overwrite_output: bool,
    use_test_transcripts: bool = False,
    do_k_fold_training: bool = False,
    freeze: int = CURRENT_FREEZE,
    gnomad_data_type: str = "exomes",
    loftee_hc_str: str = "HC",
) -> None:
    """
    Prepare Table with all possible amino acid substitutions and their missense observed to expected (OE) ratio.

    Steps:
        - Import VEP context Table and filter to keep every possible amino acid substitution
        (every codon > codon change).
        - Filter Table to rows that aren't present in gnomAD or are rare in gnomAD (using `keep_criteria`).
        - Add observed and OE annotation.
        - Filter to specified transcript set(s).
        - Write to `amino_acids_oe_path` resource path(s).

    :param overwrite_temp: Whether to overwrite intermediate temporary (OE-independent) data if it already exists.
        If False, will read existing intermediate temporary data rather than overwriting.
    :param overwrite_output: Whether to entirely overwrite final output (OE-dependent) data if it already exists.
        If False, will read and modify existing output data by adding or modifying columns rather than overwriting entirely.
        If True, will clear existing output data and write new output data.
        The output Table is the amino acid Table.
    :param use_test_transcripts: Whether to use test transcripts in calculations.
        If False, will use training transcripts only.
        If True, will use test transcripts only.
        Default is False.
    :param do_k_fold_training: Whether to generate k-fold models with the training transcripts.
        If False, will use all training transcripts in calculation of a single model.
        If True, will calculate k models corresponding to the k-folds of the training transcripts.
        Default is False.
        NOTE that `use_test_transcripts` must be False if `do_k_fold_training` is True.
    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :param gnomad_data_type: gnomAD data type. Used to retrieve public release and coverage resources.
        Must be one of "exomes" or "genomes" (check is done within `public_release`).
        Default is "exomes".
    :param loftee_hc_str: String indicating that LOFTEE a loss-of-function variant is predcited to cause
    :return: None; writes amino acid Table to resource path.
    """
    if use_test_transcripts and do_k_fold_training:
        raise DataException("Cannot generate k-fold models with test transcripts!")

    logger.info("Reading in VEP context HT...")
    # NOTE: Keeping all variant types here because need synonymous and nonsense variants to calculate missense badness
    context_ht = process_context_ht(filter_to_missense=False, add_annotations=False)

    logger.info("Selecting relevant annotations...")
    context_ht = context_ht.select(
        transcript=context_ht.transcript_consequences.transcript_id,
        consequence_terms=context_ht.transcript_consequences.consequence_terms,
        most_severe_consequence=context_ht.transcript_consequences.most_severe_consequence,
        amino_acids=context_ht.transcript_consequences.amino_acids,
        codons=context_ht.transcript_consequences.codons,
        lof=context_ht.transcript_consequences.lof,
        lof_flags=context_ht.transcript_consequences.lof_flags,
    )

    logger.info(
        "Filtering non-coding rows and rows with uninformative/unknown codons..."
    )
    context_ht = annotate_and_filter_codons(context_ht)

    logger.info("Checkpointing HT before checking LoF information...")
    context_ht = context_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/codons.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    logger.info("Checking LOFTEE LoF information...")
    loftee_lof_agg = context_ht.aggregate(hl.agg.counter(context_ht.lof))
    logger.info("LOFTEE aggregation: %s", loftee_lof_agg)

    logger.info("Filtering to LOFTEE HC pLoF without flags...")
    context_ht = context_ht.filter(
        hl.is_missing(context_ht.lof) | (context_ht.lof == loftee_hc_str)
    )
    hc_lof_flag_agg = context_ht.aggregate(hl.agg.counter(context_ht.lof_flags))
    logger.info("Flag counter for LOFTEE HC pLoF: %s", hc_lof_flag_agg)
    context_ht = context_ht.filter(hl.is_missing(context_ht.lof_flags))

    logger.info("Checkpointing HT before joining with gnomAD data...")
    context_ht = context_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/codons_hc.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    logger.info("Filtering sites using gnomAD %s...", gnomad_data_type)
    context_ht = filter_context_using_gnomad(
        context_ht,
        gnomad_data_type,
    )

    logger.info("Adding observed annotation...")
    context_ht = add_obs_annotation(context_ht)

    logger.info("Checkpointing HT after joining with gnomAD data...")
    context_ht = context_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/codons_filt.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    logger.info(
        "Getting observed to expected ratio, rekeying Table, and writing to output path..."
    )
    # Note that `get_oe_annotation` is pulling the missense OE ratio
    context_ht = get_oe_annotation(context_ht, freeze)
    context_ht = context_ht.key_by("locus", "alleles", "transcript")
    context_ht = context_ht.select(
        "ref",
        "alt",
        "observed",
        "codons",
        "amino_acids",
        "oe",
    )

    logger.info("Checkpointing HT before filtering by transcript...")
    context_ht = context_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/codons_oe.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    if use_test_transcripts or not do_k_fold_training:
        logger.info(
            "Writing out HT for %s transcripts only...",
            "test" if use_test_transcripts else "training",
        )
        filter_transcripts = hl.experimental.read_expression(
            test_transcripts_path
            if use_test_transcripts
            else training_transcripts_path()
        )
        context_ht = filter_context_to_transcript_cds(context_ht, filter_transcripts)
        context_ht.write(
            amino_acids_oe_path(is_test=use_test_transcripts, freeze=freeze),
            overwrite=overwrite_output,
        )
    else:
        logger.info(
            "Writing out amino acid OE HT for each of the %i-fold training and validation sets...",
            fold_k,
        )
        for i in range(1, fold_k + 1):
            for is_val in [True, False]:
                logger.info("Writing out amino acid OE HT for fold %i...", i)
                filter_transcripts = hl.experimental.read_expression(
                    training_transcripts_path(fold=i, is_val=is_val)
                )
                context_ht = filter_context_to_transcript_cds(
                    context_ht, filter_transcripts
                )
                context_ht.write(
                    amino_acids_oe_path(
                        is_test=use_test_transcripts,
                        fold=i,
                        is_val=is_val,
                        freeze=freeze,
                    ),
                    overwrite=overwrite_output,
                )


def variant_csq_expr(
    ref_expr: hl.expr.StringExpression, alt_expr: hl.expr.StringExpression
) -> hl.expr.StringExpression:
    """
    Determine variant consequence using reference and alternate amino acid annotations.

    Variant consequences are consistent with consequences kept in original missense badness work.
    TODO: Update variant consequences?

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


def aggregate_aa_and_filter_oe(
    ht: hl.Table,
    keep_high_oe: bool,
    oe_threshold: float = 0.6,
) -> hl.Table:
    """
    Split Table with all possible amino acid substitutions based on missense observed to expected (OE) ratio cutoff.

    Also group Table by reference and alternate amino acid, aggregate total observed and possible counts,
    and add mutation type annotation.

    :param hl.Table ht: Input Table with amino acid substitutions.
    :param bool keep_high_oe: Whether to filter to high missense OE values.
        If True, returns "boring" HT.
        If False, gets "bad" (low missense OE) Table.
    :param float oe_threshold: OE Threshold used to split Table.
        Rows with OE less than or equal to this threshold will be filtered if `keep_high_oe` is True, and
        rows with OE greater than this threshold will be kept.
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


def calculate_misbad(
    use_exac_oe_cutoffs: bool,
    overwrite_temp: bool,
    overwrite_output: bool,
    oe_threshold: float = 0.6,
    freeze: int = CURRENT_FREEZE,
) -> None:
    """
    Calculate missense badness score using Table with all amino acid substitutions and their missense observed/expected (OE) ratio.

    If `use_exac_oe_cutoffs` is set, will remove all rows with 0.6 < OE <= 0.8.

    .. note::
        Assumes table containing all possible amino acid substitutions and their missense OE ratio exists.

    :param bool use_exac_oe_cutoffs: Whether to use the same missense OE cutoffs as in ExAC missense badness calculation.
    :param bool overwrite_temp: Whether to overwrite intermediate temporary data if it already exists.
        If False, will read existing intermediate temporary data rather than overwriting.
    :param bool overwrite_output: Whether to entirely overwrite final output data if it already exists.
        If False, will read and modify existing output data by adding or modifying columns rather than overwriting entirely.
        If True, will clear existing output data and write new output data.
        The output Table is the missense badness score Table.
    :param float oe_threshold: OE Threshold used to split Table.
        Rows with OE less or equal to this threshold will be considered "low" OE, and
        rows with OE greater than this threshold will considered "high" OE.
        Default is 0.6.
    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: None; writes Table with missense badness score to resource path.
    """
    if not file_exists(amino_acids_oe.versions[freeze].path):
        raise DataException(
            "Table with all amino acid substitutions and missense OE doesn't exist!"
        )

    ht = amino_acids_oe.versions[freeze].ht()

    if use_exac_oe_cutoffs:
        logger.info("Removing rows with OE greater than 0.6 and less than 0.8...")
        ht = ht.filter((ht.oe <= 0.6) | (ht.oe > 0.8))

    logger.info(
        "Splitting input Table by OE to get synonymous and nonsense rates for high and low OE groups..."
    )
    logger.info("Creating high missense OE (OE > %s) HT...", oe_threshold)
    high_ht = aggregate_aa_and_filter_oe(ht, keep_high_oe=True)
    high_ht = high_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/amino_acids_high_oe.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    logger.info("Creating low missense OE (OE <= %s) HT...", oe_threshold)
    low_ht = aggregate_aa_and_filter_oe(ht, keep_high_oe=False)
    low_ht = low_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/amino_acids_low_oe.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

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
    ht = ht.transmute(mut_type=hl.coalesce(ht.mut_type, ht.mut_type_1))
    mb_ht = ht.group_by("ref", "alt").aggregate(
        high_low=(
            (hl.agg.sum(ht.high_obs) / hl.agg.sum(ht.high_pos))
            / (hl.agg.sum(ht.low_obs) / hl.agg.sum(ht.low_pos))
        )
    )
    mb_ht = mb_ht.annotate(mut_type=variant_csq_expr(mb_ht.ref, mb_ht.alt))

    logger.info("Calculating synonymous rates...")
    syn_obs_high = get_total_csq_count(high_ht, csq="syn", count_field="high_obs")
    syn_pos_high = get_total_csq_count(high_ht, csq="syn", count_field="high_pos")
    syn_obs_low = get_total_csq_count(low_ht, csq="syn", count_field="low_obs")
    syn_pos_low = get_total_csq_count(low_ht, csq="syn", count_field="low_pos")
    syn_rate = (syn_obs_high / syn_pos_high) / (syn_obs_low / syn_pos_low)
    logger.info("Synonymous rate: %f", syn_rate)

    logger.info("Calculating nonsense rates...")
    non_obs_high = get_total_csq_count(high_ht, csq="non", count_field="high_obs")
    non_pos_high = get_total_csq_count(high_ht, csq="non", count_field="high_pos")
    non_obs_low = get_total_csq_count(low_ht, csq="non", count_field="low_obs")
    non_pos_low = get_total_csq_count(low_ht, csq="non", count_field="low_pos")
    non_rate = (non_obs_high / non_pos_high) / (non_obs_low / non_pos_low)
    logger.info("Nonsense rate: %f", non_rate)

    logger.info("Calculating missense badness...")
    mb_ht = mb_ht.annotate(
        misbad=hl.or_missing(
            mb_ht.mut_type == "mis",
            # Cap missense badness at 1
            hl.min((mb_ht.high_low - syn_rate) / (non_rate - syn_rate), 1),
        ),
    )

    mb_ht.write(misbad.versions[freeze].path, overwrite=overwrite_output)
    logger.info("Output missense badness HT fields: %s", set(mb_ht.row))
