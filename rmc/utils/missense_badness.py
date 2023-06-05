import logging

import hail as hl
from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists

from rmc.resources.basics import TEMP_PATH_WITH_FAST_DEL
from rmc.resources.reference_data import (
    FOLD_K,
    train_val_test_transcripts_path,
)
from rmc.resources.rmc import CURRENT_FREEZE, amino_acids_oe_path, misbad_path
from rmc.utils.constraint import add_obs_annotation, get_oe_annotation
from rmc.utils.generic import (
    annotate_and_filter_codons,
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
    do_k_fold_training: bool = False,
    freeze: int = CURRENT_FREEZE,
    gnomad_data_type: str = "exomes",
    loftee_hc_str: str = "HC",
) -> None:
    """
    Prepare Table(s) with all possible amino acid substitutions and their missense observed to expected (OE) ratio.

    Steps:
        - Import VEP context Table and filter to keep every possible amino acid substitution
        (every codon > codon change).
        - Filter Table to rows that aren't present in gnomAD or are rare in gnomAD (using `keep_criteria`).
        - Add observed and OE annotation.
        - Filter to specified transcript set(s).
        - Write to `amino_acids_oe_path` resource path(s).

    :param overwrite_temp: Whether to overwrite intermediate temporary (OE-independent) data if it already exists.
        If False, will read existing intermediate temporary data rather than overwriting.
    :param do_k_fold_training: Whether to generate k-fold models with the training transcripts.
        If False, will use all training transcripts in calculation of a single model.
        If True, will calculate k models corresponding to the k-folds of the training transcripts.
        Default is False.
    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :param gnomad_data_type: gnomAD data type. Used to retrieve public release and coverage resources.
        Must be one of "exomes" or "genomes" (check is done within `public_release`).
        Default is "exomes".
    :param loftee_hc_str: String indicating that a variant is predicted to cause loss-of-function with high confidence by LOFTEE.
    :return: None; writes amino acid Table(s) to resource path(s).
    """
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
    loftee_lof_agg_path = f"{TEMP_PATH_WITH_FAST_DEL}/lof_agg.he"
    if not overwrite_temp and file_exists(loftee_lof_agg_path):
        loftee_lof_agg = hl.eval(hl.experimental.read_expression(loftee_lof_agg_path))
    else:
        loftee_lof_agg = context_ht.aggregate(hl.agg.counter(context_ht.lof))
        hl.experimental.write_expression(
            hl.literal(loftee_lof_agg), loftee_lof_agg_path
        )
    logger.info("LOFTEE aggregation: %s", loftee_lof_agg)

    logger.info("Filtering to LOFTEE HC pLoF without flags...")
    context_ht = context_ht.filter(
        hl.is_missing(context_ht.lof) | (context_ht.lof == loftee_hc_str)
    )
    hc_lof_flag_agg_path = f"{TEMP_PATH_WITH_FAST_DEL}/lof_flag_agg.he"
    if not overwrite_temp and file_exists(hc_lof_flag_agg_path):
        hc_lof_flag_agg = hl.eval(hl.experimental.read_expression(hc_lof_flag_agg_path))
    else:
        hc_lof_flag_agg = context_ht.aggregate(hl.agg.counter(context_ht.lof_flags))
        hl.experimental.write_expression(
            hl.literal(hc_lof_flag_agg), hc_lof_flag_agg_path
        )
    logger.info("Flag counter for LOFTEE HC pLoF: %s", hc_lof_flag_agg)
    context_ht = context_ht.filter(hl.is_missing(context_ht.lof_flags))

    logger.info("Checkpointing HT before joining with gnomAD data...")
    context_ht = context_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/codons_hc.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    logger.info("Filtering sites using gnomAD %s...", gnomad_data_type)
    context_ht = filter_context_using_gnomad(context_ht, gnomad_data_type)

    logger.info("Adding observed annotation...")
    context_ht = add_obs_annotation(context_ht)

    logger.info("Checkpointing HT after joining with gnomAD data...")
    context_ht = context_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/codons_filt.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    logger.info("Getting observed to expected ratio and rekeying Table...")
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

    def _filter_transcripts(
        ht: hl.Table,
        fold: int = None,
    ) -> None:
        """
        Filter amino acid OE HT to a set of training transcripts and write out.

        :param hl.Table ht: Input table.
        :param int fold: Fold number in training set to select training transcripts from.
            If not None, the Table is generated from variants in only training transcripts
            from the specified fold. If None, the Table is generated from variants in all training transcripts.
            Default is None.
        :return: None; function writes HT to specified path.
        """
        filter_transcripts = hl.experimental.read_expression(
            train_val_test_transcripts_path(fold=fold)
        )
        ht = ht.filter(filter_transcripts.contains(ht.transcript))
        ht.write(amino_acids_oe_path(fold=fold, freeze=freeze), overwrite=True)

    if not do_k_fold_training:
        logger.info("Writing out HT for training transcripts only...")
        _filter_transcripts(ht=context_ht)
    else:
        logger.info(
            "Writing out amino acid OE HT for each of the %i-fold training sets...",
            FOLD_K,
        )
        for i in range(1, FOLD_K + 1):
            logger.info("Writing out amino acid OE HT for fold %i...", i)
            _filter_transcripts(ht=context_ht, fold=i)


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
    do_k_fold_training: bool = False,
    oe_threshold: float = 0.6,
    freeze: int = CURRENT_FREEZE,
) -> None:
    """
    Calculate missense badness scores using Table with all amino acid substitutions and their missense observed/expected (OE) ratio.

    If `use_exac_oe_cutoffs` is set, will remove all rows with 0.6 < OE <= 0.8.

    .. note::
        Assumes table(s) containing all possible amino acid substitutions and their missense OE ratio exists.

    :param bool use_exac_oe_cutoffs: Whether to use the same missense OE cutoffs as in ExAC missense badness calculation.
    :param bool overwrite_temp: Whether to overwrite intermediate temporary data if it already exists.
        If False, will read existing intermediate temporary data rather than overwriting.
    :param bool do_k_fold_training: Whether to generate k-fold models with the training transcripts.
        If False, will use all training transcripts in calculation of a single model.
        If True, will calculate k models corresponding to the k-folds of the training transcripts.
        Default is False.
    :param float oe_threshold: OE Threshold used to split Table.
        Rows with OE less or equal to this threshold will be considered "low" OE, and
        rows with OE greater than this threshold will considered "high" OE.
        Default is 0.6.
    :param int freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: None; writes Table(s) with missense badness scores to resource path(s).
    """
    if not do_k_fold_training:
        if not file_exists(amino_acids_oe_path(freeze=freeze)):
            raise DataException(
                "Table with all amino acid substitutions and missense OE in training transcripts doesn't exist!"
            )
    else:
        if not all(
            [
                file_exists(amino_acids_oe_path(fold=i, freeze=freeze))
                for i in range(1, FOLD_K + 1)
            ]
        ):
            raise DataException(
                "Not all k-fold training tables with all amino acid substitutions and missense OE exist!"
            )

    def _create_misbad_model(
        temp_label: str,
        fold: int = None,
    ) -> None:
        """
        Calculate missense badness scores from a Table of amino acid substitutions and their missense OE ratios.

        :param str temp_label: Suffix to add to temporary data paths to avoid conflicting names for different models.
        :param int fold: Fold number in training set to select training transcripts from.
            If not None, the Table is generated from variants in only training transcripts
            from the specified fold. If None, the Table is generated from variants in all training transcripts.
            Default is None.
        :return: None; writes Table with missense badness scores to resource path.
        """
        ht = hl.read_table(amino_acids_oe_path(fold=fold, freeze=freeze))
        if use_exac_oe_cutoffs:
            logger.info("Removing rows with OE greater than 0.6 and less than 0.8...")
            ht = ht.filter((ht.oe <= 0.6) | (ht.oe > 0.8))

        logger.info(
            "Splitting input Table by OE to get synonymous and nonsense rates for high and low OE groups..."
        )
        oe_hts = {}
        for keep_high_oe in [True, False]:
            oe_direction = "high" if keep_high_oe else "low"
            logger.info(
                "Creating %s missense OE (OE > %s) HT...", oe_direction, oe_threshold
            )
            oe_hts[oe_direction] = aggregate_aa_and_filter_oe(
                ht, keep_high_oe=keep_high_oe
            )
            oe_hts[oe_direction] = oe_hts[oe_direction].transmute(
                **{
                    f"{oe_direction}_obs": oe_hts[oe_direction].obs,
                    f"{oe_direction}_pos": oe_hts[oe_direction].possible,
                }
            )
            oe_hts[oe_direction] = oe_hts[oe_direction].checkpoint(
                f"{TEMP_PATH_WITH_FAST_DEL}/amino_acids_{oe_direction}_oe{temp_label}.ht",
                _read_if_exists=not overwrite_temp,
                overwrite=overwrite_temp,
            )

        logger.info("Re-joining split HTs to calculate missense badness...")
        mb_ht = oe_hts["high"].join(oe_hts["low"], how="outer")
        mb_ht = mb_ht.transmute(mut_type=hl.coalesce(ht.mut_type, ht.mut_type_1))
        mb_ht = mb_ht.annotate(
            high_low=(
                (mb_ht.high_obs / mb_ht.high_pos) / (mb_ht.low_obs / mb_ht.low_pos)
            )
        )

        total_rates = {}
        for mut_type in ["syn", "non"]:
            logger.info("Calculating %s rates...", mut_type)
            total_rates[mut_type] = (
                get_total_csq_count(
                    oe_hts["high"], csq=mut_type, count_field="high_obs"
                )
                / get_total_csq_count(
                    oe_hts["high"], csq=mut_type, count_field="high_pos"
                )
            ) / (
                get_total_csq_count(oe_hts["low"], csq=mut_type, count_field="low_obs")
                / get_total_csq_count(
                    oe_hts["low"], csq=mut_type, count_field="low_pos"
                )
            )
            logger.info("%s rate: %f", mut_type, total_rates[mut_type])

        logger.info("Calculating missense badness...")
        mb_ht = mb_ht.annotate(
            misbad=hl.or_missing(
                mb_ht.mut_type == "mis",
                # Cap missense badness at 1
                hl.min(
                    (mb_ht.high_low - total_rates["syn"])
                    / (total_rates["non"] - total_rates["syn"]),
                    1,
                ),
            ),
        )

        mb_ht.write(misbad_path(fold=fold, freeze=freeze), overwrite=True)

    if not do_k_fold_training:
        _create_misbad_model(
            temp_label="_train",
        )
    else:
        for i in range(1, FOLD_K + 1):
            _create_misbad_model(temp_label=f"_train_fold{i}", fold=i)
