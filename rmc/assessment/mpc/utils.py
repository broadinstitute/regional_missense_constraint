import logging
from typing import Dict, List, Union

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists

from rmc.resources.basics import context_with_oe
from rmc.resources.grch37.reference_data import (
    cadd,
    dbnsfp,
    DE_NOVO_COUNTS,
    get_mpc_case_control_ht_path,
    get_mpc_pct_comparison_path,
    mpc_comparison,
)


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("mpc_assessment_utils")
logger.setLevel(logging.INFO)


def prep_mpc_histogram_tsv(
    case_ht: hl.Table,
    control_ht: hl.Table,
    output_tsv_path: str,
    keep_asd: bool,
    keep_dd: bool,
    case_control_field: str = "case_control",
    asd_str: str = "ASD",
    dd_str: str = "DD",
) -> None:
    """
    Create TSV of MPC scores of de novo variants from neurodevelopmental disorders (NDD) cases and controls.

    Used to create stacked histograms of MPC scores from cases and controls.

    .. note::
        Input case and control HT must be annotated with MPC and `case_control_field`.

    :param hl.Table case_ht: Table with de novo variants from NDD cases.
    :param hl.Table control_ht: Table with de novo variants from NDD controls.
    :param str output_tsv_path: Where to store output TSV.
    :param bool keep_asd: Whether to keep variants from cases with Autism Spectrum Disorder (ASD).
    :param bool keep_dd: Whether to keep variants from cases with developmental disorders (DD).
    :param str case_control_field: Field describing whether variant is from a case or control.
        Default is 'case_control'.
    :param str asd_str: String describing whether case has ASD. Default is 'ASD'.
    :param str dd_str: String describing whether case has DD. Default is 'DD'.
    :return: None; function writes TSV to specified path.
    """
    assert any(
        [keep_asd, keep_dd]
    ), "At least one of `keep_asd` and `keep_dd` must be True!"

    if not keep_asd:
        logger.info("Removing ASD cases...")
        case_ht = case_ht.filter(case_ht[case_control_field] != asd_str)
    if not keep_dd:
        logger.info("Removing DD cases...")
        case_ht = case_ht.filter(case_ht[case_control_field] != dd_str)

    logger.info("Joining case and control HTs...")
    ht = (
        case_ht.key_by("case_control")
        .select("mpc")
        .union(control_ht.key_by("case_control").select("mpc"))
    )

    logger.info("Checkpointing Table and exporting to TSV...")
    if keep_asd and keep_dd:
        asd_only = dd_only = False
    else:
        asd_only = keep_asd
        dd_only = keep_dd
    ht = ht.checkpoint(get_mpc_case_control_ht_path(asd_only, dd_only), overwrite=True)
    ht.export(output_tsv_path)


def prep_rate_ratio_tsv(
    output_tsv_path: str,
    keep_asd: bool,
    keep_dd: bool,
) -> None:
    """
    Create TSV of total number of variants from NDD cases and controls per MPC bin.

    Used as input to two-sided Poisson exact test to calculate rate ratios.

    :param str output_tsv_path: Where to store output TSV.
    :param bool keep_asd: Whether to keep variants from cases with Autism Spectrum Disorder (ASD).
    :param bool keep_dd: Whether to keep variants from cases with developmental disorders (DD).
    :return: None; function writes TSV to specified path.
    """
    logger.info("Reading in Table...")
    if keep_asd and keep_dd:
        asd_only = dd_only = False
    else:
        asd_only = keep_asd
        dd_only = keep_dd
    ht_path = get_mpc_case_control_ht_path(asd_only, dd_only)

    if not file_exists(ht_path):
        raise DataException("%s does not exist! Double check and rerun.")
    ht = hl.read_table(get_mpc_case_control_ht_path(asd_only, dd_only))

    logger.info("Getting total number of cases and controls...")
    total_cases = 0
    if keep_asd:
        total_cases += DE_NOVO_COUNTS["asd_only"]
    if keep_dd:
        total_cases += DE_NOVO_COUNTS["dd_only"]
    total_controls = DE_NOVO_COUNTS["controls"]

    logger.info("Annotating with MPC bin, grouping HT, and exporting...")
    ht = ht.annotate(
        mpc_bin=hl.case()
        .when((ht.mpc >= 1) & (ht.mpc < 2), 1)
        .when((ht.mpc >= 2) & (ht.mpc < 3), 2)
        .when(ht.mpc >= 3, 3)
        .default(0)
    )
    ht = ht.group_by("mpc_bin").aggregate(count=hl.agg.count())
    ht = ht.group_by("mpc_bin").aggregate(
        n_case=hl.agg.filter(ht.case_control != "control", hl.agg.count()),
        n_control=hl.agg.filter(ht.case_control == "control", hl.agg.count()),
        total_cases=total_cases,
        total_controls=total_controls,
    )
    ht = ht.annotate(
        rate_per_case=ht.n_case / ht.total_controls,
        rate_per_control=ht.n_control / ht.total_controls,
    )
    ht.export(output_tsv_path)


def prep_mpc_comparison_ht(
    case_ht: hl.Table,
    control_ht: hl.Table,
    temp_path_with_del: str = "gs://gnomad-tmp/mpc",
) -> None:
    """
    Create Table of de novo variants from NDD cases and controls to compare MPC performance.

    Annotations added are: Polyphen-2, SIFT, CADD, and REVEL scores.

    :param hl.Table case_ht: Table with de novo variants from NDD cases.
    :param hl.Table control_ht: Table with de novo variants from NDD controls.
    :param str temp_path_with_del: Path to bucket to store temporary data with automatic deletion policy.
        Default is 'gs://gnomad-tmp/mpc'.
        TODO: Update this to `temp_path` (and set automatic deletion policy.)
    :return: None; function writes Table to resource (`mpc_comparison`) path.
    """
    logger.info("Joining case and control HT...")
    ht = case_ht.union(control_ht)
    ht = ht.checkpoint(f"{temp_path_with_del}/all_ndd_mpc.ht", overwrite=True)

    logger.info("Adding annotations and checkpointing after each join...")
    # Read in `context_with_oe to get` Polyphen-2 and SIFT
    context_ht = context_with_oe.ht()
    context_ht = context_ht.key_by("locus", "alleles", "transcript")
    ht = ht.annotate(
        polyphen=context_ht[ht.locus, ht.alleles, ht.transcript].polyphen.score,
        sift=context_ht[ht.locus, ht.alleles, ht.transcript].sift.score,
    )
    ht = ht.checkpoint(f"{temp_path_with_del}/all_ndd_mpc_pph_sift.ht", overwrite=True)

    # CADD
    cadd_ht = cadd.ht().select("PHRED")
    ht = ht.annotate(cadd=cadd_ht[ht.key].PHRED)
    ht = ht.checkpoint(
        f"{temp_path_with_del}/all_ndd_mpc_pph_sift_cadd.ht", overwrite=True
    )

    # Read in dbNSFP Table to get REVEL score
    revel_ht = dbnsfp.ht().select("REVEL_score")
    # Split multiallelics in REVEL HT and convert score to float
    # (score is currently a string)
    revel_ht = hl.split_multi_hts(revel_ht)
    revel_ht = revel_ht.select(revel=hl.float(revel_ht.REVEL_score))
    revel_ht = revel_ht.collect_by_key()
    revel_ht = revel_ht.transmute(revel=hl.nanmax(revel_ht.values.revel))
    ht = ht.annotate(revel=revel_ht[ht.key].revel)
    ht = ht.checkpoint(
        f"{temp_path_with_del}/all_ndd_mpc_pph_sift_cadd_revel.ht", overwrite=True
    )

    logger.info("Filtering to defined annotations...")
    ht = ht.filter(
        hl.is_defined(ht.mpc)
        & hl.is_defined(ht.polyphen)
        & hl.is_defined(ht.sift)
        & hl.is_defined(ht.cadd)
        & hl.is_defined(ht.revel)
    )
    ht.write(mpc_comparison.path, overwrite=True)


def compare_frac_of_top_x_var(
    score: str,
    top_x_percent: float,
    control_str: str = "control",
    temp_path_with_del: str = "gs://gnomad-tmp/mpc",
) -> Dict[str, Union[str, float]]:
    """
    Pull the top x% of variants with highest value for desired score, determine fraction of variants from cases, and determine odds ratio.

    Used to compare MPC performance to other metrics. Odds ratio is calculated using Fisher's exact test.

    .. note::
        - MPC comparison HT must be annotated with score, and score must match the field name
            (e.g., if score is "revel", then MPC comparison HT must have the annotation "revel")

    :param str score: Score of interest.
    :param float top_x_percent: Desired percent value. E.g., top_x_percent = 5 means this function will compare
        the fraction of variants with the top 5% largest score values from cases vs controls.
    :param str control_str: String noting whether variant is from a control sample. Default is 'control'.
    :param str temp_path_with_del: Path to bucket to store temporary data with automatic deletion policy.
        Default is 'gs://gnomad-tmp/mpc'.
        TODO: Update this to `temp_path` (and set automatic deletion policy.)
    """
    ht = mpc_comparison.ht()
    total_count = ht.count()
    top_row_count = int(top_x_percent * total_count)
    bottom_row_count = total_count - top_row_count

    # Order Table by score (descending)
    ht = ht.key_by()
    ht = ht.order_by(hl.desc(ht[score]))
    ht = ht.checkpoint(f"{temp_path_with_del}/score_temp.ht", overwrite=True)

    # Split the table into the top x% and everything else (called top and bottom here)
    ht_top = ht.head(top_row_count)
    ht_top = ht_top.checkpoint(f"{temp_path_with_del}/score_top.ht", overwrite=True)
    ht_bottom = ht.tail(bottom_row_count)
    ht_bottom = ht_bottom.checkpoint(
        f"{temp_path_with_del}/score_bottom.ht", overwrite=True
    )

    # Get fraction of top x% variants that come from cases
    frac_cases = ht_top.aggregate(hl.agg.fraction(ht_top.case_control != control_str))

    # Get total counts of cases and control from both top and bottom HTs to run Fisher's exact test
    ht_n_case = ht_top.aggregate(hl.agg.count_where(ht_top.case_control != control_str))
    ht_n_control = ht_top.aggregate(
        hl.agg.count_where(ht_top.case_control == control_str)
    )
    htb_n_case = ht_bottom.aggregate(
        hl.agg.count_where(ht_bottom.case_control != control_str)
    )
    htb_n_control = ht_bottom.aggregate(
        hl.agg.count_where(ht_bottom.case_control == control_str)
    )

    # Run Fisher's exact test
    fisher_struct = hl.eval(
        hl.fisher_exact_test(ht_n_case, htb_n_case, ht_n_control, htb_n_control)
    )
    return {
        "score": score,
        "frac_cases": frac_cases,
        "odds_ratio": fisher_struct.odds_ratio,
        "p_value": fisher_struct.p_value,
    }


def compare_mpc_using_top_x_var(
    top_x_percent: float,
    scores: List[str] = ["revel", "cadd", "polyphen"],
    schema: str = "array<struct{score: str, frac_cases: float, odds_ratio: float, p_value: float}>",
) -> None:
    """
    Compare MPC score to other metrics by comparing the fraction of variants from cases in the top x% of each score.

    Function compares fraction of variants from NDD cases vs controls in two categories:
        - Variants with >= top x percent of each score
        - Variants < top x percent of each score

    .. note::
        `mpc_comparison` Table must be annotated with each of the scores provided in `scores`.

    :param float top_x_percent: Desired percent value. E.g., top_x_percent = 5 means this function will compare
        the fraction of variants with the top 5% largest score values from cases vs controls.
    :param List[str] scores: Scores to compare against MPC. Default is ['revel', 'cadd', 'polyphen'].
    :param str schema: Schema used to create comparison Table from list of dictionaries returned by
        `compare_frac_of_top_x_var`. Note that the dictionary keys must match the field names provided in this schema.
        Default is 'array<struct{score: str, frac_cases: float, odds_ratio: float, p_value: float}>'.
    :return: None; function writes Table to resource path (`get_mpc_pct_comparison_path`).
    """
    score_list = []
    ht = mpc_comparison.ht()
    for score in scores:
        score_list.append(compare_frac_of_top_x_var(ht, score, top_x_percent))

    ht = hl.Table.parallelize(hl.literal(score_list, schema))
    ht = ht.checkpoint(get_mpc_pct_comparison_path(top_x_percent), overwrite=True)
    ht.show()
