"""This script contains functions used to search for two simultaneous breaks."""
import logging
from typing import List

import hail as hl

from gnomad.utils.file_utils import file_exists, parallel_file_exists

from rmc.resources.basics import SIMUL_BREAK_TEMP_PATH, TEMP_PATH_WITH_FAST_DEL
from rmc.resources.rmc import (
    CHISQ_THRESHOLDS,
    MIN_CHISQ_THRESHOLD,
    MIN_EXP_MIS,
    simul_search_round_bucket_path,
    simul_sections_split_by_len_path,
)
from rmc.utils.constraint import (
    get_dpois_expr,
    get_obs_exp_expr,
    get_max_chisq_per_group,
)


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("search_for_two_breaks_utils")
logger.setLevel(logging.INFO)


def group_no_single_break_found_ht(
    ht_path: str,
    out_ht_path: str,
    group_str: str,
    overwrite: bool,
) -> None:
    """
    Group HT containing transcripts/transcript sections that don't have a single significant break and collect annotations into lists.

    Group HT by `group_str` and collect positions, cumulative obs, cumulative exp into lists.

    This creates the input to the two simultaneous breaks search (`search_for_two_breaks`).

    .. note::
        - Expects that input Table is keyed by locus and `group_str`.
            This is *required*, as the function expects that loci are sorted in the input Table.
        - Expects that input Table is annotated with OE for transcript/transcript section and
            this annotation is named `section_oe`.

    :param ht_path: Path to input Table with transcripts/transcript sections that didn't
        have a single significant break.
    :param out_ht_path: Path to output Table grouped by transcripts/transcript sections with lists
        of cumulative observed, expected missense counts.
    :param group_str: Field used to group observed and expected values.
    :param overwrite: Whether to overwrite existing grouped Table.
    :return: None; writes Table grouped by transcript/transcript section with cumulative observed, expected missense counts
        and all positions collected into lists to resource path.
    """
    # Aggregating values into a struct here to force the positions and observed, expected missense values to stay sorted
    # `hl.agg.collect` does not guarantee order: https://hail.is/docs/0.2/aggregators.html#hail.expr.aggregators.collect
    ht = hl.read_table(ht_path)
    group_ht = ht.group_by(group_str).aggregate(
        values=hl.sorted(
            hl.agg.collect(
                hl.struct(
                    locus=ht.locus,
                    cum_exp=ht.cumulative_exp,
                    cum_obs=ht.cumulative_obs,
                    positions=ht.locus.position,
                ),
            ),
            key=lambda x: x.locus,
        ),
        total_oe=hl.agg.take(ht.section_oe, 1)[0],
    )
    group_ht = group_ht.annotate(max_idx=hl.len(group_ht.values.positions) - 1)
    group_ht = group_ht.transmute(
        cum_obs=group_ht.values.cum_obs,
        cum_exp=group_ht.values.cum_exp,
        positions=group_ht.values.positions,
    )
    group_ht.write(out_ht_path, overwrite=overwrite)


def split_sections_by_len(
    ht_path: str,
    group_str: str,
    search_num: int,
    missense_len_threshold: int,
    overwrite: bool,
    freeze: int = CURRENT_FREEZE,
) -> None:
    """
    Split transcripts/transcript sections based on the specified number of possible missense sites.

    This is necessary because transcripts/transcript sections with more possible missense sites
    take longer to run through `hl.experimental.loop`.

    :param ht: Path to input Table (Table written using `group_no_single_break_found_ht`).
    :param group_str: Field used to group observed and expected values.
    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param missense_len_threshold: Cutoff based on possible number of missense sites in section.
    :param overwrite: Whether to overwrite existing SetExpressions.
    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: None; writes SetExpressions to resource paths.
    """
    ht = hl.read_table(ht_path)
    logger.info("Annotating HT with length of cumulative observed list annotation...")
    # This length is the number of positions with possible missense sites that need to be searched
    # Not using transcript/transcript section size here because transcript/transcript section size
    # doesn't necessarily reflect the number of positions that need to be searched
    ht = ht.annotate(missense_list_len=ht.max_idx + 1)

    logger.info(
        "Splitting sections into two categories: list length < %i and list length >= %i...",
        missense_len_threshold,
        missense_len_threshold,
    )
    under_threshold = ht.aggregate(
        hl.agg.filter(
            ht.missense_list_len < missense_len_threshold,
            hl.agg.collect_as_set(ht[group_str]),
        )
    )
    over_threshold = ht.aggregate(
        hl.agg.filter(
            ht.missense_list_len >= missense_len_threshold,
            hl.agg.collect_as_set(ht[group_str]),
        )
    )
    logger.info(
        "Found %i sections under threshold and %i sections over threshold",
        len(under_threshold),
        len(over_threshold),
    )

    if len(under_threshold) > 0:
        hl.experimental.write_expression(
            under_threshold,
            simul_sections_split_by_len_path(
                search_num=search_num,
                is_over_threshold=False,
                freeze=freeze,
            ),
            overwrite,
        )
    if len(over_threshold) > 0:
        hl.experimental.write_expression(
            over_threshold,
            simul_sections_split_by_len_path(
                search_num=search_num,
                is_over_threshold=True,
                freeze=freeze,
            ),
            overwrite,
        )


def get_sections_to_run(
    sections: List[str],
    search_num: int,
    in_parallel: bool = True,
    freeze: int = CURRENT_FREEZE,
) -> List[str]:
    """
    Check if any transcripts/sections have been previously searched by searching for success TSV existence.

    .. note::
        This step needs to be run locally due to permissions involved with `parallel_file_exists`.

    :param List[str] sections: List of transcripts/transcript sections to check.
    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param bool in_parallel: Whether to check if successful file exist in parallel.
        If True, must be run locally and not in Dataproc. Default is True.
    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: List of transcripts/sections that didn't have success TSVs and therefore still need to be processed.
    """
    logger.info("Checking if any transcripts have already been searched...")
    success_file_path = simul_search_round_bucket_path(
        search_num=search_num,
        bucket_type="success_files",
        freeze=freeze,
    )
    section_success_map = {}
    sections_to_run = []
    for section in sections:
        section_success_map[section] = f"{success_file_path}/{section}.tsv"

    if in_parallel:
        # Use parallel_file_exists if in_parallel is set to True
        success_tsvs_exist = parallel_file_exists(list(section_success_map.values()))
        for section in sections:
            if not success_tsvs_exist[section_success_map[section]]:
                sections_to_run.append(section)
    else:
        # Otherwise, use `file_exists`
        for section in sections:
            if not file_exists(section_success_map[section]):
                sections_to_run.append(section)

    return sections_to_run


def calculate_window_chisq(
    max_idx: hl.expr.Int32Expression,
    i: hl.expr.Int32Expression,
    j: hl.expr.Int32Expression,
    cum_obs: hl.expr.ArrayExpression,
    cum_exp: hl.expr.ArrayExpression,
    total_oe: hl.expr.Float64Expression,
    min_num_exp_mis: hl.expr.Float64Expression,
) -> hl.expr.Float64Expression:
    """
    Calculate chi square significance value for each possible simultaneous breaks window.

    Used only when calculating simultaneous breaks.

    Chi square formula: 2 * (hl.log(total_alt) - hl.log(total_null))

    :param hl.expr.Int32Expression max_idx: Largest list index value.
    :param hl.expr.Int32Expression i: Smaller list index value corresponding to the smaller position of the two break window.
    :param hl.expr.Int32Expression j: Larger list index value corresponding to the larger position of the two break window.
    :param hl.expr.ArrayExpression cum_obs: List containing cumulative observed missense values.
    :param hl.expr.ArrayExpression cum_exp: List containing cumulative expected missense values.
    :param hl.expr.Float64Expression total_oe: Transcript/transcript section overall observed/expected (OE) missense ratio.
    :param hl.expr.Float64Expression min_num_exp_mis: Minimum expected missense value for all three windows defined by two possible
        simultaneous breaks.
    :return: Chi square significance value.
    """
    return hl.or_missing(
        # Return missing when breakpoints do not create 3 windows
        # (return missing when index is the smallest or largest position)
        (i != 0) & (j != max_idx) &
        # Return missing when any of the windows do not have the minimum
        # number of expected missense variants
        (cum_exp[i - 1] >= min_num_exp_mis)
        & ((cum_exp[j] - cum_exp[i - 1]) >= min_num_exp_mis)
        & ((cum_exp[-1] - cum_exp[j]) >= min_num_exp_mis),
        (
            2
            * (
                # Create alt distribution
                (
                    # Create alt distribution for section [start_pos, pos[i])
                    # The missense values for this section are the cumulative values at
                    # one index smaller than index i
                    get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=get_obs_exp_expr(
                            True,
                            cum_obs[i - 1],
                            hl.max(cum_exp[i - 1], 1e-09),
                        ),
                        obs_expr=cum_obs[i - 1],
                        exp_expr=hl.max(cum_exp[i - 1], 1e-09),
                    )
                    # Create alt distribution for section [pos[i], pos[j]]
                    # The missense values for this section are the cumulative values at index j
                    # minus the cumulative values at index i -1
                    + get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=get_obs_exp_expr(
                            True,
                            (cum_obs[j] - cum_obs[i - 1]),
                            hl.max(cum_exp[j] - cum_exp[i - 1], 1e-09),
                        ),
                        obs_expr=cum_obs[j] - cum_obs[i - 1],
                        exp_expr=hl.max(cum_exp[j] - cum_exp[i - 1], 1e-09),
                    )
                    # Create alt distribution for section (pos[j], end_pos]
                    # The missense values for this section are the cumulative values at the last index
                    # minus the cumulative values at index j
                    + get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=get_obs_exp_expr(
                            True,
                            (cum_obs[-1] - cum_obs[j]),
                            hl.max(cum_exp[-1] - cum_exp[j], 1e-09),
                        ),
                        obs_expr=cum_obs[-1] - cum_obs[j],
                        exp_expr=hl.max(cum_exp[-1] - cum_exp[j], 1e-09),
                    )
                )
                # Create null distribution
                - (
                    get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=total_oe,
                        obs_expr=cum_obs[i - 1],
                        exp_expr=hl.max(cum_exp[i - 1], 1e-09),
                    )
                    + get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=total_oe,
                        obs_expr=cum_obs[j] - cum_obs[i - 1],
                        exp_expr=hl.max(cum_exp[j] - cum_exp[i - 1], 1e-09),
                    )
                    + get_dpois_expr(
                        cond_expr=True,
                        section_oe_expr=total_oe,
                        obs_expr=cum_obs[-1] - cum_obs[j],
                        exp_expr=hl.max(cum_exp[-1] - cum_exp[j], 1e-09),
                    )
                )
            )
        ),
    )


def search_for_two_breaks(
    group_ht: hl.Table,
    count: int,
    chisq_threshold: float = CHISQ_THRESHOLDS["simul"],
    min_num_exp_mis: float = 10,
    min_chisq_threshold: float = MIN_CHISQ_THRESHOLD,
    save_chisq_ht: bool = False,
    freeze: int = CURRENT_FREEZE,
) -> hl.Table:
    """
    Search for transcripts/transcript sections with simultaneous breaks.

    This function searches for breaks for all possible window sizes that exceed a minimum threshold of expected missense variants.

    :param group_ht: Input Table aggregated by transcript/transcript section with lists of cumulative observed
        and expected missense values. HT is filtered to contain only transcript/sections without
        a single significant breakpoint.
    :param count: Which transcript or transcript section group is being run (based on counter generated in `main`).
    :param chisq_threshold: Chi-square significance threshold. Default is
        CHISQ_THRESHOLDS['simul'].
        Default value used in ExAC was 13.8, which corresponds to a p-value of 0.001
        with 2 degrees of freedom.
        (https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm)
    :param min_num_exp_mis: Minimum expected missense value for all three windows defined by two possible
        simultaneous breaks.
    :param min_chisq_threshold: Minimum chi square value to emit from search.
        Default is MIN_CHISQ_THRESHOLD,
        which corresponds to a p-value of 0.025 with 2 degrees of freedom.
    :param save_chisq_ht: Whether to save HT with chi square values annotated for every locus
        (as long as chi square value is >= min_chisq_threshold).
        This saves a lot of extra data and should only occur during the initial search round.
        Default is False.
    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: Table filtered to transcript/sections with significant simultaneous breakpoints
        and annotated with breakpoint information.
    """
    # Create ArrayExpression of StructExpressions that stores each
    # pair of positions (window breakpoints) and their corresponding chi square value
    break_values_expr = hl.range(
        group_ht.start_idx.i_start, group_ht.i_max_idx
    ).flatmap(
        lambda i: hl.range(group_ht.start_idx.j_start, group_ht.j_max_idx).map(
            lambda j: hl.struct(
                i=i,
                j=j,
                chisq=hl.nanmax(
                    calculate_window_chisq(
                        group_ht.max_idx,
                        i,
                        j,
                        group_ht.cum_obs,
                        group_ht.cum_exp,
                        group_ht.total_oe,
                        min_num_exp_mis,
                    ),
                    -1,
                ),
            )
        )
    )
    # Filter ArrayExpression to only keep chi square values that are at least
    # some minimum value (to shorten this array and save on storage)
    break_values_expr = break_values_expr.filter(
        lambda x: x.chisq >= min_chisq_threshold
    )
    group_ht = group_ht.annotate(break_values=break_values_expr)
    # Extract the positions with the maximum chi square value (the "best" break)
    group_ht = group_ht.annotate(
        best_break=group_ht.break_values.fold(
            lambda x, y: hl.if_else(x.chisq >= y.chisq, x, y),
            hl.struct(i=-1, j=-1, chisq=-99),
        )
    )
    group_ht = group_ht.transmute(
        max_chisq=group_ht.best_break.chisq,
        # Adjust breakpoint inclusive/exclusiveness to be consistent with single break breakpoints, i.e.
        # so that the breakpoint site itself is the last site in the left subsection. Thus, the resulting
        # subsections will be divided as follows:
        # [section_start, first_break_pos] (first_break_pos, second_break_pos] (second_break_pos, section_end]
        breakpoints=(
            group_ht.positions[group_ht.best_break.i] - 1,
            group_ht.positions[group_ht.best_break.j],
        ),
    )
    if save_chisq_ht:
        group_ht = group_ht.checkpoint(
            f"{SIMUL_BREAK_TEMP_PATH}/freeze{freeze}_dataproc_temp_chisq_group{count}.ht",
            overwrite=True,
        )
    else:
        group_ht = group_ht.checkpoint(
            f"{TEMP_PATH_WITH_FAST_DEL}/freeze{freeze}_dataproc_temp_chisq_group{count}.ht",
            overwrite=True,
        )
    # Remove rows with maximum chi square values below the threshold
    group_ht = group_ht.filter(group_ht.max_chisq >= chisq_threshold)
    return group_ht


def process_section_group(
    ht_path: str,
    section_group: List[str],
    count: int,
    search_num: int,
    over_threshold: bool,
    output_ht_path: str,
    output_n_partitions: int = 10,
    chisq_threshold: float = CHISQ_THRESHOLDS["simul"],
    min_num_exp_mis: float = MIN_EXP_MIS,
    split_list_len: int = 500,
    read_if_exists: bool = False,
    save_chisq_ht: bool = False,
    freeze: int = CURRENT_FREEZE,
) -> None:
    """
    Run two simultaneous breaks search on a group of transcripts or transcript sections.

    Designed for use with Hail Batch.

    :param ht_path: Path to input Table (Table written using `group_no_single_break_found_ht`).
    :param section_group: List of transcripts or transcript sections to process.
    :param count: Which transcript or transcript section group is being run (based on counter generated in `main`).
    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param over_threshold: Whether input transcript/sections have more
        possible missense sites than threshold specified in `run_simultaneous_breaks`.
    :param output_ht_path: Path to output results Table.
    :param output_n_partitions: Desired number of partitions for output Table.
        Default is 10.
    :param chisq_threshold: Chi-square significance threshold. Default is
        CHISQ_THRESHOLDS['simul'].
        Default value used in ExAC was 13.8, which corresponds to a p-value of 0.001
        with 2 degrees of freedom.
        (https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm)
    :param min_num_exp_mis: Minimum expected missense value for all three windows defined by two possible
        simultaneous breaks.
        Default is MIN_EXP_MIS.
    :param split_list_len: Max length to divide transcript/sections observed or expected missense and position lists into.
        E.g., if split_list_len is 500, and the list lengths are 998, then the transcript/section will be
        split into two rows with lists of length 500 and 498.
        Only used if over_threshold is True. Default is 500.
    :param read_if_exists: Whether to read temporary Table if it already exists rather than overwrite.
        Only applies to Table that is input to `search_for_two_breaks`
        (`f"{temp_ht_path}/{transcript_group[0]}_prep.ht"`).
        Default is False.
    :param save_chisq_ht: Whether to save HT with chi square values annotated for every locus
        (as long as chi square value is >= min_chisq_threshold).
        This saves a lot of extra data and should only occur during the initial search round.
        Default is False.
    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: None; processes Table and writes to path. Also writes success TSV to path.
    """
    ht = hl.read_table(ht_path)
    ht = ht.filter(hl.literal(section_group).contains(ht.section))

    if over_threshold:
        # If transcript/sections longer than threshold, split into multiple rows
        # Each row conducts an independent search within a non-overlapping portion of the transcript/section
        # E.g., if a transcript has 1003 possible missense sites, and the `split_list_len` is 500,
        # then this section will split that transcript into 9 rows, with the following windows:
        # [i_start=0, j_start=0], [i_start=0, j_start=500], [i_start=0, j_start=1000],
        # [i_start=500, j_start=0], [i_start=500, j_start=500], [i_start=500, j_start=1000],
        # [i_start=1000, j_start=0], [i_start=1000, j_start=500], [i_start=1000, j_start=1000]
        # The only windows that should be kept and processed are:
        # [i_start=0, j_start=0], [i_start=0, j_start=500], [i_start=0, j_start=1000],
        # [i_start=500, j_start=500], [i_start=500, j_start=1000], [i_start=1000, j_start=1000]
        ht = ht.annotate(
            start_idx=hl.flatmap(
                lambda i: hl.map(
                    lambda j: hl.struct(i_start=i, j_start=j),
                    hl.range(0, ht.max_idx + 1, split_list_len),
                ),
                hl.range(0, ht.max_idx + 1, split_list_len),
            )
        )
        # Remove entries in `start_idx` struct where j_start is smaller than i_start
        ht = ht.annotate(
            start_idx=hl.filter(lambda x: x.j_start >= x.i_start, ht.start_idx)
        )
        ht = ht.explode("start_idx")
        ht = ht.annotate(i=ht.start_idx.i_start, j=ht.start_idx.j_start)
        ht = ht._key_by_assert_sorted("section", "i", "j")
        ht = ht.annotate(
            # NOTE: i_max_idx needs to be adjusted here to be one smaller than the max
            # This is because we don't need to check the situation where i is the last index in a list
            # For example, if the section has 1003 possible missense sites,
            # (1002 is the largest list index)
            # we don't need to check the scenario where i = 1002
            i_max_idx=hl.min(ht.i + split_list_len, ht.max_idx) - 1,
            j_max_idx=hl.min(ht.j + split_list_len, ht.max_idx),
        )
        # Adjust j_start in rows where j_start is the same as i_start
        ht = ht.annotate(
            start_idx=ht.start_idx.annotate(
                j_start=hl.if_else(
                    ht.start_idx.i_start == ht.start_idx.j_start,
                    ht.start_idx.j_start + 1,
                    ht.start_idx.j_start,
                ),
            ),
        )
        n_rows = ht.count()
        # NOTE: added repartition here because repartioning on read did not
        # keep the desired number of partitions
        # (would sometimes repartition to a lower number of partitions)
        ht = ht.repartition(n_rows)
        prep_path = simul_search_round_bucket_path(
            search_num=search_num,
            bucket_type="prep",
            freeze=freeze,
        )
        ht = ht.checkpoint(
            f"{prep_path}/{section_group[0]}.ht",
            overwrite=not read_if_exists,
            _read_if_exists=read_if_exists,
        )
    else:
        # Add start_idx struct with i_start, j_start, i_max_idx, j_max_idx annotations
        # (these are expected by `search_for_two_breaks`)
        ht = ht.annotate(
            start_idx=hl.struct(i_start=0, j_start=1),
            # Adjusting i_max_idx here to be ht.max_idx - 1 (see note above)
            i_max_idx=ht.max_idx - 1,
            j_max_idx=ht.max_idx,
        )

    # Search for two simultaneous breaks
    ht = search_for_two_breaks(
        group_ht=ht,
        count=count,
        chisq_threshold=chisq_threshold,
        min_num_exp_mis=min_num_exp_mis,
        save_chisq_ht=save_chisq_ht,
    )

    # If over threshold, checkpoint HT and check if there were any breaks
    if over_threshold:
        raw_path = simul_search_round_bucket_path(
            search_num=search_num,
            bucket_type="raw_results",
            freeze=freeze,
        )
        ht = ht.checkpoint(f"{raw_path}/{section_group[0]}.ht", overwrite=True)
        # If any rows had a significant breakpoint,
        # find the one "best" breakpoint (breakpoint with largest chi square value)
        if ht.count() > 0:
            ht = get_max_chisq_per_group(ht, "section", "max_chisq")
            ht = ht.filter(ht.max_chisq == ht.section_max_chisq)

    ht = ht.annotate_globals(chisq_threshold=chisq_threshold)
    ht = ht.naive_coalesce(output_n_partitions)
    # TODO: Restructure ht to match locus-level formats from single breaks
    ht.write(output_ht_path, overwrite=True)
    # TODO: Consider whether we want to write out temp information on chisq values for each potential break combination
    #   like in single breaks

    success_tsvs_path = simul_search_round_bucket_path(
        search_num=search_num,
        bucket_type="success_files",
        freeze=freeze,
    )
    for section in section_group:
        tsv_path = f"{success_tsvs_path}/{section}.tsv"
        with hl.hadoop_open(tsv_path, "w") as o:
            o.write(f"{section}\n")
