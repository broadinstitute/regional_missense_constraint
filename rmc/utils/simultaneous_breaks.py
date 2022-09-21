"""This script contains functions used to search for two simultaneous breaks."""
from collections.abc import Callable
import logging
from typing import List, Optional, Tuple

import hail as hl

from gnomad.utils.file_utils import file_exists, parallel_file_exists

from rmc.resources.basics import SIMUL_BREAK_TEMP_PATH
from rmc.resources.rmc import sections_to_simul_by_threshold_path, simul_search_round_bucket_path
from rmc.utils.constraint import get_dpois_expr, get_obs_exp_expr


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
    group_ht.write(out_ht_path, overwrite=True)


def split_sections_by_len(
    ht_path: str,
    group_str: str,
    search_num: int,
    is_rescue: bool,
    missense_len_threshold: int,
    ttn_id: str,
    overwrite: bool,
) -> None:
    """
    Split transcripts/transcript sections based on the specified number of possible missense sites.

    This is necessary because transcripts/transcript sections with more possible missense sites
    take longer to run through `hl.experimental.loop`.

    :param ht: Path to input Table (Table written using `group_no_single_break_found_ht`).
    :param group_str: Field used to group observed and expected values.
    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param is_rescue: Whether current search is in rescue pathway.
    :param missense_len_threshold: Cutoff based on possible number of missense sites in section.
    :param ttn_id: TTN transcript ID. TTN is large and needs to be processed separately.
    :param overwrite: Whether to overwrite existing SetExpressions.
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
    # TODO: Re-evaluate if this check is necessary for v4
    if ttn_id in list(over_threshold):
        logger.warning(
            "TTN is present in input transcripts! It will need to be run separately."
        )
        over_threshold = list(over_threshold)
        over_threshold.remove(ttn_id)
        over_threshold = set(over_threshold)
    hl.experimental.write_expression(
        under_threshold,
        sections_to_simul_by_threshold_path(
            search_num=search_num,
            is_rescue=is_rescue,
            is_over_threshold=False,
        ),
        overwrite
    )
    hl.experimental.write_expression(
        over_threshold,
        sections_to_simul_by_threshold_path(
            search_num=search_num,
            is_rescue=is_rescue,
            is_over_threshold=True,
        ),
        overwrite
    )


def check_for_successful_transcripts(
    transcripts: List[str],
    search_num: int,
    is_rescue: bool,
    in_parallel: bool = True,
) -> List[str]:
    """
    Check if any transcripts have been previously searched by searching for success TSV existence.

    .. note::
        This step needs to be run locally due to permissions involved with `parallel_file_exists`.

    :param List[str] transcripts: List of transcripts to check.
    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param is_rescue: Whether search is in rescue pathway.
    :param bool in_parallel: Whether to check if successful file exist in parallel.
        If True, must be run locally and not in Dataproc. Default is True.
    :return: List of transcripts didn't have success TSVs and therefore still need to be processed.
    """
    logger.info("Checking if any transcripts have already been searched...")
    success_file_path = simul_search_round_bucket_path(
        is_rescue=is_rescue,
        search_num=search_num,
        bucket_type="success_files",
    )
    transcript_success_map = {}
    transcripts_to_run = []
    for transcript in transcripts:
        transcript_success_map[transcript] = f"{success_file_path}/{transcript}.tsv"

    if in_parallel:
        # Use parallel_file_exists if in_parallel is set to True
        success_tsvs_exist = parallel_file_exists(list(transcript_success_map.values()))
        for transcript in transcripts:
            if not success_tsvs_exist[transcript_success_map[transcript]]:
                transcripts_to_run.append(transcript)
    else:
        # Otherwise, use file_exists
        for transcript in transcripts:
            if not file_exists(transcript_success_map[transcript]):
                transcripts_to_run.append(transcript)

    return transcripts_to_run


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

    Chi square formula: 2 * (hl.log10(total_alt) - hl.log10(total_null))

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
    chisq_threshold: float = 9.2,
    min_num_exp_mis: float = 10,
) -> hl.Table:
    """
    Search for transcripts/transcript sections with simultaneous breaks.

    This function searches for breaks for all possible window sizes that exceed a minimum threshold of expected missense variants.

    :param hl.Table group_ht: Input Table aggregated by transcript/transcript section with lists of cumulative observed
        and expected missense values. HT is filtered to contain only transcript/sections without
        a single significant breakpoint.
    :param float chisq_threshold:  Chi-square significance threshold. Default is 9.2.
        This value corresponds to a p-value of 0.01 with 2 degrees of freedom.
        (https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm)
        Default value used in ExAC was 13.8, which corresponds to a p-value of 0.001.
    :param float min_num_exp_mis: Minimum expected missense value for all three windows defined by two possible
        simultaneous breaks.
    :return: Table filtered to transcript/sections with significant simultaneous breakpoints
        and annotated with breakpoint information.
    :rtype: hl.Table
    """

    def _simul_break_loop(
        loop_continue: Callable,
        i: int,
        j: int,
        start_idx_j: int,
        max_idx_i: int,
        max_idx_j: int,
        cur_max_chisq: float,
        cur_best_i: int,
        cur_best_j: int,
    ) -> Tuple[float, int, int]:
        """
        Iterate over each possible pair of indices in a transcript's cumulative value lists to find the optimum two break window.

        :param Callable[float, int, int] loop_continue: Function to restart hail loop.
            First argument to `hl.experimental.loop` must be a function (`_simul_break_loop` in this case),
            and the first argument to that function must be another function.
            Calling `loop_continue` tells hail to go back to the top of the loop with loop variables updated.
        :param int i: Smaller list index value. This index defines the current position of the first break.
            It's the `i` in 3 windows defined by intervals: [start, i), [i, j], (j, end].
        :param int j: Larger list index value. This index defines the current position of the second break.
            It's the `j` in 3 windows defined by intervals: [start, i), [i, j], (j, end].
        :param int start_idx_j: Smallest list index for larger list index value `j`.
        :param int max_idx_i: Largest list index for smaller list index value.
        :param int max_idx_j: Largest list index for larger list index value.
        :param float cur_max_chisq: Current maximum chi square value.
        :param int cur_best_i: Current best index i.
        :param int cur_best_j: Current best index j.
        :return: Maximum chi square significance value and optimum index pair i, j.
        """
        # Calculate chi squared value associated with transcript subsections created using this index pair i, j
        chisq = calculate_window_chisq(
            group_ht.max_idx,
            i,
            j,
            group_ht.cum_obs,
            group_ht.cum_exp,
            group_ht.total_oe,
            min_num_exp_mis,
        )

        # Make sure chi square isn't NaN
        chisq = hl.nanmax(chisq, -1)

        # Update current best indices and chi square if new chi square (calculated above)
        # is better than the current stored value (`cur_max_chisq`)
        cur_best_i = hl.if_else(chisq > cur_max_chisq, i, cur_best_i)
        cur_best_j = hl.if_else(chisq > cur_max_chisq, j, cur_best_j)
        cur_max_chisq = hl.max(chisq, cur_max_chisq)

        return hl.if_else(
            # Return the best indices at the end of the iteration through the position list
            # Note that max_idx_i has been adjusted to be ht.max_idx - 1 (or i + window_size - 1):
            # see note in `process_section_group`
            # Also note that j needs to be checked here to ensure that j is also at the end of its loop
            # (This check is necessary when transcripts have been split into multiple i, j windows
            # across multiple rows)
            (i == max_idx_i) & (j == max_idx_j),
            (cur_max_chisq, cur_best_i, cur_best_j),
            # If we haven't reached the end of the position list with index i,
            # continue with the loop
            hl.if_else(
                j == max_idx_j,
                # At end of j iteration, continue to next i index
                # Set i to i + 1
                # and set j to the larger value between i + 2 and start index value for j
                # This is to avoid redundant work in larger transcripts; e.g.:
                # start_idx_i = 0, start_idx_j = 50 ->
                # using `hl.max()` here means that j will be reset to 50 rather than 2 on the second
                # iteration of the loop will restart at 50
                # Note that the j index should always be larger than the i index
                loop_continue(
                    i + 1,
                    hl.min(hl.max(i + 2, start_idx_j), max_idx_j),
                    start_idx_j,
                    max_idx_i,
                    max_idx_j,
                    cur_max_chisq,
                    cur_best_i,
                    cur_best_j,
                ),
                # Otherwise, if j hasn't gotten to the maximum index,
                # continue to the next j value for current i
                loop_continue(
                    i,
                    j + 1,
                    start_idx_j,
                    max_idx_i,
                    max_idx_j,
                    cur_max_chisq,
                    cur_best_i,
                    cur_best_j,
                ),
            ),
        )

    group_ht = group_ht.annotate(
        max_break=hl.experimental.loop(
            _simul_break_loop,
            hl.ttuple(hl.tfloat, hl.tint, hl.tint),
            group_ht.start_idx.i_start,
            group_ht.start_idx.j_start,
            group_ht.start_idx.j_start,
            group_ht.i_max_idx,
            group_ht.j_max_idx,
            0.0,
            0,
            0,
        )
    )
    group_ht = group_ht.transmute(
        max_chisq=group_ht.max_break[0],
        # Adjust breakpoint inclusive/exclusiveness to be consistent with single break breakpoints, i.e.
        # so that the breakpoint site itself is the last site in the left subsection. Thus, the resulting
        # subsections will be divided as follows:
        # [section_start, first_break_pos] (first_break_pos, second_break_pos] (second_break_pos, section_end]
        breakpoints=(
            group_ht.positions[group_ht.max_break[1]] - 1,
            group_ht.positions[group_ht.max_break[2]],
        ),
    )
    # Remove rows with maximum chi square values below the threshold
    group_ht = group_ht.filter(group_ht.max_chisq >= chisq_threshold)
    return group_ht


def process_section_group(
    ht_path: str,
    section_group: List[str],
    over_threshold: bool,
    output_ht_path: str,
    success_tsvs_path: str,
    temp_ht_path: Optional[str] = None,
    chisq_threshold: float = 9.2,
    min_num_exp_mis: float = 10,
    split_list_len: int = 500,
    read_if_exists: bool = False,
) -> None:
    """
    Run two simultaneous breaks search on a group of transcripts or transcript sections.

    Designed for use with Hail Batch.

    :param str ht_path: Path to input Table (Table written using `group_no_single_break_found_ht`).
    :param List[str] section_group: List of transcripts or transcript sections to process.
    :param bool over_threshold: Whether input transcript/sections have more
        possible missense sites than threshold specified in `run_simultaneous_breaks`.
    :param str output_ht_path: Path to output results Table.
    :param str success_tsvs_path: Path to success TSV bucket.
    :param Optional[str] temp_ht_path: Path to temporary Table. Required only if over_threshold is True.
    :param float chisq_threshold: Chi-square significance threshold. Default is 9.2.
        This value corresponds to a p-value of 0.01 with 2 degrees of freedom.
        (https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm)
        Default value used in ExAC was 13.8, which corresponds to a p-value of 0.001.
    :param float min_num_exp_mis: Minimum expected missense value for all three windows defined by two possible
        simultaneous breaks.
    :param int split_list_len: Max length to divide transcript/sections observed or expected missense and position lists into.
        E.g., if split_list_len is 500, and the list lengths are 998, then the transcript/section will be
        split into two rows with lists of length 500 and 498.
        Only used if over_threshold is True. Default is 500.
    :param bool read_if_exists: Whether to read temporary Table if it already exists rather than overwrite.
        Only applies to Table that is input to `search_for_two_breaks`
        (`f"{temp_ht_path}/{transcript_group[0]}_prep.ht"`).
        Default is False.
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
        ht = ht.checkpoint(
            f"{temp_ht_path}/{section_group[0]}_prep.ht",
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
    ht = search_for_two_breaks(ht, chisq_threshold, min_num_exp_mis)

    # If over threshold, checkpoint HT and check if there were any breaks
    if over_threshold:
        ht = ht.checkpoint(f"{temp_ht_path}/{section_group[0]}_res.ht", overwrite=True)
        # If any rows had a significant breakpoint,
        # find the one "best" breakpoint (breakpoint with largest chi square value)
        if ht.count() > 0:
            group_ht = ht.group_by("section").aggregate(
                section_max_chisq=hl.agg.max(ht.max_chisq)
            )
            ht = ht.annotate(section_max_chisq=group_ht[ht.section].section_max_chisq)
            ht = ht.filter(ht.max_chisq == ht.section_max_chisq)

    ht = ht.annotate_globals(chisq_threshold=chisq_threshold)
    # TODO: Restructure ht to match locus-level formats from single breaks
    ht.write(output_ht_path, overwrite=True)
    # TODO: Consider whether we want to write out temp information on chisq values for each potential break combination
    #   like in single breaks

    for section in section_group:
        success_tsvs_path = f"{success_tsvs_path}/{section}.tsv"
        with hl.hadoop_open(success_tsvs_path, "w") as o:
            o.write(f"{section}\n")
