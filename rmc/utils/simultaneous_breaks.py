"""This script contains functions used to search for two simultaneous breaks."""
from collections.abc import Callable
import logging
from typing import Optional, Tuple

import hail as hl

from gnomad.utils.reference_genome import get_reference_genome
from gnomad.resources.resource_utils import DataException

from rmc.resources.basics import not_one_break_grouped
from rmc.utils.constraint import get_dpois_expr, get_obs_exp_expr
from rmc.utils.generic import get_avg_bases_between_mis


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("search_for_two_breaks_utils")
logger.setLevel(logging.INFO)


def group_not_one_break_ht(
    ht: hl.Table,
    transcript_ht: hl.Table,
    get_min_window_size: bool = True,
    get_total_exome_bases: bool = False,
    get_total_gnomad_missense: bool = False,
    min_window_size: Optional[int] = None,
    min_num_obs: Optional[int] = None,
) -> None:
    """
    Group HT with transcripts that don't have a single significant break by transcript and collect annotations into lists.

    This creates the input to the two simultaneous breaks search (`search_for_two_breaks`).

    .. note::
        Expects that input Table is keyed by locus and transcript.
        This is *required*, as the function expects that loci are sorted in the input Table.

    :param hl.Table ht: Input Table with transcript that didn't have a single significant break.
    :param hl.Table transcript_ht: Table with start and end positions per transcript.
    :param bool get_min_window_size: Determine minimum window size for two simultaneous breaks search.
        Default is True. Must be True if min_window_size is None.
    :param bool get_total_exome_bases: Get total number of bases in the exome.
        If True, will pull default value from TOTAL_EXOME_BASES (in basics.py).
        Default is False.
    :param bool get_total_gnomad_missense: Get total number of missense variants in gnomAD.
        If True, will pull default value from TOTAL_GNOMAD_MISSENSE (in basics.py).
        Default is False.
    :param Optional[int] min_window_size: Minimum window size for two simultaneous breaks search.
        Must be specified if get_min_window_size is False.
        Default is None.
    :param Optional[int] min_num_obs: Number of observed variants. Used when determining the smallest possible window size.
        Must be specified if get_min_window_size is False.
        Default is None.
    :return: None; writes Table grouped by transcript with cumulative observed, expected missense counts
        and all positions collected into lists to resource path.
    """
    if not get_min_window_size and not min_num_obs:
        raise DataException(
            "min_num_obs must be specified if get_min_window_size is True!"
        )
    if not min_window_size and not get_min_window_size:
        raise DataException(
            "min_window_size must be specified if get_min_window_size is False!"
        )
    if get_min_window_size and min_window_size:
        logger.warning(
            "get_min_window_size is True but min_window_size was also specified. Proceeding with specified min_window_size value..."
        )

    # Get number of base pairs needed to observe `num` number of missense variants (on average)
    # This number is used to determine the min_window_size - which is the smallest allowed distance between simultaneous breaks
    min_window_size = (
        (
            get_avg_bases_between_mis(
                get_reference_genome(ht.locus).name,
                get_total_exome_bases,
                get_total_gnomad_missense,
            )
            * min_num_obs
        )
        if get_min_window_size
        else min_window_size
    )
    logger.info(
        "Minimum window size (window size needed to observe %i missense variants on average): %i",
        min_num_obs,
        min_window_size,
    )

    # Aggregating values into a struct here to force the positions and observed, expected missense values to stay sorted
    # `hl.agg.collect` does not guarantee order: https://hail.is/docs/0.2/aggregators.html#hail.expr.aggregators.collect
    group_ht = ht.group_by("transcript").aggregate(
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
        total_oe=hl.agg.take(ht.overall_oe, 1)[0],
    )
    group_ht = group_ht.annotate_globals(min_window_size=min_window_size)
    group_ht = group_ht.annotate(max_idx=hl.len(group_ht.values.positions) - 1)
    group_ht = group_ht.annotate(
        transcript_start=transcript_ht[group_ht.key].start,
        transcript_end=transcript_ht[group_ht.key].stop,
    )
    group_ht = group_ht.transmute(
        cum_obs=group_ht.values.cum_obs,
        cum_exp=group_ht.values.cum_exp,
        positions=group_ht.values.positions,
    )
    group_ht.write(not_one_break_grouped.path, overwrite=True)


def calculate_window_chisq(
    max_idx: int,
    i: int,
    j: int,
    cum_obs: hl.expr.ArrayExpression,
    cum_exp: hl.expr.ArrayExpression,
    total_oe: hl.expr.Float64Expression,
) -> hl.expr.Float64Expression:
    """
    Calculate chi square significance value for each possible simultaneous breaks window.

    Used only when calculating simultaneous breaks.

    Chi square formula: 2 * (hl.log10(total_alt) - hl.log10(total_null))

    :param int max_idx: Largest list index value.
    :param int i: Smaller list index value corresponding to the smaller position of the two break window.
    :param int j: Larger list index value corresponding to the larger position of the two break window.
    :param hl.expr.ArrayExpression cum_obs: List containing cumulative observed missense values.
    :param hl.expr.ArrayExpression cum_exp: List containing cumulative expected missense values.
    :param expr.Float64Expression total_oe: Transcript overall observed/expected (OE) missense ratio.
    :return: Chi square significance value.
    """
    return (
        hl.case()
        .when(
            # Return -1 when the window spans the entire transcript
            (i == 0) & (j == max_idx),
            -1,
        )
        .when(
            # If i index is the smallest position (anchored at one end of transcript),
            # there are only two transcript subsections: [start_pos, pos[j]], (pos[j], end_pos]
            i == 0,
            (
                2
                * (
                    # Create alt distribution
                    hl.log10(
                        # Create alt distribution for section [start_pos, pos[j]]
                        # The missense values for this section are just the cumulative values at index j
                        get_dpois_expr(
                            cond_expr=True,
                            section_oe_expr=get_obs_exp_expr(
                                True, cum_obs[j], cum_exp[j]
                            ),
                            obs_expr=cum_obs[j],
                            exp_expr=cum_exp[j],
                        )
                        # Create alt distribution for section (pos[j], end_pos]
                        # The missense values for this section are the cumulative values at the last index
                        # minus the values at index j
                        * get_dpois_expr(
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
                    - hl.log10(
                        get_dpois_expr(
                            cond_expr=True,
                            section_oe_expr=total_oe,
                            obs_expr=cum_obs[j],
                            exp_expr=cum_exp[j],
                        )
                        * get_dpois_expr(
                            cond_expr=True,
                            section_oe_expr=total_oe,
                            obs_expr=cum_obs[-1] - cum_obs[j],
                            exp_expr=hl.max(cum_exp[-1] - cum_exp[j], 1e-09),
                        )
                    )
                )
            ),
        )
        .when(
            # If j index is anchored at the largest position, there are two transcript subsections:
            # [start_pos, pos[i]), [pos[i], end_pos]
            j == max_idx,
            (
                2
                * (
                    # Create alt distribution
                    hl.log10(
                        # Create alt distribution for section [start_pos, pos[i])
                        # The missense values for this section are the cumulative values at
                        # one index smaller than index i
                        get_dpois_expr(
                            cond_expr=True,
                            section_oe_expr=get_obs_exp_expr(
                                True, cum_obs[i - 1], cum_exp[i - 1]
                            ),
                            obs_expr=cum_obs[i - 1],
                            exp_expr=cum_exp[i - 1],
                        )
                        # Create alt distribution for section [pos[i], end_pos]
                        # The missense values for this section are the cumulative values at
                        # the last index minus the cumulative values at index i - 1
                        * get_dpois_expr(
                            cond_expr=True,
                            section_oe_expr=get_obs_exp_expr(
                                True,
                                (cum_obs[-1] - cum_obs[i - 1]),
                                hl.max(cum_exp[-1] - cum_exp[i - 1], 1e-09),
                            ),
                            obs_expr=cum_obs[-1] - cum_obs[i - 1],
                            exp_expr=hl.max(cum_exp[-1] - cum_exp[i - 1], 1e-09),
                        )
                    )
                    # Create null distribution
                    - hl.log10(
                        get_dpois_expr(
                            cond_expr=True,
                            section_oe_expr=total_oe,
                            obs_expr=cum_obs[i - 1],
                            exp_expr=cum_exp[i - 1],
                        )
                        * get_dpois_expr(
                            cond_expr=True,
                            section_oe_expr=total_oe,
                            obs_expr=cum_obs[-1] - cum_obs[i - 1],
                            exp_expr=hl.max(cum_exp[-1] - cum_exp[i - 1], 1e-09),
                        )
                    )
                )
            ),
        )
        .default(
            # Neither index is the smallest or largest position,
            # so there are three transcript subsections:
            # [start_pos, pos[i]), [pos[i], pos[j]], (pos[j], end_pos]
            (
                2
                * (
                    # Create alt distribution
                    hl.log10(
                        # Create alt distribution for section [start_pos, pos[i])
                        # The missense values for this section are the cumulative values at
                        # one index smaller than index i
                        get_dpois_expr(
                            cond_expr=True,
                            section_oe_expr=get_obs_exp_expr(
                                True, cum_obs[i - 1], cum_exp[i - 1]
                            ),
                            obs_expr=cum_obs[i - 1],
                            exp_expr=cum_exp[i - 1],
                        )
                        # Create alt distribution for section [pos[i], pos[j]]
                        # The missense values for this section are the cumulative values at index j
                        # minus the cumulative values at index i -1
                        * get_dpois_expr(
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
                        * get_dpois_expr(
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
                    - hl.log10(
                        get_dpois_expr(
                            cond_expr=True,
                            section_oe_expr=total_oe,
                            obs_expr=cum_obs[i - 1],
                            exp_expr=cum_exp[i - 1],
                        )
                        * get_dpois_expr(
                            cond_expr=True,
                            section_oe_expr=total_oe,
                            obs_expr=cum_obs[j] - cum_obs[i - 1],
                            exp_expr=hl.max(cum_exp[j] - cum_exp[i - 1], 1e-09),
                        )
                        * get_dpois_expr(
                            cond_expr=True,
                            section_oe_expr=total_oe,
                            obs_expr=cum_obs[-1] - cum_obs[j],
                            exp_expr=hl.max(cum_exp[-1] - cum_exp[j], 1e-09),
                        )
                    )
                )
            )
        )
    )


def search_for_two_breaks(
    group_ht: hl.Table, chisq_threshold: float = 9.2,
) -> hl.Table:
    """
    Search for windows of constraint in transcripts with simultaneous breaks.

    This function searches for breaks for all possible window sizes but only keeps break sizes >= `min_window_size`.
    `min_window_size` is the number of base pairs needed, on average, to see 10 missense variants (by default).
    For gnomAD v2.1, `min_window_size` is 100bp.

    :param hl.Table ht: Input Table aggregated by transcript with lists of cumulative observed and expected
        missense values. HT is filtered to contain only transcripts with simultaneous breaks.
    :param float chisq_threshold:  Chi-square significance threshold. Default is 9.2.
        This value corresponds to a p-value of 0.99 with 2 degrees of freedom.
        (https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm)
        Default value used in ExAC was 13.8, which corresponds to a p-value of 0.999.
    :return: Table with largest simultaneous break window size annotated per transcript.
    :rtype: hl.Table
    """

    def _simul_break_loop(
        loop_continue: Callable,
        i: int,
        j: int,
        max_idx: hl.expr.Int32Expression,
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
        :param int j: Larger list index value. This index defines the current position of the first break.
            It's the `j` in 3 windows defined by intervals: [start, i), [i, j], (j, end].
        :param hl.expr.Int32Expression: Largest list index for transcript.
        :param float cur_max_chisq: Current maximum chi square value.
        :param int cur_best_i: Current best index i.
        :param int cur_best_j: Current best index j.
        :return: Maximum chi square significance value and optimum index pair i, j.
        """
        # Calculate chi squared value associated with transcript subections created using this index pair i, j
        chisq = calculate_window_chisq(
            max_idx, i, j, group_ht.cum_obs, group_ht.cum_exp, group_ht.total_oe
        )

        # Update current best indices and chi square if new chi square (calculated above)
        # is better than the current stored value (`cur_max_chisq`)
        cur_best_i = hl.if_else(chisq > cur_max_chisq, i, cur_best_i)
        cur_best_j = hl.if_else(chisq > cur_max_chisq, j, cur_best_j)
        cur_max_chisq = hl.max(chisq, cur_max_chisq)

        return hl.if_else(
            # At the end of the iteration through the position list
            # (when index i is at the second to last index of the list),
            # return the best indices.
            # Note that this is the second to last index because all of the windows created where i is the last index
            # were already checked in previous iterations of the loop
            i == (max_idx - 1),
            (cur_max_chisq, cur_best_i, cur_best_j),
            # If we haven't reached the end of the position list with index i,
            # continue with the loop
            hl.if_else(
                j == max_idx,
                # At end of j iteration, continue to next i index
                # Set i to i+1 and j to i+2 (so that the j index is always greater than the i index)
                loop_continue(
                    i + 1, i + 2, max_idx, cur_max_chisq, cur_best_i, cur_best_j
                ),
                # Otherwise, if j hasn't gotten to the maximum index,
                # continue to the next j value for current i
                loop_continue(i, j + 1, max_idx, cur_max_chisq, cur_best_i, cur_best_j),
            ),
        )

    group_ht = group_ht.annotate(
        max_break=hl.experimental.loop(
            _simul_break_loop,
            hl.ttuple(hl.tfloat, hl.tint, hl.tint),
            0,
            1,
            # Have hail log transcript at start of each loop
            hl._console_log(group_ht.transcript, group_ht.max_idx),
            0.0,
            0,
            0,
        )
    )
    group_ht = group_ht.transmute(
        max_chisq=group_ht.max_break[0],
        start_pos=group_ht.positions[group_ht.max_break[1]],
        end_pos=group_ht.positions[group_ht.max_break[2]],
    )
    # Remove rows with maximum chi square values below the threshold
    # or rows where none of the transcript sections is the minimum window size
    group_ht = group_ht.filter(
        # Remove rows with maximum chi square values below the threshold
        (group_ht.max_chisq >= chisq_threshold)
        & (
            (group_ht.end_pos - group_ht.start_pos > group_ht.min_window_size)
            | (group_ht.transcript_end - group_ht.end_pos > group_ht.min_window_size)
            | (
                group_ht.start_pos - group_ht.transcript_start
                > group_ht.min_window_size
            )
        )
    )
    return group_ht
