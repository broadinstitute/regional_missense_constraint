"""
This script calls various functions to run the simultaneous two break search.

This search is done in transcripts that did not have one significant breakpoint.

Note that a couple functions have been copied into this script from `constraint.py`:
- `get_obs_exp_expr`
- `get_dpois_expr`

Note also that a couple functions are contained within this script and not within `simultaneous_breaks.py`:
- `calculate_window_chisq`
- `search_for_two_breaks`
- `process_transcript_group`

This is because python imports do not work in Hail Batch PythonJobs unless
the python scripts are included within the provided Dockerfile, and the scripts within the RMC repo are
too interdependent (would have to copy the entire repo into the Dockerfile).
"""
import argparse
import logging
import subprocess
from tqdm import tqdm

from collections.abc import Callable
from typing import Dict, List, Optional, Tuple, Union

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import (
    LOGGING_PATH,
    no_breaks,
    not_one_break,
    not_one_break_grouped,
    simul_break,
    simul_break_over_threshold,
    simul_break_temp,
    simul_break_under_threshold,
)
from rmc.resources.grch37.reference_data import gene_model
from rmc.slack_creds import slack_token
from rmc.utils.simultaneous_breaks import (
    check_for_successful_transcripts,
    group_not_one_break_ht,
    split_transcripts_by_len,
)


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("search_for_two_breaks")
logger.setLevel(logging.INFO)


def get_obs_exp_expr(
    cond_expr: hl.expr.BooleanExpression,
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
) -> hl.expr.Float64Expression:
    """
    Return observed/expected annotation based on inputs.

    Typically imported from `constraint.py`. See `constraint.py` for full docstring.

    :param hl.expr.BooleanExpression cond_expr: Condition to check prior to adding obs/exp expression.
    :param hl.expr.Int64Expression obs_expr: Expression containing number of observed variants.
    :param hl.expr.Float64Expression exp_expr: Expression containing number of expected variants.
    :return: Observed/expected expression.
    :rtype: hl.expr.Float64Expression
    """
    return hl.or_missing(cond_expr, hl.min(obs_expr / exp_expr, 1))


def get_dpois_expr(
    cond_expr: hl.expr.BooleanExpression,
    section_oe_expr: hl.expr.Float64Expression,
    obs_expr: Union[
        Dict[hl.expr.StringExpression, hl.expr.Int64Expression], hl.expr.Int64Expression
    ],
    exp_expr: Union[
        Dict[hl.expr.StringExpression, hl.expr.Float64Expression],
        hl.expr.Float64Expression,
    ],
) -> hl.expr.StructExpression:
    """
    Calculate null and alt values in preparation for chi-squared test to find significant breaks.

    Typically imported from `constraint.py`. See `constraint.py` for full docstring.

    :param hl.expr.BooleanExpression cond_expr: Conditional expression to check before calculating null and alt values.
    :param hl.expr.Float64Expression section_oe_expr: Expression of section observed/expected value.
    :param Union[Dict[hl.expr.StringExpression, hl.expr.Int64Expression], hl.expr.Int64Expression] obs_expr: Expression containing observed variants count.
    :param Union[Dict[hl.expr.StringExpression, hl.expr.Float64Expression], hl.expr.Float64Expression] exp_expr: Expression containing expected variants count.
    :return: Struct containing forward or reverse null and alt values (either when searching for first or second break).
    :rtype: hl.expr.StructExpression
    """
    return hl.or_missing(cond_expr, hl.dpois(obs_expr, exp_expr * section_oe_expr))


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
        max_idx_i: hl.expr.Int32Expression,
        max_idx_j: hl.expr.Int32Expression,
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
        :param hl.expr.Int32Expression max_idx_i: Largest list index for smaller list index value.
        :param hl.expr.Int32Expression max_idx_j: Largest list index for larger list index value.
        :param float cur_max_chisq: Current maximum chi square value.
        :param int cur_best_i: Current best index i.
        :param int cur_best_j: Current best index j.
        :return: Maximum chi square significance value and optimum index pair i, j.
        """
        # Calculate chi squared value associated with transcript subections created using this index pair i, j
        chisq = calculate_window_chisq(
            group_ht.max_idx,
            i,
            j,
            group_ht.cum_obs,
            group_ht.cum_exp,
            group_ht.total_oe,
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
            # Note that this is the second to last index (max_idx_i - 1) because all of the windows created where i is the last index
            # were already checked in previous iterations of the loop
            i == (max_idx_i - 1),
            (cur_max_chisq, cur_best_i, cur_best_j),
            # If we haven't reached the end of the position list with index i,
            # continue with the loop
            hl.if_else(
                j == max_idx_j,
                # At end of j iteration, continue to next i index
                # Set i to i+1 and j to i+2 (so that the j index is always greater than the i index)
                loop_continue(
                    i + 1,
                    i + 2,
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
            group_ht.i_max_idx,
            group_ht.j_max_idx,
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


def process_transcript_group(
    ht_path: str,
    transcript_group: List[str],
    over_threshold: bool,
    output_ht_path: str,
    output_tsv_path: str,
    temp_ht_path: Optional[str] = None,
    chisq_threshold: float = 9.2,
    split_window_size: int = 500,
) -> None:
    """
    Run two simultaneous breaks search on a group of transcripts.

    Designed for use with Hail Batch.

    :param str ht_path: Path to input Table (Table written using `group_not_one_break_ht`).
    :param List[str] transcript_group: List of transcripts to process.
    :param bool over_threshold: Whether input transcripts have more
        possible missense variants than threshold specified in `run_simultaneous_breaks`.
    :param str output_ht_path: Path to output results Table.
    :param str output_tsv_path: Path to success TSV bucket.
    :param Optional[str] temp_ht_path: Path to temporary Table. Required only if over_thresold is True.
    :param float chisq_threshold: Chi-square significance threshold. Default is 9.2.
        This value corresponds to a p-value of 0.99 with 2 degrees of freedom.
        (https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm)
        Default value used in ExAC was 13.8, which corresponds to a p-value of 0.999.
    :param int split_window_size: Window size to search for transcripts that have more
        possible missense variants than threshold. Only used if over_threshold is True.
    :return: None; processes Table and writes to path. Also writes success TSV to path.
    """
    ht = hl.read_table(ht_path)
    ht = ht.filter(hl.literal(transcript_group).contains(ht.transcript))

    if over_threshold:
        # If transcripts longer than threshold, split transcripts into multiple rows
        # Each row has a window to search
        # E.g., if a transcript has 1000 possible missense variants, and the `split_window_size` is 500,
        # then this section will split that transcript into 9 rows, with the following windows:
        # [i_start=0, j_start=0], [i_start=0, j_start=500], [i_start=0, j_start=1000],
        # [i_start=500, j_start=0], [i_start=500, j_start=500], [i_start=500, j_start=1000],
        # [i_start=1000, j_start=0], [i_start=1000, j_start=500], [i_start=1000, j_start=1000]
        ht = ht.annotate(
            start_idx=hl.flatmap(
                lambda i: hl.map(
                    lambda j: hl.struct(i_start=i, j_start=j),
                    hl.range(0, ht.list_len, split_window_size),
                ),
                hl.range(0, ht.list_len, split_window_size),
            )
        )
        ht = ht.explode("start_idx")
        ht = ht.annotate(i=ht.start_idx.i_start, j=ht.start_idx.j_start)
        ht = ht._key_by_assert_sorted("transcript", "i", "j")
        ht = ht.filter(ht.j >= ht.i)
        ht = ht.annotate(
            i_max_idx=hl.min(ht.i + split_window_size, ht.list_len - 1),
            j_max_idx=hl.min(ht.j + split_window_size, ht.list_len - 1),
        )
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
        ht = ht.repartition(n_rows)
        ht = ht.checkpoint(
            f"{temp_ht_path}/{transcript_group[0]}_prep.ht", overwrite=True
        )
    else:
        # Add start_idx struct with i_start, j_start, i_max_idx, j_max_idx annotations
        # (these are expected by `search_for_two_breaks`)
        ht = ht.annotate(
            start_idx=hl.struct(i_start=0, j_start=0),
            i_max_idx=ht.max_idx,
            j_max_idx=ht.max_idx,
        )

    # Search for two simultaneous breaks
    ht = search_for_two_breaks(ht, chisq_threshold)

    # If over threshold, checkpoint HT and check if there were any breaks
    if over_threshold:
        ht = ht.checkpoint(
            f"{temp_ht_path}/{transcript_group[0]}_res.ht", overwrite=True
        )
        # If any rows had a significant breakpoint,
        # find the one "best" breakpoint (breakpoint with largest chi square value)
        if ht.count() > 0:
            max_chisq = ht.aggregate(hl.agg.max(ht.max_chisq))
            ht = ht.filter(ht.max_chisq == max_chisq)

    ht.write(output_ht_path, overwrite=True)

    for transcript in transcript_group:
        success_tsv_path = f"{output_tsv_path}/{transcript}.tsv"
        with hl.hadoop_open(success_tsv_path, "w") as o:
            o.write(f"{transcript}\n")


def main(args):
    """Search for two simultaneous breaks in transcripts without evidence of a single significant break."""
    if not args.command:
        raise DataException("Please specify command for this script!")

    try:

        if args.command == "create-grouped-ht":
            logger.warning("This step should be run in Dataproc!")
            hl.init(log="/search_for_two_breaks_create_grouped_ht.log")

            if not file_exists(not_one_break_grouped.path) or args.create_grouped_ht:
                if args.min_num_obs == 0:
                    # Make sure user didn't specify a min obs of 0
                    raise DataException(
                        "Minimum number of observed variants must be greater than zero!"
                    )

                logger.info(
                    "Creating grouped HT with lists of cumulative observed and expected missense values..."
                )
                group_not_one_break_ht(
                    ht=not_one_break.ht(),
                    transcript_ht=gene_model.ht(),
                    get_min_window_size=args.get_min_window_size,
                    get_total_exome_bases=args.get_total_exome_bases,
                    get_total_gnomad_missense=args.get_total_gnomad_missense,
                    min_window_size=args.min_window_size,
                    min_num_obs=args.min_num_obs,
                )

        if args.command == "split-transcripts":
            logger.warning("This step should be run in Dataproc!")
            hl.init(log="/search_for_two_breaks_split_Transcripts.log")
            split_transcripts_by_len(
                ht=not_one_break_grouped.ht(),
                transcript_len_threshold=args.transcript_len_threshold,
                ttn_id=args.ttn,
            )

        if args.command == "run-batches":
            import hailtop.batch as hb

            logger.warning("This step should be run locally!")
            hl.init(log="search_for_two_breaks_run_batches.log")

            # Make sure custom machine wasn't specified with under threshold
            if args.under_threshold and args.use_custom_machine:
                raise DataException(
                    "Do not specify --use-custom-machine when transcripts are --under-threshold size!"
                )

            logger.info("Importing SetExpression with transcripts...")
            if not args.under_threshold and not args.over_threshold:
                raise DataException(
                    "Must specify if transcript sizes are --under-threshold or --over-threshold!"
                )

            transcripts_to_run = check_for_successful_transcripts(
                transcripts=(
                    list(
                        hl.eval(
                            hl.experimental.read_expression(simul_break_under_threshold)
                        )
                    )
                    if args.under_threshold
                    else list(
                        hl.eval(
                            hl.experimental.read_expression(simul_break_over_threshold)
                        )
                    )
                ),
            )
            logger.info("Found %i transcripts to search...", len(transcripts_to_run))

            logger.info("Setting up Batch parameters...")
            backend = hb.ServiceBackend(
                billing_project=args.billing_project,
                remote_tmpdir=args.batch_bucket,
                google_project=args.google_project,
            )

            if args.use_custom_machine:
                b = hb.Batch(
                    name="simul_breaks",
                    backend=backend,
                    default_python_image=args.docker_image,
                )
            else:
                b = hb.Batch(
                    name="simul_breaks",
                    backend=backend,
                    default_memory=args.batch_memory,
                    default_cpu=args.batch_cpu,
                    default_storage=args.batch_storage,
                    default_python_image=args.docker_image,
                )

            if args.under_threshold:
                transcript_groups = [
                    transcripts_to_run[x : x + args.group_size]
                    for x in range(0, len(transcripts_to_run), args.group_size)
                ]
                split_window_size = None
                count = 1
                for group in tqdm(transcript_groups, unit="transcript group"):
                    logger.info("Working on group number %s...", count)
                    logger.info(group)
                    job_name = (
                        f'group{count}{"over" if args.over_threshold else "under"}'
                    )
                    j = b.new_python_job(name=job_name)
                    j.call(
                        process_transcript_group,
                        not_one_break_grouped.path,
                        group,
                        args.over_threshold,
                        f"{simul_break_temp}/hts/simul_break_{job_name}.ht",
                        f"{simul_break_temp}/success_files",
                        None,
                        args.chisq_threshold,
                        split_window_size,
                    )
                    count += 1

            else:
                # transcript_groups = [[transcript] for transcript in transcripts_to_run]
                transcript_groups = [["ENST00000301030"]]
                for group in transcript_groups:
                    j = b.new_python_job(name=group[0])
                    if args.use_custom_machine:
                        j = b.new_python_job(name=job_name)
                        j._machine_type = "n1-highmem-32"
                        j._preemptible = True
                        j.storage("100Gi")
                    j.call(
                        process_transcript_group,
                        "gs://gnomad-tmp/kc/test_100_over5k.ht",
                        # not_one_break_grouped.path,
                        group,
                        args.over_threshold,
                        f"{simul_break_temp}/hts/simul_break_{group[0]}.ht",
                        f"{simul_break_temp}/success_files",
                        f"{simul_break_temp}",
                        args.chisq_threshold,
                        args.group_size,
                    )

            b.run(wait=False)

        if args.command == "run-batches-dataproc":
            logger.warning(
                "This step should be run in Dataproc on an autoscaling cluster!"
            )
            hl.init(log="/search_for_two_breaks_run_batches_dataproc.log")
            transcript_groups = [
                [transcript] for transcript in args.dataproc_transcripts.split(",")
            ]

            for group in transcript_groups:
                process_transcript_group(
                    ht_path=not_one_break_grouped.path,
                    transcript_group=group,
                    over_threshold=True,
                    output_ht_path=f"{simul_break_temp}/hts/simul_break_{group[0]}.ht",
                    output_tsv_path=f"{simul_break_temp}/success_files",
                    temp_ht_path=f"{simul_break_temp}",
                    chisq_threshold=args.chisq_threshold,
                    split_window_size=args.group_size,
                )

        if args.command == "verify-transcripts":
            logger.warning("This step should be run locally!")
            hl.init(log="search_for_two_breaks_verify_transcripts.log")

            logger.info("Verifying that all transcripts were processed...")
            transcripts = list(
                hl.eval(
                    hl.experimental.read_expression(simul_break_under_threshold).union(
                        hl.experimental.read_expression(simul_break_over_threshold)
                    )
                )
            )
            missing_transcripts = check_for_successful_transcripts(transcripts)
            if len(missing_transcripts) > 0:
                logger.error(missing_transcripts)
                raise DataException(
                    f"{len(missing_transcripts)} are missing! Please rerun."
                )

            # Check if TTN was run and print a warning if it wasn't
            # TTN ID isn't included in `simul_break_under_threshold` or `simul_break_over_threshold`
            # It needs to be run separately due to its size
            logger.info("Checking if TTN was processed...")
            ttn_missing = check_for_successful_transcripts(transcripts=[args.ttn])
            if len(ttn_missing) > 0:
                logger.warning(
                    "TTN wasn't processed successfully. Double check whether this is expected!"
                )
            logger.info("Done searching for transcript success TSVS!")

        if args.command == "merge-hts":
            logger.warning("This step should be run in Dataproc!")
            hl.init(log="/search_for_two_breaks_merge_hts.log")

            logger.info("Collecting all HT paths...")
            intermediate_hts = []
            ht_bucket = f"{simul_break_temp}/hts/"
            temp_ht_paths = (
                subprocess.check_output(["gsutil", "ls", ht_bucket])
                .decode("utf8")
                .strip()
                .split("\n")
            )
            for ht_path in temp_ht_paths:
                ht_path = ht_path.strip("/")
                temp = hl.read_table(ht_path)
                intermediate_hts.append(temp)

            ht = intermediate_hts[0].union(*intermediate_hts[1:])
            ht = ht.checkpoint(simul_break.path, overwrite=args.overwrite)
            logger.info("Wrote simultaneous breaks HT with %i lines", ht.count())

            # Collect all transcripts with two simultaneous breaks
            simul_break_transcripts = ht.aggregate(
                hl.agg.collect_as_set(ht.transcript),
            )
            logger.info(
                "%i transcripts had two simultaneous breaks",
                len(simul_break_transcripts),
            )
            simul_break_transcripts = hl.literal(simul_break_transcripts)

            logger.info(
                "Getting transcripts with no evidence of regional missense constraint..."
            )
            context_ht = not_one_break.ht()
            context_ht = context_ht.filter(
                ~simul_break_transcripts.contains(context_ht.transcript)
            )
            context_ht.write(no_breaks.path, overwrite=args.overwrite)

    finally:
        # Don't copy log if running locally because copy_log only operates on files in gcloud
        if args.command != "run-batches" and args.command != "verify-transcripts":
            logger.info("Copying hail log to logging bucket...")
            hl.copy_log(LOGGING_PATH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This regional missense constraint script two simultaneous breaks in transcripts without evidence of a single significant break."
    )
    parser.add_argument(
        "--chisq-threshold",
        help="Chi-square significance threshold. Value should be 9.2 (value adjusted from ExAC code due to discussion with Mark).",
        type=float,
        default=9.2,
    )
    parser.add_argument(
        "--overwrite", help="Overwrite existing data.", action="store_true"
    )
    parser.add_argument(
        "--slack-channel",
        help="Send message to Slack channel/user.",
        default="@kc (she/her)",
    )
    parser.add_argument(
        "--TTN",
        help="TTN transcript ID. TTN is so large that it needs to be treated separately.",
        default="ENST00000589042",
    )

    # Create subparsers for each step
    # Need to specify `dest` to be able to check which subparser is being invoked
    # `dest`: https://docs.python.org/3/library/argparse.html#dest
    subparsers = parser.add_subparsers(title="command", dest="command")

    create_grouped_ht = subparsers.add_parser(
        "create-grouped-ht",
        help="""
        Create hail Table grouped by transcript with cumulative observed and expected missense values collected into lists.
        This step should be run in Dataproc.
        """,
    )
    min_window_group = create_grouped_ht.add_mutually_exclusive_group()
    min_window_group.add_argument(
        "--min-window-size",
        help="Smallest possible window size for simultaneous breaks. Determined by running --get-min-window-size.",
        type=int,
    )
    min_window_group.add_argument(
        "--get-min-window-size",
        help="Determine smallest possible window size for simultaneous breaks.",
        action="store_true",
    )
    create_grouped_ht.add_argument(
        "--get-total-exome-bases",
        help="Get total number of bases in the exome. If not set, will pull default value from TOTAL_EXOME_BASES.",
        action="store_true",
    )
    create_grouped_ht.add_argument(
        "--get-total-gnomad-missense",
        help="Get total number of missense variants in gnomAD. If not set, will pull default value from TOTAL_GNOMAD_MISSENSE.",
        action="store_true",
    )
    create_grouped_ht.add_argument(
        "--min-num-obs",
        help="Number of observed variants. Used when determining the smallest possible window size for simultaneous breaks.",
        type=int,
        default=10,
    )

    split_transcripts = subparsers.add_parser(
        "split-transcripts",
        help="""
        Split transcripts based on number of possible missense positions.
        This is used to create batches of transcripts to run through search for two breaks code.
        This step should be run in Dataproc.
        """,
    )
    split_transcripts.add_argument(
        "--transcript-len-threshold",
        help="Cutoff for number of possible missense positions in transcript. Used to create batches of transcripts.",
        type=int,
        default=5000,
    )

    run_batches = subparsers.add_parser(
        "run-batches",
        help="Run batches of transcripts using Hail Batch. This step should be run locally.",
    )
    transcript_size = run_batches.add_mutually_exclusive_group()
    transcript_size.add_argument(
        "--under-threshold",
        help="Transcripts in batch should have less than --transcript-len-threshold possible missense positions.",
        action="store_true",
    )
    transcript_size.add_argument(
        "--over-threshold",
        help="Transcripts in batch should greater than or equal to --transcript-len-threshold possible missense positions.",
        action="store_true",
    )
    run_batches.add_argument(
        "--group-size",
        help="""
        Number of transcripts to include in each group of transcripts to be submitted to Hail Batch if --under-threshold.
        Size of windows to split transcripts if --over-threshold.
        Default is 100.
        """,
        type=int,
        default=100,
    )
    run_batches.add_argument(
        "--billing-project",
        help="Billing project to use with hail batch.",
        default="gnomad-production",
    )
    run_batches.add_argument(
        "--batch-bucket",
        help="Bucket provided to hail batch for temporary storage.",
        default="gs://gnomad-tmp/kc/",
    )
    run_batches.add_argument(
        "--google-project",
        help="Google cloud project provided to hail batch for storage objects access.",
        default="broad-mpg-gnomad",
    )
    run_batches.add_argument(
        "--use-custom-machine",
        help="Use custom machine for hail batch rather than setting batch memory, cpu, and storage. Only necessary if --over-threshold.",
        action="store_true",
    )
    run_batches.add_argument(
        "--batch-memory",
        help="Amount of memory to request for hail batch jobs.",
        default="standard",
    )
    run_batches.add_argument(
        "--batch-cpu",
        help="Number of CPUs to request for hail batch jobs.",
        default=8,
        type=int,
    )
    run_batches.add_argument(
        "--batch-storage",
        help="Amount of disk storage to request for hail batch jobs.",
        default="10Gi",
    )
    run_batches.add_argument(
        "--docker-image",
        help="""
        Docker image to provide to hail Batch. Must have dill, hail, and python installed.
        Suggested image with --under-threshold: gcr.io/broad-mpg-gnomad/tgg-methods-vm:20220302.

        If running with --over-threshold, Docker image must also contain this line:
        `ENV PYSPARK_SUBMIT_ARGS="--driver-memory 8g --executor-memory 8g pyspark-shell"`
        to make sure the job allocates memory correctly.
        Default image (gcr.io/broad-mpg-gnomad/rmc:20220304) has this line.
        """,
        default="gcr.io/broad-mpg-gnomad/rmc:20220304",
    )

    run_batches_dataproc = subparsers.add_parser(
        "run-batches-dataproc",
        help="""
        Run batches of transcripts using Google Cloud Dataproc.
        This should only be used for transcripts that are over --transcript-len-threshold if they are too slow or getting preempted in
        hail Batch. All transcripts under --transcript-len-threshold should be run in hail Batch.

        If using this step to run TTN, use a large autoscaling cluster (highmem-8, scales to 100 preemptibles).
        (This step should be run in Dataproc.)
        """,
    )
    run_batches_dataproc.add_argument(
        "--dataproc-transcripts", help="Comma separated list of transcript IDs to run."
    )

    verify_transcripts = subparsers.add_parser(
        "verify-transcripts",
        help="Verify that all transcripts were processed using hail Batch. This step should be run locally.",
    )

    merge_hts = subparsers.add_parser(
        "merge-hts",
        help="Merge all intermediate result Tables into a single Table. This step should be run in Dataproc.",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
