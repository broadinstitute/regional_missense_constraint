# TODO: sync code with changes in rmc.utils.simultaneous_breaks (once those changes are finalized)
"""
This script searches for two simultaneous breaks in groups of transcripts using Hail Batch.

Note that a couple functions have been copied into this script from `constraint.py`:
- `get_obs_exp_expr`
- `get_dpois_expr`

Note also that a few functions have been copied into this script from `simultaneous_breaks.py`:
- `calculate_window_chisq`
- `search_for_two_breaks`
- `process_transcript_group`

This is because python imports do not work in Hail Batch PythonJobs unless
the python scripts are included within the provided Dockerfile, and the scripts within the RMC repo are
too interdependent (would have to copy the entire repo into the Dockerfile).

This script should be run locally because it submits jobs to the Hail Batch service.
"""
import argparse
import logging
from tqdm import tqdm

from collections.abc import Callable
from typing import Dict, List, Optional, Tuple, Union

import hail as hl
import hailtop.batch as hb

from gnomad.resources.resource_utils import DataException
from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import SIMUL_BREAK_TEMP_PATH, TEMP_PATH_WITH_DEL
from rmc.resources.rmc import (
    not_one_break_grouped,
    sections_to_simul_by_threshold_path,
)
from rmc.slack_creds import slack_token
from rmc.utils.simultaneous_breaks import check_for_successful_transcripts


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("run_batches")
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
    # Cap the o/e ratio at 1 to avoid pulling out regions that are enriched for missense variation
    # Code is looking for missense constraint, so regions with a ratio of >= 1.0 can be grouped together
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
    Calculate probabilities (in log10 space) of the observed values under a Poisson model.

    Typically imported from `constraint.py`. See `constraint.py` for full docstring.

    :param cond_expr: Conditional expression to check before calculating probability.
    :param section_oe_expr: Expression of section observed/expected value.
    :param obs_expr: Expression containing observed variants count.
    :param exp_expr: Expression containing expected variants count.
    :return: log10 of the probability under Poisson model.
    """
    # log_p = True returns the natural logarithm of the probability density
    # Divide this value by hl.log(10) to convert to log base 10
    return hl.or_missing(
        cond_expr,
        hl.dpois(obs_expr, exp_expr * section_oe_expr, log_p=True) / hl.log(10),
    )


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
    output_tsv_path: str,
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
    :param str output_tsv_path: Path to success TSV bucket.
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
    ht.write(output_ht_path, overwrite=True)
    # TODO: Consider whether we want to write out temp information on chisq values for each potential break combination
    #   like in single breaks

    for section in section_group:
        success_tsv_path = f"{output_tsv_path}/{section}.tsv"
        with hl.hadoop_open(success_tsv_path, "w") as o:
            o.write(f"{section}\n")


def main(args):
    """Search for two simultaneous breaks in transcripts without evidence of a single significant break."""
    hl.init(log="search_for_two_breaks_run_batches.log", tmp_dir=TEMP_PATH_WITH_DEL)

    # Make sure custom machine wasn't specified with under threshold
    if args.under_threshold and args.use_custom_machine:
        raise DataException(
            "Do not specify --use-custom-machine when transcripts are --under-threshold size!"
        )

    logger.info("Importing SetExpression with transcripts...")
    transcripts_to_run = check_for_successful_transcripts(
        transcripts=(
            list(
                hl.eval(
                    hl.experimental.read_expression(
                        sections_to_simul_by_threshold_path(
                            is_rescue=args.is_rescue,
                            search_num=args.search_num,
                            is_over_threshold=False,
                        )
                    )
                )
            )
            if args.under_threshold
            else list(
                hl.eval(
                    hl.experimental.read_expression(
                        sections_to_simul_by_threshold_path(
                            is_rescue=args.is_rescue,
                            search_num=args.search_num,
                            is_over_threshold=True,
                        )
                    )
                )
            )
        ),
        is_rescue=args.is_rescue,
        search_num=args.search_num,
    )
    logger.info("Found %i transcripts to search...", len(transcripts_to_run))

    if not args.docker_image:
        logger.info("Picking default docker image...")
        # Use a docker image that specifies spark memory allocation if --use-custom-machine was specified
        if args.use_custom_machine:
            args.docker_image = "gcr.io/broad-mpg-gnomad/rmc:20220304"
        # Otherwise, use the default docker image
        else:
            args.docker_image = "gcr.io/broad-mpg-gnomad/tgg-methods-vm:20220302"

    logger.info("Setting up Batch parameters...")
    backend = hb.ServiceBackend(
        billing_project=args.billing_project,
        remote_tmpdir=args.batch_bucket,
        google_project=args.google_project,
    )
    b = hb.Batch(
        name="simul_breaks",
        backend=backend,
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
            job_name = f'group{count}{"over" if args.over_threshold else "under"}'
            j = b.new_python_job(name=job_name)
            j.memory(args.batch_memory)
            j.cpu(args.batch_cpu)
            j.storage(args.batch_storage)
            j.call(
                process_section_group,
                not_one_break_grouped.path,
                group,
                args.over_threshold,
                f"{SIMUL_BREAK_TEMP_PATH}/hts/round{args.search_num}/simul_break_{job_name}.ht",
                f"{SIMUL_BREAK_TEMP_PATH}/success_files",
                None,
                args.chisq_threshold,
                split_window_size,
            )
            count += 1

    else:
        transcript_groups = [[transcript] for transcript in transcripts_to_run]
        for group in transcript_groups:
            if args.use_custom_machine:
                # NOTE: you do not specify memory and cpu when specifying a custom machine
                j = b.new_python_job(name=job_name)
                j._machine_type = "n1-highmem-32"
                j._preemptible = True
                j.storage("100Gi")
            else:
                j = b.new_python_job(name=group[0])
                j.memory(args.batch_memory)
                j.cpu(args.batch_cpu)
                j.storage(args.batch_storage)
            j.call(
                process_section_group,
                not_one_break_grouped.path,
                group,
                args.over_threshold,
                f"{SIMUL_BREAK_TEMP_PATH}/hts/round{args.search_num}/simul_break_{group[0]}.ht",
                f"{SIMUL_BREAK_TEMP_PATH}/success_files",
                f"{SIMUL_BREAK_TEMP_PATH}",
                args.chisq_threshold,
                args.group_size,
            )

    b.run(wait=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        This regional missense constraint script searches for two simultaneous breaks in transcripts without evidence
        of a single significant break.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--chisq-threshold",
        help="Chi-square significance threshold. Value should be 9.2 (value adjusted from ExAC code due to discussion with Mark).",
        type=float,
        default=9.2,
    )
    parser.add_argument(
        "--slack-channel",
        help="Send message to Slack channel/user.",
    )
    parser.add_argument(
        "--search-num",
        help="Search iteration number (e.g., second round of searching for two simultaneous breaks would be 2).",
        type=int,
    )
    parser.add_argument(
        "--is-rescue",
        help="""
        Whether search is part of the 'rescue' pathway (pathway
        with lower chi square significance cutoff).
        """,
        action="store_true",
    )

    transcript_size = parser.add_mutually_exclusive_group(required=True)
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
    parser.add_argument(
        "--group-size",
        help="""
        Number of transcripts to include in each group of transcripts to be submitted to Hail Batch if --under-threshold.
        Size of windows to split transcripts if --over-threshold.
        Default is 100.
        """,
        type=int,
        default=100,
    )
    parser.add_argument(
        "--billing-project",
        help="Billing project to use with hail batch.",
        default="gnomad-production",
    )
    parser.add_argument(
        "--batch-bucket",
        help="Bucket provided to hail batch for temporary storage.",
        default="gs://gnomad-tmp/kc/",
    )
    parser.add_argument(
        "--google-project",
        help="Google cloud project provided to hail batch for storage objects access.",
        default="broad-mpg-gnomad",
    )
    parser.add_argument(
        "--use-custom-machine",
        help="""
        Use custom large machine for hail batch rather than setting batch memory, cpu, and storage.
        Used when the batch defaults are not enough and jobs run into out of memory errors.
        Only necessary if --over-threshold.
        """,
        action="store_true",
    )
    parser.add_argument(
        "--batch-memory",
        help="Amount of memory to request for hail batch jobs.",
        default="standard",
    )
    parser.add_argument(
        "--batch-cpu",
        help="Number of CPUs to request for hail batch jobs.",
        default=8,
        type=int,
    )
    parser.add_argument(
        "--batch-storage",
        help="Amount of disk storage to request for hail batch jobs.",
        default="10Gi",
    )
    parser.add_argument(
        "--docker-image",
        help="""
        Docker image to provide to hail Batch. Must have dill, hail, and python installed.
        Suggested image: gcr.io/broad-mpg-gnomad/tgg-methods-vm:20220302.

        If running with --use-custom-machine, Docker image must also contain this line:
        `ENV PYSPARK_SUBMIT_ARGS="--driver-memory 8g --executor-memory 8g pyspark-shell"`
        to make sure the job allocates memory correctly.
        Suggested image: gcr.io/broad-mpg-gnomad/rmc:20220304.

        Default is None -- script will select default image if not specified on the command line.
        """,
        default=None,
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
