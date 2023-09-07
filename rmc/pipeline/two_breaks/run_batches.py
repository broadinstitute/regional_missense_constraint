"""
This script searches for two simultaneous breaks in groups of transcripts using Hail Batch.

# TODO: Add RMC repo to a docker image and make sure its updated

Note that a couple functions have been copied into this script from `constraint.py`:
- `get_obs_exp_expr`
- `get_dpois_expr`

Note also that a few functions have been copied into this script from `simultaneous_breaks.py`:
- `calculate_window_chisq`
- `search_for_two_breaks`
- `process_section_group`

This is because python imports do not work in Hail Batch PythonJobs unless
the python scripts are included within the provided Dockerfile, and the scripts within the RMC repo are
too interdependent (would have to copy the entire repo into the Dockerfile).

This script should be run locally because it submits jobs to the Hail Batch service.
"""
import argparse
import logging
from typing import List, Optional, Set

import hail as hl
import hailtop.batch as hb
import scipy
from gnomad.resources.resource_utils import DataException
from gnomad.utils.slack import slack_notifications
from tqdm import tqdm

from rmc.resources.basics import TEMP_PATH_WITH_FAST_DEL
from rmc.resources.rmc import (
    CURRENT_FREEZE,
    MIN_CHISQ_THRESHOLD,
    MIN_EXP_MIS,
    P_VALUE,
    grouped_single_no_break_ht_path,
    simul_sections_split_by_len_path,
)
from rmc.slack_creds import slack_token

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("run_batches")
logger.setLevel(logging.INFO)


RMC_PREFIX = "gs://regional_missense_constraint"
"""
Path to bucket attached to regional missense constraint (RMC) project.

Adding this constant here to avoid ModuleNotFound errors in the
PythonJobs. See `rmc.resources.basics` for full docstring.
"""

TEMP_PATH_WITH_FAST_DEL = "gs://gnomad-tmp-4day/rmc"
"""
Path to bucket for temporary files.

Adding this constant here to avoid ModuleNotFound errors in the
PythonJobs. See `rmc.resources.basics` for full docstring.
"""

SIMUL_BREAK_TEMP_PATH = f"{RMC_PREFIX}/temp/simul_breaks"
"""
Path to bucket to store temporary results for simultaneous searches.

Adding this constant here to avoid ModuleNotFound errors in the
PythonJobs. See `rmc.resources.basics` for full docstring.
"""


SIMUL_SEARCH_BUCKET_NAMES = {"prep", "raw_results", "final_results", "success_files"}
"""
Names of buckets nested within round bucket of `SIMUL_BREAK_TEMP_PATH`.

Adding this constant here to avoid ModuleNotFound errors in the
PythonJobs. See `rmc.resources.rmc` for full docstring.
"""


def simul_search_bucket_path(
    search_num: Optional[int] = None,
    freeze: int = CURRENT_FREEZE,
) -> str:
    """
    Return path to bucket associated with simultaneous break search inputs and results.

    Adding this constant here to avoid ModuleNotFound errors in the
    PythonJobs. See `rmc.resources.rmc` for full docstring.

    :param search_num: Search iteration number
        (e.g., second round of searching for simultaneous break would be 2).
        Default is None.
    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: Path to simultaneous break search round bucket.
    """
    return (
        f"{SIMUL_BREAK_TEMP_PATH}/{freeze}/round{search_num}"
        if search_num
        else f"{SIMUL_BREAK_TEMP_PATH}/{freeze}"
    )


def simul_search_round_bucket_path(
    search_num: int,
    bucket_type: str,
    bucket_names: Set[str] = SIMUL_SEARCH_BUCKET_NAMES,
    freeze: int = CURRENT_FREEZE,
) -> str:
    """
    Return path to bucket with  Tables resulting from a specific round of simultaneous break search.

    Adding this constant here to avoid ModuleNotFound errors in the
    PythonJobs. See `rmc.resources.rmc` for full docstring.

    :param search_num: Search iteration number
        (e.g., second round of searching for single break would be 2).
    :param bucket_type: Bucket type.
        Must be in `bucket_names`.
    :param bucket_names: Possible bucket names for simultaneous search bucket type.
        Default is `SIMUL_SEARCH_BUCKET_NAMES`.
    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: Path to a bucket in the simultaneous break search round bucket.
    """
    assert bucket_type in bucket_names, f"Bucket type must be one of {bucket_names}!"
    return f"{simul_search_bucket_path(search_num, freeze)}/{bucket_type}"


def get_obs_exp_expr(
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
) -> hl.expr.Float64Expression:
    """
    Return observed/expected annotation based on inputs.

    Typically imported from `constraint.py`. See `constraint.py` for full docstring.

    :param hl.expr.Int64Expression obs_expr: Expression containing number of observed variants.
    :param hl.expr.Float64Expression exp_expr: Expression containing number of expected variants.
    :return: Observed/expected expression.
    :rtype: hl.expr.Float64Expression
    """
    # Cap the o/e ratio at 1 to avoid pulling out regions that are enriched for missense variation
    # Code is looking for missense constraint, so regions with a ratio of >= 1.0 can be grouped together
    return hl.nanmin(obs_expr / exp_expr, 1)


def get_dpois_expr(
    oe_expr: hl.expr.Float64Expression,
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
) -> hl.expr.StructExpression:
    """
    Calculate probability densities (natural log) of the observed values under a Poisson model.

    Typically imported from `constraint.py`. See `constraint.py` for full docstring.

    :param oe_expr: Expression of observed/expected value.
    :param obs_expr: Expression containing observed variants count.
    :param exp_expr: Expression containing expected variants count.
    :return: Natural log of the probability density under Poisson model.
    """
    # log_p = True returns the natural logarithm of the probability density
    return hl.dpois(obs_expr, exp_expr * oe_expr, log_p=True)


def calculate_window_chisq(
    max_idx: hl.expr.Int32Expression,
    i: hl.expr.Int32Expression,
    j: hl.expr.Int32Expression,
    cum_obs: hl.expr.ArrayExpression,
    cum_exp: hl.expr.ArrayExpression,
    section_oe: hl.expr.Float64Expression,
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
    :param hl.expr.Float64Expression section_oe: Transcript/transcript section overall observed/expected (OE) missense ratio.
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
                        oe_expr=get_obs_exp_expr(
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
                        oe_expr=get_obs_exp_expr(
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
                        oe_expr=get_obs_exp_expr(
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
                        oe_expr=section_oe,
                        obs_expr=cum_obs[i - 1],
                        exp_expr=hl.max(cum_exp[i - 1], 1e-09),
                    )
                    + get_dpois_expr(
                        oe_expr=section_oe,
                        obs_expr=cum_obs[j] - cum_obs[i - 1],
                        exp_expr=hl.max(cum_exp[j] - cum_exp[i - 1], 1e-09),
                    )
                    + get_dpois_expr(
                        oe_expr=section_oe,
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
    chisq_threshold: float = scipy.stats.chi2.ppf(1 - P_VALUE, 2),
    min_num_exp_mis: float = MIN_EXP_MIS,
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
    :param count: Which transcript group is being run (based on counter generated in `main`).
    :param chisq_threshold: Chi-square significance threshold. Default is
        `scipy.stats.chi2.ppf(1 - P_VALUE, 2)`.
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
                        group_ht.section_oe,
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
        chisq=group_ht.best_break.chisq,
        # Adjust breakpoint inclusive/exclusiveness to be consistent with single break breakpoints, i.e.
        # so that the breakpoint site itself is the last site in the left subsection. Thus, the resulting
        # subsections will be divided as follows:
        # [section_start, first_break_pos] (first_break_pos, second_break_pos] (second_break_pos, section_end]
        breakpoints=(
            group_ht.positions[group_ht.best_break.i] - 1,
            group_ht.positions[group_ht.best_break.j],
        ),
    )
    group_ht_path = (
        f"{TEMP_PATH_WITH_FAST_DEL}/freeze{freeze}_batch_temp_chisq_group{count}.ht"
    )
    if save_chisq_ht:
        group_ht_path = (
            f"{SIMUL_BREAK_TEMP_PATH}/freeze{freeze}/batch_temp_chisq_group{count}.ht"
        )
    group_ht = group_ht.checkpoint(group_ht_path, overwrite=True)

    # Remove rows with maximum chi square values below the threshold
    group_ht = group_ht.filter(group_ht.chisq >= chisq_threshold)
    return group_ht


def process_section_group(
    ht_path: str,
    section_group: List[str],
    count: int,
    search_num: int,
    over_threshold: bool,
    output_ht_path: str,
    output_n_partitions: int = 10,
    chisq_threshold: float = scipy.stats.chi2.ppf(1 - P_VALUE, 2),
    min_num_exp_mis: float = MIN_EXP_MIS,
    split_list_len: int = 500,
    read_if_exists: bool = False,
    save_chisq_ht: bool = False,
    requester_pays_bucket: str = RMC_PREFIX,
    google_project: str = "broad-mpg-gnomad",
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
        `scipy.stats.chi2.ppf(1 - P_VALUE, 2)`.
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
        This saves a lot of extra data and should only occur once.
        Default is False.
    :param requester_pays_bucket: Requester-pays bucket to initialize hail configuration with.
        Default is `RMC_PREFIX`.
    :param google_project: Google project used to read and write data to requester-pays bucket.
    :param freeze: RMC freeze number. Default is CURRENT_FREEZE.
    :return: None; processes Table and writes to path. Also writes success TSV to path.
    """
    # Initialize hail to read from requester-pays correctly
    hl.init(
        spark_conf={
            "spark.hadoop.fs.gs.requester.pays.mode": "CUSTOM",
            "spark.hadoop.fs.gs.requester.pays.buckets": (
                f"{requester_pays_bucket.lstrip('gs:/')}"
            ),
            "spark.hadoop.fs.gs.requester.pays.project.id": f"{google_project}",
        },
        tmp_dir=TEMP_PATH_WITH_FAST_DEL,
    )
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
            _read_if_exists=read_if_exists,
            overwrite=not read_if_exists,
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
        freeze=freeze,
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
            group_ht = ht.group_by("section").aggregate(
                section_max_chisq=hl.agg.max(ht.max_chisq)
            )
            ht = ht.annotate(section_max_chisq=group_ht[ht.section].section_max_chisq)
            ht = ht.filter(ht.chisq == ht.section_max_chisq)

    ht = ht.annotate_globals(chisq_threshold=chisq_threshold)
    ht = ht.naive_coalesce(output_n_partitions)
    ht.write(output_ht_path, overwrite=True)

    success_tsvs_path = simul_search_round_bucket_path(
        search_num=search_num,
        bucket_type="success_files",
        freeze=freeze,
    )
    for section in section_group:
        tsv_path = f"{success_tsvs_path}/{section}.tsv"
        with hl.hadoop_open(tsv_path, "w") as o:
            o.write(f"{section}\n")


def main(args):
    """Search for two simultaneous breaks in transcripts without evidence of a single significant break."""
    hl.init(
        spark_conf={
            "spark.hadoop.fs.gs.requester.pays.mode": "CUSTOM",
            "spark.hadoop.fs.gs.requester.pays.buckets": f"{RMC_PREFIX.lstrip('gs:/')}",
            "spark.hadoop.fs.gs.requester.pays.project.id": f"{args.google_project}",
        },
        log="search_for_two_breaks_run_batches.log",
        tmp_dir=TEMP_PATH_WITH_FAST_DEL,
    )

    # Make sure custom machine wasn't specified with under threshold
    if args.under_threshold and args.use_custom_machine:
        raise DataException(
            "Do not specify --use-custom-machine when transcripts/sections are"
            " --under-threshold size!"
        )

    chisq_threshold = hl.eval(hl.qchisqtail(P_VALUE, 2))
    if args.p_value:
        chisq_threshold = hl.eval(hl.qchisqtail(args.p_value, 2))

    logger.info("Importing SetExpression with transcripts or transcript sections...")
    sections_to_run = (
        list(
            hl.eval(
                hl.experimental.read_expression(
                    simul_sections_split_by_len_path(
                        search_num=args.search_num,
                        is_over_threshold=False,
                        freeze=args.freeze,
                    )
                )
            )
        )
        if args.under_threshold
        else list(
            hl.eval(
                hl.experimental.read_expression(
                    simul_sections_split_by_len_path(
                        search_num=args.search_num,
                        is_over_threshold=True,
                        freeze=args.freeze,
                    )
                )
            )
        )
    )
    logger.info(
        "Found %i transcripts or transcript sections to search...", len(sections_to_run)
    )

    raw_path = simul_search_round_bucket_path(
        search_num=args.search_num,
        bucket_type="raw_results",
        freeze=args.freeze,
    )

    if not args.docker_image:
        logger.info("Picking default docker image...")
        # Use a docker image that specifies spark memory allocation if --use-custom-machine was specified
        if args.use_custom_machine:
            args.docker_image = "gcr.io/broad-mpg-gnomad/rmc:20220930"
            raise DataException(f"Hail in {args.docker_image} image is out of date!")
        # Otherwise, use the default docker image
        else:
            args.docker_image = (
                "us-central1-docker.pkg.dev/broad-mpg-gnomad/images/rmc_simul_search"
            )
            logger.warning(
                "Using %s image; please make sure Hail version in image is up to date",
                args.docker_image,
            )

    logger.info("Setting up Batch parameters...")
    backend = hb.ServiceBackend(
        billing_project=args.billing_project,
        remote_tmpdir=args.batch_bucket,
        gcs_requester_pays_configuration=args.google_project,
    )
    b = hb.Batch(
        name="simul_breaks",
        backend=backend,
        default_python_image=args.docker_image,
        requester_pays_project=args.google_project,
    )
    # Check if user specified list of numbers for batches
    # These numbers are used to write output files for batch jobs
    count_list = None
    if args.counter:
        count_list = list(map(int, args.counter.split(",")))

    if args.under_threshold:
        section_groups = [
            sections_to_run[x : x + args.group_size]
            for x in range(0, len(sections_to_run), args.group_size)
        ]
        count = 1

        for group in tqdm(section_groups, unit="section group"):
            if count_list:
                group_num = count_list[count - 1]
                assert len(count_list) == len(
                    section_groups
                ), "Number of section groups doesn't match specified number of batches!"
            else:
                group_num = count

            logger.info("Working on group number %s...", group_num)
            logger.info(group)
            job_name = f"group{group_num}under"
            j = b.new_python_job(name=job_name)
            j.memory(args.batch_memory)
            j.cpu(args.batch_cpu)
            j.storage(args.batch_storage)
            # NOTE: This section was commented out due to a Hail bug
            # The bug made us unable to pass keyword args to a PythonJob
            # See: https://hail.zulipchat.com/#narrow/stream/223457-Hail-Batch-support/topic/j.2Ecall.28.29.20ValueError.3A.20too.20many.20values.20to.20unpack.20.28expected.202.29
            # Revert to keyword args since the above bug is fixed in 0.2.122
            j.call(
                process_section_group,
                ht_path=grouped_single_no_break_ht_path(args.search_num, args.freeze),
                section_group=group,
                count=group_num,
                search_num=args.search_num,
                over_threshold=False,
                output_ht_path=f"{raw_path}/simul_break_{job_name}.ht",
                output_n_partitions=args.output_n_partitions,
                chisq_threshold=chisq_threshold,
                split_list_len=args.group_size,
                save_chisq_ht=args.save_chisq_ht,
                google_project=args.google_project,
                freeze=args.freeze,
            )
            count += 1

    else:
        section_groups = [[section] for section in sections_to_run]
        for group in section_groups:
            if count_list:
                group_num = count_list[count - 1]
                assert len(count_list) == len(
                    section_groups
                ), "Number of section groups doesn't match specified number of batches!"
            else:
                group_num = count
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
                ht_path=grouped_single_no_break_ht_path(args.search_num, args.freeze),
                section_group=group,
                count=group_num,
                search_num=args.search_num,
                over_threshold=True,
                output_ht_path=f"{raw_path}/simul_break_{group[0]}.ht",
                output_n_partitions=args.output_n_partitions,
                chisq_threshold=chisq_threshold,
                split_list_len=args.group_size,
                save_chisq_ht=args.save_chisq_ht,
                google_project=args.google_project,
                freeze=args.freeze,
            )
            count += 1
    b.run(wait=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        This regional missense constraint script searches for two simultaneous breaks in transcripts/transcript sections without evidence
        of a single significant break.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--p-value",
        help="""
        p-value significance threshold for single break search.
        Used to determine chi square threshold for likelihood ratio test.

        If not specified, script will default to threshold set
        in `P_VALUE`.
        """,
        type=float,
    )
    parser.add_argument(
        "--output-n-partitions",
        help="Number of desired partitions for output Tables. Default is 10.",
        type=int,
        default=10,
    )
    parser.add_argument(
        "--slack-channel",
        help="Send message to Slack channel/user.",
    )
    parser.add_argument(
        "--search-num",
        help=(
            "Search iteration number (e.g., second round of searching for two"
            " simultaneous breaks would be 2)."
        ),
        type=int,
    )
    parser.add_argument(
        "--counter",
        help="""
        Comma separated string of counter numbers, e.g. '31,32,40'.
        These numbers correspond to batch numbers and should only be
        specified for any batches that failed the first submission.
        """,
        action="store_true",
    )
    parser.add_argument(
        "--freeze",
        help="RMC data freeze number",
        type=int,
        default=CURRENT_FREEZE,
    )
    parser.add_argument(
        "--save-chisq-ht",
        help="""
        Save temporary Table that contains chi square significance values
        for all possible loci. Note that chi square values will be missing for
        any loci that would divide a transcript into subsections with fewer than
        `MIN_EXP_MIS` expected missense variants.

        NOTE that this temporary Table should only get saved once per gnomAD version,
        (save when running RMC pipeline at strictest p-value threshold).
        """,
        action="store_true",
    )

    section_size = parser.add_mutually_exclusive_group(required=True)
    section_size.add_argument(
        "--under-threshold",
        help=(
            "Transcripts/sections in batch should have less than"
            " --section-len-threshold possible missense positions."
        ),
        action="store_true",
    )
    section_size.add_argument(
        "--over-threshold",
        help=(
            "Transcripts/sections in batch should have greater than or equal to"
            " --section-len-threshold possible missense positions."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--group-size",
        help="""
        Number of transcripts/transcript sections to include in each group to be submitted to Hail Batch if --under-threshold.
        Size of windows to split transcripts/sections if --over-threshold.
        Default is 100.
        """,
        type=int,
        default=100,
    )
    parser.add_argument(
        "--billing-project",
        help="Billing project to use with hail batch.",
        default="rmc-production",
    )
    parser.add_argument(
        "--batch-bucket",
        help="Bucket provided to hail batch for temporary storage.",
        default=f"{TEMP_PATH_WITH_FAST_DEL}/rmc/",
    )
    parser.add_argument(
        "--google-project",
        help="""
            Google cloud project provided to hail batch for storage objects access.
            Also used to read and write to requester-pays buckets.
            """,
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
        Suggested image: us-central1-docker.pkg.dev/broad-mpg-gnomad/images/rmc_simul_search.

        If running with --use-custom-machine, Docker image must also contain this line:
        `ENV PYSPARK_SUBMIT_ARGS="--driver-memory 8g --executor-memory 8g pyspark-shell"`
        to make sure the job allocates memory correctly.
        Example image: gcr.io/broad-mpg-gnomad/rmc:20220930.

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
