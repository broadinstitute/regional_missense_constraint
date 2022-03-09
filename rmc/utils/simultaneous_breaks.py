"""This script contains functions used to search for two simultaneous breaks."""
import logging
from typing import List, Optional

import hail as hl

from gnomad.utils.file_utils import file_exists, parallel_file_exists
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.resources.resource_utils import DataException

from rmc.resources.basics import (
    not_one_break_grouped,
    simul_break_over_threshold,
    simul_break_temp,
    simul_break_under_threshold,
)
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


def split_transcripts_by_len(
    ht: hl.Table, transcript_len_threshold: int, ttn_id: str
) -> None:
    """
    Split transcripts based on the specified number of possible missense variants.

    This is necessary because transcripts with more possible missense variants take longer to run through `hl.experimental.loop`.

    :param hl.Table ht: Input Table (Table written using `group_not_one_break_ht`).
    :param int transcript_len_threshold: Possible number of missense variants cutoff.
    :param str ttn_id: TTN transcript ID. TTN is large and needs to be processed separately.
    :return: None; writes SetExpressions to resource paths (`simul_break_under_threshold`, `simul_break_over_threshold`).
    """
    logger.info("Annotating HT with length of cumulative observed list annotation...")
    # This length is the number of positions with possible missense variants that need to be searched
    # Not using transcript size here because transcript size
    # doesn't necessarily reflect the number of positions that need to be searched
    ht = ht.annotate(list_len=hl.len(ht.cum_obs))

    logger.info(
        "Splitting transcripts into two categories: list length < %i and list length >= %i...",
        transcript_len_threshold,
        transcript_len_threshold,
    )
    under_threshold = ht.aggregate(
        hl.agg.filter(
            ht.list_len < transcript_len_threshold,
            hl.agg.collect_as_set(ht.transcript),
        )
    )
    over_threshold = ht.aggregate(
        hl.agg.filter(
            ht.list_len >= transcript_len_threshold,
            hl.agg.collect_as_set(ht.transcript),
        )
    )
    if over_threshold.contains(ttn_id):
        logger.warning(
            "TTN is present in input transcripts! It will need to be run separately."
        )
        over_threshold = over_threshold.remove(ttn_id)
    hl.experimental.write_expression(under_threshold, simul_break_under_threshold)
    hl.experimental.write_expression(over_threshold, simul_break_over_threshold)


def check_for_successful_transcripts(
    transcripts: List[str], in_parallel: bool = True
) -> List[str]:
    """
    Check if any transcripts have been previously searched by searching for success TSV existence.

    .. note::
        This step needs to be run locally due to permissions involved with `parallel_file_exists`.

    :param List[str] transcripts: List of transcripts to check.
    :param bool in_parallel: Whether to check if successful file exist in parallel.
        If True, must be run locally and not in Dataproc. Default is True.
    :return: List of transcripts didn't have success TSVs and therefore still need to be processed.
    """
    logger.info("Checking if any transcripts have already been searched...")
    success_file_path = f"{simul_break_temp}/success_files"
    transcript_success_map = {}
    transcripts_to_run = []
    for transcript in transcripts:
        transcript_success_map[
            transcript
        ] = f"{success_file_path}/{transcript}_success.tsv"

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
