import argparse
import logging
import numpy as np
import pandas as pd
import sys
from tqdm import tqdm
from typing import List

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists, parallel_file_exists
from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import (
    not_one_break,
    temp_path,
)
from rmc.slack_creds import slack_token
from rmc.utils.generic import get_outlier_transcripts


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("search_for_two_breaks")
logger.setLevel(logging.INFO)


def agg_metrics(ht: hl.Table, pos_ht: hl.Table) -> None:
    """
    Group input Table by transcript and collect annotations required to search for two simultaneous breaks.

    Annotations are:
        - cumulative observed missense variant counts
        - cumulative expected missense variant counts
        - total observed missense variant count
        - total expected missense variant count
        - transcript start position
        - transcript end position
        - list of positions per transcript

    :param hl.Table ht: Input Table with potential two breaks transcripts
        (Table with transcripts that didn't have one significant break).
    :param hl.Table pos_ht: Table that has all positions per transcript..
    :return: None; writes GroupedTable to temporary path.
    """
    group_ht = ht.group_by("transcript").aggregate(
        cum_obs=hl.agg.collect(ht.cumulative_obs),
        cum_exp=hl.agg.collect(ht.cumulative_exp),
        total_obs=hl.agg.take(ht.total_obs, 1)[0],
        total_exp=hl.agg.take(ht.total_exp, 1)[0],
    )
    group_ht = group_ht.annotate(positions=pos_ht[group_ht.transcript].positions,)
    group_ht.write(f"{temp_path}/no_break_grouped.ht", overwrite=True)


def get_chunks(
    anchor_position: int,
    start_position: int,
    end_position: int,
    position_list: List[int],
    min_window_size: int,
) -> List[str]:
    """
    Create list of transcript subsections.

    :param int anchor_position: Transcript subsection start position.
    :param int start_position: Start position for transcript.
    :param int end_position: End position for transcript.
    :param List[int] position_list: List of positions in transcript.
    :param int min_window_size: Minimum window size to search for two simultaneous breaks.
    :return: List of possible transcript subsections starting at `anchor_position`.
    """
    return [
        f"{anchor_position}-{x}"
        for x in position_list
        if (
            x - anchor_position >= min_window_size
            and (anchor_position != start_position or x != end_position)
        )
    ]


def search_for_two_breaks(
    ht: hl.GroupedTable,
    transcript: str,
    success_ht_path: str,
    success_tsv_path: str,
    min_window_size: int,
    chisq_threshold: float = 13.8,
) -> None:
    """
    Search for two simultaneous breaks in a transcript.

    :param hl.GroupedTable ht: Input GroupedTable (grouped by transcript).
        Contains all annotations required to search for two simultaneous breaks.
    :param str transcript: Transcript to earch.
    :param str success_ht_path: Path to output HT (if transcript has evidence of two simultaneous breaks).
    :param str success_tsv_path: Path to output TSV (to indicate transcript has successfully been searched).
    :param int min_window_size: Minimum window size to search for two simultaneous breaks.
    :param float chisq_threshold: Chi-square significance threshold. Default is 13.8.
        Default is from ExAC RMC code and corresponds to a p-value of 0.999 with 2 degrees of freedom.
        (https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm)
    :return: None
    """
    # Filter HT to transcript
    ht = ht.filter(ht.transcript == transcript)

    # Grab necessary annotations from HT
    # Grab list of positions and first/last position seen in HT
    positions = ht.aggregate(hl.agg.take(ht.positions, 1)[0])
    start_pos = positions[0]
    end_pos = positions[-1]
    transcript_len = end_pos - start_pos
    cum_obs = ht.aggregate(hl.agg.take(ht.cum_obs, 1)[0])
    cum_exp = ht.aggregate(hl.agg.take(ht.cum_exp, 1)[0])
    total_obs = ht.aggregate(hl.agg.take(ht.total_obs, 1)[0])
    total_exp = ht.aggregate(hl.agg.take(ht.total_exp, 1)[0])
    total_oe = min(total_obs / total_exp, 1)

    # Double check positions are sorted
    sorted_positions = sorted(positions)
    for count, pos in enumerate(positions):
        for c2, p2 in enumerate(sorted_positions):
            if p2 == pos:
                if count != c2:
                    sys.exit(
                        f"{pos} is index {count} in positions but index {c2} in sorted positions"
                    )

    # Cycle through all positions in transcript and create all possible subsections starting at each position
    # These subsections used to create the two break windows
    for count, pos in enumerate(positions):
        if count == 0:
            sections = get_chunks(
                pos, start_pos, end_pos, positions[count + 1 :], min_window_size
            )
        else:
            sections.extend(
                get_chunks(
                    pos, start_pos, end_pos, positions[count + 1 :], min_window_size
                )
            )

    # Iterate over all possible transcript subsections and divide transcripts into two or three pieces
    # Store each piece as a group, where the group = (left piece, middle piece, right piece)
    groups = []
    for chunk in sections:
        # Check if section is most of transcript
        start, end = map(int, chunk.split("-"))
        length = end - start
        if transcript_len - length < min_window_size:
            if end == end_pos:
                # Create transcript subsection that is to the "left" of the current subsection
                # This subsection should be from the first position in seen in the transcript to
                # subsection start
                left = f"{start_pos}-{start}"
                # Add None for third section because this transcript subsection only
                # divides the transcript into two pieces
                groups.append((chunk, left, None))
            else:
                # Create transcript subsection that is to the "right" of the current subsection
                # This subsection should be the subsection end pos to the last position seen in the transcript
                right = f"{end}-{end_pos}"
                # Add None for third section because this transcript subsection only
                # divides the transcript into two pieces
                groups.append((chunk, right, None))
        else:
            # Create boolean for whether this subsection has a pair subsection in `sections`
            # "pair" subsection means that section plus this additional section equals the whole transcript
            found = False
            for chunk2 in sections:
                if chunk2 != chunk:
                    start2, end2 = map(int, chunk2.split("-"))
                    if start2 >= end and end2 >= end:
                        groups.append((chunk, chunk2, None))
                        found = True
            # If the subsection doesn't have any matching pair subsections in `sections`,
            # create the matching "left" and "right" sections
            if not found:
                right = None
                left = None
                if start != start_pos:
                    left = f"{start_pos}-{start}"
                if end != end_pos:
                    right = f"{end}-{end_pos}"
                # Add all three transcript subsections to list
                groups.append((chunk, left, right))

    # Remove any redundant groups (multiples of the same matches)
    # i.e., group 1 was found twice: once when subsection 1 matched with subsection 5
    # and again when subsection 5 was matched with subsection 1
    final_sections = []
    seen = []
    for values in groups:
        if set(values) not in seen:
            seen.append(set(values))
            final_sections.append(values)

    # Create dict that has this structure:
    # subsection: (section_obs, section_exp, section_oe, null_pois, alt_pois)
    value_map = {}
    for group in final_sections:
        for section in group:
            if not section:
                continue
            start, end = map(int, section.split("-"))
            start_idx = positions.index(start)
            end_idx = positions.index(end)
            if start_idx == 0:
                section_obs = cum_obs[end_idx]
                section_exp = cum_exp[end_idx]
            else:
                section_obs = cum_obs[end_idx] - cum_obs[start_idx]
                section_exp = cum_exp[end_idx] - cum_exp[start_idx]
            section_oe = min(section_obs / section_exp, 1)
            value_map[section] = (
                section_obs,
                section_exp,
                section_oe,
                hl.eval(hl.dpois(section_obs, section_exp * total_oe)),
                hl.eval(hl.dpois(section_obs, section_exp * section_oe)),
            )

    # Convert dict to hail literal
    value_map = hl.literal(value_map)

    # Convert list of transcript subsection groups into numpy ndarray
    # (this is for easy conversion into hail Table format)
    data = np.array(final_sections)
    df = pd.DataFrame.from_records(data)
    ht = hl.Table.from_pandas(df)

    # Multiply all null and alt distributions for each transcript subsection group
    ht = ht.annotate(
        null=hl.if_else(
            hl.is_defined(ht["2"]),
            value_map[ht["0"]][3] * value_map[ht["1"]][3] * value_map[ht["2"]][3],
            value_map[ht["0"]][3] * value_map[ht["1"]][3],
        ),
        alt=hl.if_else(
            hl.is_defined(ht["2"]),
            value_map[ht["0"]][4] * value_map[ht["1"]][4] * value_map[ht["2"]][4],
            value_map[ht["0"]][4] * value_map[ht["1"]][4],
        ),
    )

    # Add chi square annotation
    ht = ht.annotate(chisq=2 * (hl.log10(ht.alt) - hl.log10(ht.null)))

    # Get maximum chi square value and save transcript as HT if
    # it has two significant simultaneous breaks
    max_chisq = ht.aggregate(hl.agg.max(ht.chisq))
    max_chisq_ht = ht.filter(ht.chisq == max_chisq & ht.chisq >= chisq_threshold)
    if max_chisq_ht.count() != 0:
        max_chisq_ht.write(success_ht_path, overwrite=True)

    # Create success file to track that this transcript has been processed
    with hl.hadoop_open(success_tsv_path, "w") as o:
        o.write("")


def main(args):
    """Search for two simultaneous breaks in transcripts without evidence of a single significant break."""
    transcript_tsv_path = args.transcript_tsv

    if args.get_no_break_transcripts:
        hl.init(log="/RMC_simul_breaks_get_transcripts.log")
        logger.info("Reading in not one break HT...")
        context_ht = not_one_break.ht()
        context_ht = context_ht.drop("_obs_scan")

        if args.remove_outlier_transcripts:
            outlier_transcripts = get_outlier_transcripts()
            context_ht = context_ht.filter(
                ~outlier_transcripts.contains(context_ht.transcript)
            )

        # Converting to list here because `hl.agg.collect_as_set` will return `frozenset`
        transcripts = list(
            context_ht.aggregate(hl.agg.collect_as_set(context_ht.transcript))
        )
        with hl.hadoop_open(transcript_tsv_path, "w") as o:
            for transcript in transcripts:
                o.write(f"{transcript}\n")

    if args.aggregate_no_break_ht:
        ht = not_one_break.ht()
        # Table stored in temp bucket that is used to calculate constraint but can be deleted afterwards
        # This Table is grouped by transcript and has all positions in each transcript collected into a list
        pos_ht_path = f"{temp_path}/pos_per_transcript.ht"
        if not file_exists(pos_ht_path) or args.overwrite_pos_ht:
            ht = not_one_break.ht()
            pos_ht = ht.group_by("transcript").aggregate(
                positions=hl.sorted(hl.agg.collect(ht.locus.position)),
            )
            pos_ht.write(f"{temp_path}/pos_per_transcript.ht", overwrite=True)
        pos_ht = hl.read_table(f"{temp_path}/pos_per_transcript.ht")
        agg_metrics(ht, pos_ht)

    if args.run_simul_breaks_batch:
        import hailtop.batch as hb

        if not file_exists(transcript_tsv_path):
            raise DataException(
                f"{transcript_tsv_path} doesn't exist. Please rerun with --get-no-break-transcripts!"
            )
        transcripts = []
        with hl.hadoop_open(transcript_tsv_path) as i:
            for line in i:
                transcripts.append(line.strip())

        logger.info("Setting up batch parameters...")
        backend = hb.ServiceBackend(
            billing_project=args.billing_project,
            remote_tmpdir=args.batch_bucket,
            google_project=args.google_project,
        )
        b = hb.Batch(
            name="simul_breaks",
            backend=backend,
            default_memory=args.batch_memory,
            default_cpu=args.batch_cpu,
            default_storage=args.batch_storage,
            default_python_image=args.docker_image,
        )

        logger.info("Checking for output file existence...")
        transcript_success_map = {}
        transcript_ht_map = {}
        for transcript in transcripts:
            transcript_success_map[transcript] = f"{temp_path}/{transcript}_success.txt"
            transcript_ht_map[transcript] = f"{temp_path}/simul_breaks_{transcript}.ht"
        success_files_exist = parallel_file_exists(
            list(transcript_success_map.values())
        )
        hts_exist = parallel_file_exists(list(transcript_ht_map.values()))

        simul_break_transcripts = []
        transcripts_to_run = []
        for transcript in transcripts:
            output_tsv = f"{temp_path}/{transcript}_success.txt"
            # Adding this code for one time use to skip transcripts that finished while running serially
            if hts_exist[transcript_ht_map[transcript]]:
                continue
            if not success_files_exist[transcript_success_map[transcript]]:
                transcripts_to_run.append(transcript)
        logger.info("Found %i transcripts to search...", len(transcripts))

        for transcript in tqdm(transcripts, unit="transcripts"):
            logger.info("Working on %s...", transcript)
            j = b.new_python_job(name=transcript)
            j.call(
                search_for_two_breaks,
                f"{temp_path}/no_break_grouped.ht",
                transcript,
                f"{temp_path}/simul_breaks_{transcript}.ht",
                transcript_success_map[transcript],
                args.min_window_size,
                args.chisq_threshold,
            )

        b.run(wait=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script searches for regional missense constraint in gnomAD"
    )
    parser.add_argument(
        "--transcript-tsv",
        help="Path to store transcripts to search for two simultaneous breaks. Path should be to a file in Google cloud storage.",
        default=f"{temp_path}/no_break_transcripts.tsv",
    )
    parser.add_argument(
        "--get-no-break-transcripts",
        help="Get all transcripts that don't have one significant break.",
        action="store_true",
    )
    parser.add_argument(
        "--remove-outlier-transcripts",
        help="Remove outlier transcripts (transcripts with too many/few LoF, synonymous, or missense variants)",
        action="store_true",
    )
    parser.add_argument(
        "--aggregate-no-break-ht",
        help="Group HT with transcripts to search by transcript and aggregate annotations required to search for two breaks.",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite-pos-ht",
        help="Overwrite the positions per transcript HT (HT keyed by transcript with a list of positiosn per transcript), even if it already exists.",
        action="store_true",
    )
    parser.add_argument(
        "--run-simul-breaks-batch",
        help="Submit hail batch job for each transcript to search.",
        action="store_true",
    )
    parser.add_argument(
        "--chisq-threshold",
        help="Chi-square significance threshold. Value should be 10.8 (single break) and 13.8 (two breaks) (values from ExAC RMC code).",
        type=float,
        default=10.8,
    )
    parser.add_argument(
        "--min-window-size",
        help="Smallest possible window size for simultaneous breaks. Determined by running --get-min-window-size.",
        type=int,
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
        "--batch-memory",
        help="Amount of memory to request for hail batch jobs.",
        default=15,
        type=int,
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
        default=10,
        type=int,
    )
    parser.add_argument(
        "--docker-image",
        help="Docker image to provide to hail batch. Must have dill, hail, and python installed.",
        default="gcr.io/broad-mpg-gnomad/tgg-methods-vm:20211130",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
