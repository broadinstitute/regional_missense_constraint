import argparse
import logging
import os
import subprocess

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import (
    LOGGING_PATH,
    not_one_break,
    not_one_break_grouped,
    simul_break_over_10k,
    simul_break_temp,
    simul_break_under_10k,
)
from rmc.resources.grch37.reference_data import gene_model
from rmc.slack_creds import slack_token
from rmc.utils.constraint import search_for_two_breaks
from rmc.utils.generic import get_avg_bases_between_mis


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("search_for_two_breaks")
logger.setLevel(logging.INFO)


def main(args):
    """Search for two simultaneous breaks in transcripts without evidence of a single significant break."""
    try:
        if args.command == "create_grouped_ht":
            hl.init(log="/search_for_two_breaks_prep_batches.log")

            if not file_exists(not_one_break_grouped.path) or args.create_grouped_ht:
                logger.info("Creating grouped version of not one break HT...")
                ht = not_one_break.ht()

                if args.min_num_obs == 0:
                    # Make sure user didn't specify a min obs of 0
                    raise DataException(
                        "Minimum number of observed variants must be greater than zero!"
                    )

                logger.info(
                    "Getting start and end positions and total size for each transcript..."
                )
                transcript_ht = gene_model.ht()

                # Get number of base pairs needed to observe `num` number of missense variants (on average)
                # This number is used to determine the min_window_size - which is the smallest allowed distance between simultaneous breaks
                min_window_size = (
                    (
                        get_avg_bases_between_mis(
                            get_reference_genome(ht.locus).name,
                            args.get_total_exome_bases,
                            args.get_total_gnomad_missense,
                        )
                        * args.min_num_obs
                    )
                    if args.get_min_window_size
                    else args.min_window_size
                )
                logger.info(
                    "Minimum window size (window size needed to observe %i missense variants on average): %i",
                    args.min_num_obs,
                    min_window_size,
                )

                logger.info(
                    "Creating grouped HT with lists of cumulative observed and expected missense values..."
                )
                # Aggregating values into a struct here to force the positions and observed, expected missense values to stay sorted
                # `hl.agg.collect` does not guarantee order: https://hail.is/docs/0.2/aggregators.html#hail.expr.aggregators.collect
                # not_one_break HT is keyed by locus and transcript, so these values are sorted in the input HT
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
                group_ht = group_ht.annotate(
                    max_idx=hl.len(group_ht.values.positions) - 1
                )
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

            ht = not_one_break_grouped.ht()
            logger.info(
                "Annotating HT with length of cumulative observed list annotation..."
            )
            # This length is the number of positions with possible missense variants that need to be searched
            # Not using transcript size here because transcript size
            # doesn't necessarily reflect the number of positions that need to be searched
            ht = ht.annotate(list_len=hl.len(ht.cum_obs))

            logger.info(
                "Splitting transcripts into two categories: list length < 10k and list length >= 10k..."
            )
            under_10k = ht.aggregate(
                hl.agg.filter(ht.list_len < 10000, hl.agg.collect_as_set(ht.transcript))
            )
            over_10k = ht.aggregate(
                hl.agg.filter(
                    ht.list_len >= 10000, hl.agg.collect_as_set(ht.transcript)
                )
            )
            hl.experimental.write_expression(under_10k, simul_break_under_10k)
            hl.experimental.write_expression(over_10k, simul_break_over_10k)

        if args.command == "run-batches":
            hl.init(log="/search_for_two_breaks_run_batches.log")
            logger.info("Checking for output from any previous runs...")
            success_file_path = f"{simul_break_temp}/success_files"
            success_files = (
                subprocess.check_output(["gsutil", "ls", success_file_path])
                .decode("utf8")
                .strip()
                .split("\n")
            )
            tsvs = [os.path.split(f)[-1] for f in success_files if f.endswith(".tsv")]
            if tsvs:
                for tsv in tsvs:
                    with hl.hadoop_open(tsv) as t:
                        for line in t:
                            transcript = line.strip()
                            if transcript not in successful_transcripts:
                                successful_transcripts.append(transcript)
                successful_transcripts = set(successful_transcripts)

            if args.under_10k:
                transcripts = hl.eval(
                    hl.experimental.read_expression(simul_break_under_10k)
                )
                if tsvs:
                    transcripts_to_run = transcripts.difference(successful_transcripts)
                else:
                    transcripts_to_run = transcripts
                logger.info(
                    "Found %i transcripts to search...", len(transcripts_to_run)
                )

                group_size = args.group_size
                transcript_groups = [
                    transcripts_to_run[x : x + group_size]
                    for x in range(0, len(transcripts_to_run), group_size)
                ]
                ht = not_one_break_grouped.ht()
                temp_hts = []
                for count, group in enumerate(transcript_groups):
                    logger.info("Working on transcript group number: %i", count)
                    temp_ht = ht.filter(hl.literal(group).contains(ht.transcript))
                    temp_ht = search_for_two_breaks(temp_ht, args.chisq_threshold)
                    temp_ht = temp_ht.checkpoint(
                        f"{simul_break_temp}/hts/simul_break_group{count}.ht",
                        overwrite=args.overwrite,
                    )
                    temp_hts.append(temp_ht)
                    with hl.hadoop_open(
                        f"{success_file_path}/group{count}.tsv", "w"
                    ) as o:
                        o.write("\n".split(group) + "\n")

                logger.info("Joining and writing...")
                ht = temp_hts[0].union(*temp_hts[1:])
                ht.write(
                    f"{simul_break_temp}/under_10k/under_10k.ht",
                    overwrite=args.overwrite,
                )

            elif args.over_10k:
                transcripts = hl.eval(
                    hl.experimental.read_expression(simul_break_over_10k)
                )
            else:
                raise DataException(
                    "Must specify if transcript sizes are --under-10k or --over-10k!"
                )

    finally:
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

    # Create subparsers for each step
    # Need to specify `dest` to be able to check which subparser is being invoked
    # `dest`: https://docs.python.org/3/library/argparse.html#dest
    subparsers = parser.add_subparsers(title="command", dest="command")

    create_grouped_ht = subparsers.add_parser(
        "create-grouped-ht",
        help="Create hail Table grouped by transcript with cumulative observed and expected missense values collected into lists.",
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
        default=10,
        type=int,
    )

    run_batches = subparsers.add_parser(
        "run-batches", help="Run batches of transcripts using Hail Batch."
    )
    transcript_size = run_batches.add_mutually_exclusive_group()
    transcript_size.add_argument(
        "--under-10k",
        help="Transcripts in batch should have <10,000 possible missense positions.",
        action="store_true",
    )
    transcript_size.add_argument(
        "--over-10k",
        help="Transcripts in batch should have >=10,000 possible missense positions.",
        action="store_true",
    )
    run_batches.add_argument(
        "--group-size",
        help="Number of transcripts to include in each group of transcripts to be submitted to Hail Batch. Default is 100.",
        type=int,
        default=100,
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
