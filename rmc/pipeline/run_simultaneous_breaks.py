import argparse
import logging

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists, parallel_file_exists
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import (
    LOGGING_PATH,
    not_one_break,
    not_one_break_grouped,
    simul_break_over_5k,
    simul_break_temp,
    simul_break_under_5k,
)
from rmc.resources.grch37.reference_data import gene_model
from rmc.slack_creds import slack_token
from rmc.utils.constraint import group_not_one_break_ht, search_for_two_breaks
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
            logger.warning("This step should be run in Dataproc!")
            hl.init(log="/search_for_two_breaks_prep_batches.log")

            if not file_exists(not_one_break_grouped.path) or args.create_grouped_ht:
                if args.min_num_obs == 0:
                    # Make sure user didn't specify a min obs of 0
                    raise DataException(
                        "Minimum number of observed variants must be greater than zero!"
                    )

                # Get number of base pairs needed to observe `num` number of missense variants (on average)
                # This number is used to determine the min_window_size - which is the smallest allowed distance between simultaneous breaks
                ht = not_one_break.ht()
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
                group_not_one_break_ht(
                    ht=ht,
                    transcript_ht=gene_model.ht(),
                    min_window_size=min_window_size,
                )

            ht = not_one_break_grouped.ht()
            logger.info(
                "Annotating HT with length of cumulative observed list annotation..."
            )
            # This length is the number of positions with possible missense variants that need to be searched
            # Not using transcript size here because transcript size
            # doesn't necessarily reflect the number of positions that need to be searched
            ht = ht.annotate(list_len=hl.len(ht.cum_obs))

            logger.info(
                "Splitting transcripts into two categories: list length < 5k and list length >= 5k..."
            )
            under_5k = ht.aggregate(
                hl.agg.filter(ht.list_len < 5000, hl.agg.collect_as_set(ht.transcript))
            )
            over_5k = ht.aggregate(
                hl.agg.filter(ht.list_len >= 5000, hl.agg.collect_as_set(ht.transcript))
            )
            hl.experimental.write_expression(under_5k, simul_break_under_5k)
            hl.experimental.write_expression(over_5k, simul_break_over_5k)

        if args.command == "run-batches":
            logger.warning("This step should be run locally!")
            hl.init(log="search_for_two_breaks_run_batches.log")

            logger.info("Importing SetExpression with transcripts...")
            if not args.under_5k and not args.over_5k:
                raise DataException(
                    "Must specify if transcript sizes are --under-5k or --over-5k!"
                )
            transcripts = (
                list(hl.eval(hl.experimental.read_expression(simul_break_under_5k)))
                if args.under_5k
                else list(hl.eval(hl.experimental.read_expression(simul_break_over_5k)))
            )

            logger.info("Checking if any transcripts have already been searched...")
            success_file_path = f"{simul_break_temp}/success_files"
            transcript_success_map = {}
            for transcript in transcripts:
                transcript_success_map[
                    transcript
                ] = f"{success_file_path}/{transcript}_success.txt"
            success_tsvs_exist = parallel_file_exists(
                list(transcript_success_map.values())
            )

            transcripts_to_run = []
            for transcript in transcripts:
                if not success_tsvs_exist[transcript_success_map[transcript]]:
                    transcripts_to_run.append(transcript)
            logger.info("Found %i transcripts to search...", len(transcripts_to_run))

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
        "--under-5k",
        help="Transcripts in batch should have <5,000 possible missense positions.",
        action="store_true",
    )
    transcript_size.add_argument(
        "--over-5k",
        help="Transcripts in batch should have >=5,000 possible missense positions.",
        action="store_true",
    )
    run_batches.add_argument(
        "--group-size",
        help="Number of transcripts to include in each group of transcripts to be submitted to Hail Batch. Default is 100.",
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
        "--batch-memory",
        help="Amount of memory to request for hail batch jobs.",
        default=15,
        type=int,
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
        default=10,
        type=int,
    )
    run_batches.add_argument(
        "--docker-image",
        help="Docker image to provide to hail batch. Must have dill, hail, and python installed.",
        default="gcr.io/broad-mpg-gnomad/tgg-methods-vm:latest",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
