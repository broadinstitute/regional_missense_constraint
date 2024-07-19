"""
This script searches for two simultaneous breaks in groups of transcripts using Hail Batch.

# TODO: Add RMC repo to docker image and make sure it's up to date
NOTE: `process_section_group` has been copied here to make sure Hail initializes properly
(and can read from requester-pays buckets) in Batch.

This script should be run locally because it submits jobs to the Hail Batch service.
"""
import argparse
import logging
from typing import List

import hail as hl
import hailtop.batch as hb
import scipy
from gnomad.utils.slack import slack_notifications
from tqdm import tqdm

from rmc.resources.basics import RMC_PREFIX, TEMP_PATH_WITH_FAST_DEL
from rmc.resources.rmc import (
    CURRENT_FREEZE,
    MIN_EXP_MIS,
    P_VALUE,
    grouped_single_no_break_ht_path,
    simul_search_round_bucket_path,
    simul_sections_split_by_len_path,
)
from rmc.slack_creds import slack_token
from rmc.utils.simultaneous_breaks import search_for_two_breaks

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("run_batches")
logger.setLevel(logging.INFO)


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
        args.docker_image = (
            "us-central1-docker.pkg.dev/broad-mpg-gnomad/images/rmc_simul_search"
        )
        # NOTE: Python version on your local machine must match the Python version in this image
        logger.warning(
            "Using %s image; please make sure Hail and Python versions in image are"
            " up to date",
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
            # NOTE: There was a Hail bug that prevented passing keyword args to a PythonJob in RMC freeze 7
            # See: https://hail.zulipchat.com/#narrow/stream/223457-Hail-Batch-support/topic/j.2Ecall.28.29.20ValueError.3A.20too.20many.20values.20to.20unpack.20.28expected.202.29
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
