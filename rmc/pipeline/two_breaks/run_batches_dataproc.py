"""
This script searches for two simultaneous breaks in groups of transcripts using Hail Query within Google Cloud Dataproc.

This script should be run only on transcripts that are greater than or equal
to the --transcript-len-threshold specified in `prepare_transcripts.py` if these transcripts are too slow or getting preempted in Hail Batch.
Transcripts smaller than --transcript-len-threshold  should be run using `run_batches.py` as they run quickly and inexpensively in Hail Batch.

If using this step to run TTN, use a large autoscaling cluster (highmem-8, scales to 100 preemptibles).
Otherwise, an autoscaling cluster of highmem-8s that scales to 50 preemptibles should suffice.
"""
import argparse
import logging

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import LOGGING_PATH, TEMP_PATH_WITH_DEL
from rmc.resources.rmc import (
    grouped_single_no_break_ht_path,
    simul_search_round_bucket_path,
    simul_sections_split_by_len_path,
)
from rmc.slack_creds import slack_token
from rmc.utils.simultaneous_breaks import process_section_group

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("run_batches_dataproc")
logger.setLevel(logging.INFO)


def main(args):
    """Search for two simultaneous breaks in transcripts without evidence of a single significant break."""
    try:
        logger.warning("This step should be run on an autoscaling cluster!")
        hl.init(
            log=f"/round{args.search_num}search_for_two_breaks_run_batches_dataproc.log",
            tmp_dir=TEMP_PATH_WITH_DEL,
        )
        save_chisq_ht = False
        if args.search_num == 1 and not args.is_rescue:
            save_chisq_ht = True

        if args.run_ttn:
            section_groups = [[args.ttn_id]]

        if args.run_sections_over_threshold:
            sections_to_run = list(
                hl.eval(
                    hl.experimental.read_expression(
                        simul_sections_split_by_len_path(
                            is_rescue=args.is_rescue,
                            search_num=args.search_num,
                            is_over_threshold=True,
                        )
                    )
                )
            )
            if args.group_size:
                logger.info(
                    "Splitting transcripts/transcript sections into groups of %i",
                    args.group_size,
                )
                section_groups = [
                    sections_to_run[x : x + args.group_size]
                    for x in range(0, len(sections_to_run), args.group_size)
                ]
            else:
                logger.info("Running transcripts/transcript sections one at a time...")
                section_groups = [[section] for section in sections_to_run]

        raw_path = simul_search_round_bucket_path(
            is_rescue=args.is_rescue,
            search_num=args.search_num,
            bucket_type="raw_results",
        )
        section_groups = [
            [
                "ENST00000262294_57059999_57157187",
                "ENST00000347433_100472733_100481731",
                "ENST00000292385_176893927_176901402",
                "ENST00000457091_6542743_6543129",
                "ENST00000374694_35928685_35930221",
                "ENST00000371873_47799469_47799717",
                "ENST00000439576_89346674_89386817",
                "ENST00000311234_51952584_52028400",
                "ENST00000362003_699537_731308",
                "ENST00000263331_113299492_113332535",
                "ENST00000307063_159520750_159520839",
                "ENST00000309955_202000828_202041410",
                "ENST00000263073_1963133_1968913",
                "ENST00000479441_50400233_50404348",
                "ENST00000301691_41833059_41833128",
                "ENST00000453960_153287024_153296771",
                "ENST00000403437_10304023_10304048",
                "ENST00000339475_57270975_57277197",
                "ENST00000404971_50847246_51259674",
                "ENST00000370187_103310517_103317078",
                "ENST00000358752_94278375_94283063",
                "ENST00000360565_50155961_50161899",
                "ENST00000335712_60688108_60780702",
                "ENST00000281701_224455788_224518089",
                "ENST00000302101_160340702_160342638",
                "ENST00000533486_82684175_82698664",
                "ENST00000287878_151253197_151329181",
                "ENST00000261799_149493400_149498403",
                "ENST00000602142_6828482_6828860",
                "ENST00000233557_27650657_27656281",
                "ENST00000360565_50145382_50155960",
                "ENST00000294016_4003388_4164092",
                "ENST00000282412_44395108_44429023",
                "ENST00000538426_63755837_63769283",
                "ENST00000346541_32232259_32237842",
                "ENST00000544216_34663409_34687543",
                "ENST00000371015_58349315_58422766",
                "ENST00000581347_5882071_5891846",
                "ENST00000265990_93757459_93790082",
                "ENST00000343537_223967601_223983578",
                "ENST00000349769_145626882_145634753",
                "ENST00000393845_113005777_113092327",
                "ENST00000338101_117645947_117723965",
                "ENST00000359655_70939988_70945757",
                "ENST00000325222_148395006_148454109",
                "ENST00000315073_180651587_180662809",
                "ENST00000308406_19281034_19285417",
                "ENST00000311457_24641062_24641949",
                "ENST00000370397_101948055_101953165",
                "ENST00000370508_101163602_101190381",
                "ENST00000504595_15936727_15939900",
                "ENST00000356805_54683422_54843406",
                "ENST00000367435_193091147_193104690",
                "ENST00000598441_49588676_49611495",
                "ENST00000399788_461490_498620",
                "ENST00000409235_38834265_38861589",
                "ENST00000261731_113905192_113906072",
                "ENST00000356840_16129356_16130895",
                "ENST00000268171_91411822_91422922",
                "ENST00000008391_50681541_50696607",
                "ENST00000392476_75202825_75213179",
                "ENST00000369076_106632351_106634491",
                "ENST00000336395_35608176_35610038",
                "ENST00000252599_17666403_17688030",
                "ENST00000262032_56420719_56432219",
                "ENST00000254337_14067011_14070291",
                "ENST00000481195_133530025_133541803",
                "ENST00000429205_7465192_7468815",
                "ENST00000346872_37949104_38020441",
                "ENST00000309657_66612376_66614017",
                "ENST00000427805_100645812_100646746",
                "ENST00000446108_6130950_6157226",
                "ENST00000339582_126173333_126173669",
                "ENST00000393892_40120569_40129749",
                "ENST00000295770_31574130_31666530",
                "ENST00000340635_58511694_59817947",
                "ENST00000324103_24615892_24621023",
                "ENST00000350997_61584182_61596790",
                "ENST00000371621_49126891_49195828",
                "ENST00000360652_2818891_2821836",
                "ENST00000327064_11031521_11033453",
                "ENST00000379491_10873456_10876184",
                "ENST00000323116_182511288_182538142",
                "ENST00000342310_165182967_165325952",
                "ENST00000378501_9389302_9402050",
                "ENST00000542735_140050379_140052207",
                "ENST00000282249_106680971_106681109",
                "ENST00000536769_125396668_125397000",
                "ENST00000374259_70316047_70316666",
                "ENST00000266126_75475828_75476292",
                "ENST00000394976_134181370_134190851",
                "ENST00000248248_77224732_77227503",
                "ENST00000367506_185087220_185119561",
                "ENST00000377807_17594323_17614213",
                "ENST00000529230_46831032_46867847",
                "ENST00000396352_27145803_27148122",
                "ENST00000295108_182537815_182543243",
                "ENST00000369144_120817735_120840316",
                "ENST00000267197_122242086_122265451",
                "ENST00000376874_124017956_124018265",
            ]
        ]
        for counter, group in enumerate(section_groups):
            # Double check TTN has been removed
            if args.ttn_id in group:
                group.remove(args.ttn_id)

            """output_ht_path = (
                f"{raw_path}/simul_break_dataproc_ttn.ht"
                if args.run_ttn
                else f"{raw_path}/simul_break_dataproc_{counter}.ht"
            )"""
            output_ht_path = "gs://regional_missense_constraint/temp/simul_breaks/initial/round2/raw_results/simul_break_group32under.ht"
            if file_exists(output_ht_path):
                raise DataException(
                    f"Output already exists at {output_ht_path}! Double check before running script again."
                )

            process_section_group(
                ht_path=grouped_single_no_break_ht_path(
                    args.is_rescue, args.search_num
                ),
                section_group=group,
                # count=counter,
                count=32,
                is_rescue=args.is_rescue,
                search_num=args.search_num,
                over_threshold=True,
                output_ht_path=output_ht_path,
                chisq_threshold=args.chisq_threshold,
                split_list_len=args.split_list_len,
                read_if_exists=args.read_if_exists,
                save_chisq_ht=save_chisq_ht,
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(LOGGING_PATH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        This regional missense constraint script searches for two simultaneous breaks in transcripts without evidence
        of a single significant break.
        """,
        # Add default values for args to help message
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
    parser.add_argument(
        "--group-size",
        help="""
        Number of transcripts/transcript sections to include in each group to be run.
        """,
        type=int,
    )
    parser.add_argument(
        "--split-list-len",
        help="Max length to divide transcript/sections observed or expected missense and position lists into.",
        type=int,
        default=500,
    )
    section_ids = parser.add_mutually_exclusive_group()
    section_ids.add_argument(
        "--run-sections-over-threshold",
        help="Search for simultaneous breaks in sections that are over length cutoff.",
        action="store_true",
    )
    section_ids.add_argument(
        "--run-ttn",
        help="Run TTN. TTN is so large that it needs to be treated separately.",
        action="store_true",
    )
    parser.add_argument(
        "--ttn-id",
        help="TTN transcript ID. TTN is so large that it needs to be treated separately.",
        default="ENST00000589042",
    )
    parser.add_argument(
        "--read-if-exists",
        help="Use temporary Tables if they already exist.",
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
