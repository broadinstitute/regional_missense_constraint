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
                "ENST00000268864_34068082_34068224",
                "ENST00000496391_102040054_102067123",
                "ENST00000382438_13543829_13543895",
                "ENST00000221462_45594654_45648838",
                "ENST00000383791_15296360_15300407",
                "ENST00000264183_235294949_235409762",
                "ENST00000356421_100210409_100759201",
                "ENST00000263370_41235128_41246765",
                "ENST00000342694_35801123_35809729",
                "ENST00000531206_45247287_45266788",
                "ENST00000556492_57858311_57882635",
                "ENST00000293677_10434100_10443978",
                "ENST00000333188_139074442_139078149",
                "ENST00000395749_44256749_44279233",
                "ENST00000367463_149262461_149398126",
                "ENST00000331738_122989190_123003479",
                "ENST00000400930_1191482_1209265",
                "ENST00000382368_5541404_5541484",
                "ENST00000374080_70338406_70356352",
                "ENST00000248668_39804701_39805976",
                "ENST00000233535_27476552_27481141",
                "ENST00000300131_57485130_57489259",
                "ENST00000338056_87329765_87407176",
                "ENST00000357398_166231424_166245748",
                "ENST00000423656_51433298_51455640",
                "ENST00000330315_138664598_138665410",
                "ENST00000530339_140186659_140188134",
                "ENST00000389194_30365139_30365270",
                "ENST00000218075_14627177_14708862",
                "ENST00000302904_42282401_42285127",
                "ENST00000334512_81065908_81076276",
                "ENST00000518783_153077690_153193429",
                "ENST00000547517_40282456_40307035",
                "ENST00000268712_15932471_16012162",
                "ENST00000359013_30647994_30713784",
                "ENST00000344922_222800943_222841354",
                "ENST00000216297_21819631_21831442",
                "ENST00000263265_49340354_49362317",
                "ENST00000264010_67596310_67645369",
                "ENST00000370685_84662411_84704181",
                "ENST00000218652_80107471_80130210",
                "ENST00000233078_1407568_1432576",
                "ENST00000228437_108137049_108155049",
                "ENST00000291900_131497668_131534693",
                "ENST00000373644_70320413_70406525",
                "ENST00000332958_117239241_117253326",
                "ENST00000367204_204042243_204092830",
                "ENST00000264824_13211475_13213975",
                "ENST00000398246_12594255_12613582",
                "ENST00000366810_226413246_226414755",
                "ENST00000296233_126061478_126070754",
                "ENST00000367178_155054459_155141387",
                "ENST00000441003_1092686_1093153",
                "ENST00000367387_197896791_197898170",
                "ENST00000394166_96869167_96875568",
                "ENST00000406925_44092964_44095241",
                "ENST00000307063_159488808_159520749",
                "ENST00000278840_61630801_61634826",
                "ENST00000323686_51428731_51429748",
                "ENST00000407418_39257455_39262774",
                "ENST00000275493_55269027_55324313",
                "ENST00000332822_54671060_54671973",
                "ENST00000370025_109242109_109244425",
                "ENST00000376148_31515971_31515995",
                "ENST00000361952_33289314_33290205",
                "ENST00000374811_64749131_64754655",
                "ENST00000309042_57774075_57796260",
                "ENST00000382276_841690_842051",
                "ENST00000551150_120636479_120639038",
                "ENST00000298282_125281728_125299894",
                "ENST00000335255_80541985_80602538",
                "ENST00000568956_56126877_56128635",
                "ENST00000310992_130417984_130418133",
                "ENST00000267294_100615218_100622611",
                "ENST00000371610_49366565_49373332",
                "ENST00000585527_2291557_2308156",
                "ENST00000535094_113548692_113739209",
                "ENST00000309035_126692062_126849739",
                "ENST00000257290_55095264_55151611",
                "ENST00000315987_112313284_112323398",
                "ENST00000327337_50746315_50790405",
                "ENST00000307677_30309786_30311792",
                "ENST00000355630_40820272_41032706",
                "ENST00000251453_39924396_39926588",
                "ENST00000425394_56712123_56713689",
                "ENST00000309660_12836837_12859192",
                "ENST00000555818_72055941_72138143",
                "ENST00000375370_114288320_114295785",
                "ENST00000303077_45235958_45236569",
                "ENST00000358365_30528036_30535283",
                "ENST00000378069_43698155_43741693",
                "ENST00000349499_36041688_36042418",
                "ENST00000409544_71503691_71591335",
                "ENST00000262584_146015150_146016793",
                "ENST00000262126_9258759_9285983",
                "ENST00000325455_100996888_101001255",
                "ENST00000261205_79257773_79685907",
                "ENST00000278174_13409548_13441082",
                "ENST00000244020_42088758_42092245",
                "ENST00000335749_76000249_76063979",
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
            output_ht_path = "gs://regional_missense_constraint/temp/simul_breaks/initial/round2/raw_results/simul_break_group40under.ht"
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
                count=40,
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
