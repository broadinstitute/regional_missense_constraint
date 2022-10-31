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
                "ENST00000373078_130823512_130829093",
                "ENST00000229595_119228672_119230332",
                "ENST00000292778_21983609_21984353",
                "ENST00000325404_181430494_181432221",
                "ENST00000358536_48459924_48465182",
                "ENST00000379177_24225534_24234372",
                "ENST00000396946_2945775_2977546",
                "ENST00000375571_31434453_31438211",
                "ENST00000298229_71934745_71946982",
                "ENST00000326793_194995465_195013088",
                "ENST00000398557_140953084_140954661",
                "ENST00000240328_59483171_59486827",
                "ENST00000262710_23538726_23564823",
                "ENST00000282633_51889192_51889245",
                "ENST00000294725_196295940_196342307",
                "ENST00000406360_108984754_109005977",
                "ENST00000414982_7600851_7626650",
                "ENST00000278193_27523113_27528320",
                "ENST00000261745_112464500_112486109",
                "ENST00000299138_46711243_46715383",
                "ENST00000389629_41771365_41775761",
                "ENST00000324225_117049449_117059481",
                "ENST00000263980_27429181_27480531",
                "ENST00000360537_144261437_144263348",
                "ENST00000321582_124584207_124855058",
                "ENST00000375650_31795987_31797348",
                "ENST00000218089_123094062_123164851",
                "ENST00000286827_32624145_32624172",
                "ENST00000310441_153220792_153237258",
                "ENST00000306732_111538579_111539768",
                "ENST00000406246_65423374_65430565",
                "ENST00000264344_89647106_89671689",
                "ENST00000393232_129251555_129297318",
                "ENST00000358102_149075870_149365850",
                "ENST00000222573_20370325_20418730",
                "ENST00000325404_181430272_181430493",
                "ENST00000261712_30771279_30771613",
                "ENST00000261745_112486110_112546826",
                "ENST00000585527_2289774_2290947",
                "ENST00000288699_27165545_27173219",
                "ENST00000262525_30789778_30793871",
                "ENST00000338458_56771302_57113357",
                "ENST00000379483_41767808_41768225",
                "ENST00000259605_36336393_36353189",
                "ENST00000367797_169484687_169487750",
                "ENST00000433060_94439597_94532822",
                "ENST00000395842_50609160_50616591",
                "ENST00000539097_62204809_62214976",
                "ENST00000305123_227596033_227662731",
                "ENST00000319296_30429153_30440920",
                "ENST00000341776_15497154_15522252",
                "ENST00000374685_33161365_33167016",
                "ENST00000273480_141457046_141457271",
                "ENST00000360131_27265232_27278681",
                "ENST00000362042_46125691_46136512",
                "ENST00000398238_44668035_44782128",
                "ENST00000380494_74664311_74806944",
                "ENST00000373812_29063133_29069824",
                "ENST00000359947_44064549_44072551",
                "ENST00000295888_85590704_85658342",
                "ENST00000367658_139487617_139487689",
                "ENST00000263381_15536535_15560762",
                "ENST00000373019_38446275_38456593",
                "ENST00000258428_100058870_100106497",
                "ENST00000336314_154092462_154173440",
                "ENST00000419308_22561643_22563090",
                "ENST00000394170_99104039_99132323",
                "ENST00000400485_37734526_37758446",
                "ENST00000267889_40650436_40660715",
                "ENST00000394670_30677136_30687691",
                "ENST00000376511_30568177_30571941",
                "ENST00000406403_54641444_54649651",
                "ENST00000253108_10227642_10230596",
                "ENST00000421367_101912106_101938978",
                "ENST00000263239_118572226_118579600",
                "ENST00000302472_40679600_40692169",
                "ENST00000374542_33288237_33297046",
                "ENST00000329101_77211752_77289325",
                "ENST00000287461_30613879_30616726",
                "ENST00000261858_69560977_69564556",
                "ENST00000265077_82837957_82878122",
                "ENST00000267807_56974499_57210769",
                "ENST00000329363_50969721_50971009",
                "ENST00000260970_170440850_170492941",
                "ENST00000285398_128014866_128046283",
                "ENST00000295797_169940153_169981187",
                "ENST00000299164_130318869_130342967",
                "ENST00000358933_67914095_67918406",
                "ENST00000352632_40978652_40998948",
                "ENST00000377705_6605214_6614595",
                "ENST00000581347_5891847_5895954",
                "ENST00000283254_165944032_165947239",
                "ENST00000485317_111146052_111174096",
                "ENST00000293677_10426888_10434099",
                "ENST00000284274_14690297_14699820",
                "ENST00000375440_113897410_113919399",
                "ENST00000375482_95709733_95766331",
                "ENST00000337907_8412457_8421385",
                "ENST00000417133_39751949_39755581",
                "ENST00000284957_66273946_66276451",
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
            output_ht_path = "gs://regional_missense_constraint/temp/simul_breaks/initial/round2/raw_results/simul_break_group31under.ht"
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
                count=31,
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
