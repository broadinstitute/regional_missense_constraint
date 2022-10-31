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
                "ENST00000468418_177001340_177034319",
                "ENST00000378679_132157833_132161153",
                "ENST00000518448_124219705_124222314",
                "ENST00000225235_96162261_96234425",
                "ENST00000370685_84649805_84662410",
                "ENST00000222330_42737276_42744204",
                "ENST00000430629_27741373_27741426",
                "ENST00000373036_38275239_38297975",
                "ENST00000322941_1957448_1961544",
                "ENST00000371497_51588946_51870574",
                "ENST00000396029_49505585_49508757",
                "ENST00000261667_50273447_50279789",
                "ENST00000337872_4304747_4323513",
                "ENST00000302326_28144265_28195171",
                "ENST00000396884_38374005_38383429",
                "ENST00000266022_49977440_50095299",
                "ENST00000392425_123636867_123678972",
                "ENST00000229239_6643093_6646797",
                "ENST00000202556_104206806_104313927",
                "ENST00000361923_7844380_7889936",
                "ENST00000260323_54542483_54557683",
                "ENST00000335953_114112899_114121398",
                "ENST00000372970_42825136_42825822",
                "ENST00000367858_134490384_134494633",
                "ENST00000391884_48958766_48965661",
                "ENST00000402010_63662733_63668142",
                "ENST00000216268_50247490_50279084",
                "ENST00000435706_25639475_25657066",
                "ENST00000447092_50355221_50357410",
                "ENST00000427836_208811178_208890284",
                "ENST00000372348_133589333_133755966",
                "ENST00000360284_52416616_52448857",
                "ENST00000403092_40366600_40838193",
                "ENST00000283943_230628554_230652364",
                "ENST00000373510_33005028_33059033",
                "ENST00000392485_74730197_74732960",
                "ENST00000216267_50166931_50216724",
                "ENST00000317025_38146175_38239790",
                "ENST00000263991_62900986_63221111",
                "ENST00000264122_105374305_105422873",
                "ENST00000284551_228581374_228593984",
                "ENST00000374656_33176272_33178954",
                "ENST00000370397_101953166_101977865",
                "ENST00000381638_4009075_4017636",
                "ENST00000306960_62677913_62689279",
                "ENST00000550722_112642374_112645755",
                "ENST00000425432_71894382_71929239",
                "ENST00000287295_129281817_129299861",
                "ENST00000263640_158630608_158630679",
                "ENST00000252674_6212966_6230607",
                "ENST00000357234_128577666_128586030",
                "ENST00000262483_6354584_6364798",
                "ENST00000359933_96747595_96761320",
                "ENST00000545128_78505560_78794585",
                "ENST00000366937_216676588_216850544",
                "ENST00000383202_136196251_136240123",
                "ENST00000370793_78195550_78225537",
                "ENST00000281419_9346894_9496207",
                "ENST00000568956_56126626_56126876",
                "ENST00000374888_57619317_57623906",
                "ENST00000371817_137716574_137736686",
                "ENST00000266744_103351464_103352366",
                "ENST00000360304_123328896_123401085",
                "ENST00000290894_45470370_45493373",
                "ENST00000340438_51486481_51487348",
                "ENST00000263388_15272343_15311792",
                "ENST00000260526_94645492_94654441",
                "ENST00000254854_7906931_7923657",
                "ENST00000248668_39804680_39804700",
                "ENST00000550785_113707592_113736390",
                "ENST00000259335_114122972_114246405",
                "ENST00000323686_51430165_51435330",
                "ENST00000314499_166535347_166545917",
                "ENST00000542965_30662038_30666322",
                "ENST00000434752_138730689_138730885",
                "ENST00000426804_20458506_20461786",
                "ENST00000413988_45715879_45719399",
                "ENST00000507955_176940771_176944470",
                "ENST00000395699_44915896_44924357",
                "ENST00000601417_53552985_53553711",
                "ENST00000392589_110422829_110501207",
                "ENST00000542845_109292400_109292444",
                "ENST00000420765_43825660_43851107",
                "ENST00000304718_34900737_34937842",
                "ENST00000354258_32821386_32821755",
                "ENST00000562955_43881638_43888194",
                "ENST00000371471_52183604_52193706",
                "ENST00000301764_61083827_61110068",
                "ENST00000342628_2977870_2986206",
                "ENST00000296861_47221208_47254136",
                "ENST00000370185_99730095_99775146",
                "ENST00000246117_3201615_3209573",
                "ENST00000265560_49315264_49362352",
                "ENST00000340281_90460671_90487899",
                "ENST00000506030_65222303_65344547",
                "ENST00000202556_104206565_104206805",
                "ENST00000322623_127816194_127816337",
                "ENST00000251582_178634621_178772431",
                "ENST00000420765_43851108_43882430",
                "ENST00000558170_145141648_145154127",
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
            output_ht_path = "gs://regional_missense_constraint/temp/simul_breaks/initial/round2/raw_results/simul_break_group1under.ht"
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
                count=1,
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
