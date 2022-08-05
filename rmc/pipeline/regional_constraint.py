import argparse
import logging

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications

from rmc.resources.basics import (
    LOGGING_PATH,
    multiple_breaks,
    one_break,
    simul_break,
    temp_path,
)
from rmc.resources.grch37.gnomad import constraint_ht
from rmc.resources.grch37.reference_data import filtered_context
from rmc.resources.resource_utils import CURRENT_VERSION
from rmc.slack_creds import slack_token
from rmc.utils.constraint import process_additional_breaks


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_missense_constraint")
logger.setLevel(logging.INFO)


def get_constraint_transcripts(outlier: bool = True) -> hl.expr.SetExpression:
    """
    Read in LoF constraint HT results to get set of transcripts.

    Return either set of transcripts to keep (transcripts that passed transcript QC)
    or outlier transcripts.

    Transcripts are removed for the reasons detailed here:
    https://gnomad.broadinstitute.org/faq#why-are-constraint-metrics-missing-for-this-gene-or-annotated-with-a-note

    :param bool outlier: Whether to filter LoF constraint HT to outlier transcripts (if True),
        or QC-pass transcripts (if False). Default is True.
    :return: Set of outlier transcripts or transcript QC pass transcripts.
    :rtype: hl.expr.SetExpression
    """
    logger.warning(
        "Assumes LoF constraint has been separately calculated and that constraint HT exists..."
    )
    if not file_exists(constraint_ht.path):
        raise DataException("Constraint HT not found!")

    constraint_transcript_ht = constraint_ht.ht().key_by("transcript")
    constraint_transcript_ht = constraint_transcript_ht.filter(
        constraint_transcript_ht.canonical
    ).select("constraint_flag")
    if outlier:
        constraint_transcript_ht = constraint_transcript_ht.filter(
            hl.len(constraint_transcript_ht.constraint_flag) > 0
        )
    else:
        constraint_transcript_ht = constraint_transcript_ht.filter(
            hl.len(constraint_transcript_ht.constraint_flag) == 0
        )
    return hl.literal(
        constraint_transcript_ht.aggregate(
            hl.agg.collect_as_set(constraint_transcript_ht.transcript)
        )
    )


def main(args):
    """Call functions from `constraint.py` to calculate regional missense constraint."""
    exac = args.exac

    try:
        if args.search_for_additional_breaks:
            hl.init(log="/RMC_additional_breaks.log")

            # Set hail flag to avoid method too large and out of memory errors
            hl._set_flags(no_whole_stage_codegen="1")

            logger.info(
                "Searching for additional breaks in transcripts with at least one significant break..."
            )
            context_ht = one_break.ht()

            # Remove outlier transcripts
            outlier_transcripts = get_constraint_transcripts(outlier=True)
            context_ht = context_ht.filter(
                ~outlier_transcripts.contains(context_ht.transcript)
            )

            # Add break_list annotation to context HT
            context_ht = context_ht.annotate(break_list=[context_ht.is_break])
            break_ht = context_ht

            # Start break number counter at 2
            break_num = 2

            while True:
                # Search for additional breaks
                # This technically should search for two additional breaks at a time:
                # this calls `process_sections`, which checks each section of the transcript for a break
                # sections are transcript section pre and post first breakpoint
                break_ht = process_additional_breaks(
                    break_ht, break_num, args.chisq_threshold
                )
                break_ht = break_ht.checkpoint(
                    f"{temp_path}/break_{break_num}.ht", overwrite=True
                )

                # Filter context HT to lines with break and check for transcripts with at least one additional break
                is_break_ht = break_ht.filter(break_ht.is_break)
                group_ht = is_break_ht.group_by("transcript").aggregate(
                    n=hl.agg.count()
                )
                group_ht = group_ht.filter(group_ht.n >= 1)

                # Exit loop if no additional breaks are found for any transcripts
                if group_ht.count() == 0:
                    break

                # Otherwise, pull transcripts and annotate context ht
                break_ht = break_ht.key_by("locus", "transcript")
                transcripts = group_ht.aggregate(
                    hl.agg.collect_as_set(group_ht.transcript), _localize=False
                )
                globals_annot_expr = {f"break_{break_num}_transcripts": transcripts}
                context_ht = context_ht.annotate_globals(**globals_annot_expr)
                annot_expr = {
                    f"break_{break_num}_chisq": break_ht[context_ht.key].chisq,
                    f"break_{break_num}_max_chisq": break_ht[context_ht.key].max_chisq,
                    f"break_{break_num}_null": break_ht[context_ht.key].total_null,
                    f"break_{break_num}_alt": break_ht[context_ht.key].total_alt,
                    "is_break": break_ht[context_ht.key].is_break,
                }
                context_ht = context_ht.annotate(**annot_expr)
                context_ht = context_ht.annotate(
                    break_list=context_ht.break_list.append(context_ht.is_break)
                )

                break_ht = break_ht.filter(transcripts.contains(break_ht.transcript))
                break_num += 1

            # context_ht.write(multiple_breaks.path, overwrite=args.overwrite)
            context_ht.write(f"{temp_path}/multiple_breaks.ht", overwrite=True)

        if args.finalize:
            hl.init(log="/RMC_finalize.log")

            logger.info(
                "Getting start and end positions and total size for each transcript..."
            )
            if not file_exists(f"{temp_path}/transcript.ht"):
                raise DataException(
                    "Transcript HT doesn't exist. Please double check and recreate!"
                )

            if args.remove_outlier_transcripts:
                outlier_transcripts = get_constraint_transcripts(outlier=True)

            logger.info("Reading in context HT...")
            # Drop extra annotations from context HT
            context_ht = filtered_context.ht().drop(
                "_obs_scan", "_mu_scan", "forward_oe", "values"
            )
            context_ht = context_ht.filter(
                ~outlier_transcripts.contains(context_ht.transcript)
            )

            logger.info("Reading in results HTs...")
            multiple_breaks_ht = multiple_breaks.ht()
            multiple_breaks_ht = multiple_breaks_ht.filter(
                ~outlier_transcripts.contains(multiple_breaks_ht.transcript)
            )
            # TODO: Update simul breaks reformatting
            simul_breaks_ht = simul_break.ht()
            simul_breaks_ht = simul_breaks_ht.filter(
                ~outlier_transcripts.contains(simul_breaks_ht.transcript)
            )

            logger.info("Getting simultaneous breaks transcripts...")
            simul_break_transcripts = simul_breaks_ht.aggregate(
                hl.agg.collect_as_set(simul_breaks_ht.transcript), _localize=False
            )
            rmc_transcripts = simul_break_transcripts

            logger.info("Getting transcript information for breaks HT...")
            global_fields = multiple_breaks_ht.globals
            n_breaks = [
                int(x.split("_")[1]) for x in global_fields if "break" in x
            ].sort()
            rmc_transcripts = []
            for break_num in n_breaks:
                ht = hl.read_table(f"{temp_path}/break_{break_num}.ht")
                ht = ht.filter(ht.is_break)
                rmc_transcripts.append(
                    ht.aggregate(hl.agg.collect_as_set(ht.transcript), _localize=False)
                )

            logger.info("Removing overlapping transcript information...")
            # Extracting the unique set of transcripts for each break number
            # e.g., keeping only transcripts with a single break
            # and not transcripts also that had additional breaks in "break_1_transcripts"
            filtered_transcripts = {}

            logger.info(
                "Cycling through each set of transcripts with break information..."
            )
            for index, transcripts in enumerate(rmc_transcripts):
                # Cycle through the sets of transcripts for every number of breaks larger than current number of breaks
                for t in rmc_transcripts[index + 1 :]:
                    transcripts = transcripts.difference(t)

                # Separate transcripts with a single break from transcripts with multiple breaks
                # This is because transcripts with a single break are annotated into a separate struct
                # in the release HT globals
                if index == 0:
                    one_break_transcripts = transcripts
                else:
                    filtered_transcripts[f"break_{index + 1}_transcripts"] = transcripts

            # Getting total number of transcripts with evidence of rmc
            total_rmc_transcripts = simul_break_transcripts
            for break_num in filtered_transcripts:
                total_rmc_transcripts = total_rmc_transcripts.union(
                    filtered_transcripts[break_num]
                )
            logger.info(
                "Number of transcripts with evidence of RMC: %s",
                hl.eval(hl.len(total_rmc_transcripts)),
            )

            logger.info("Creating breaks HT...")
            breaks_ht = context_ht.filter(
                total_rmc_transcripts.contains(context_ht.transcript)
            )

            logger.info("Adding simultaneous breaks information to breaks HT...")
            breaks_ht = breaks_ht.annotate(
                simul_break_info=hl.struct(
                    is_break=simul_breaks_ht[breaks_ht.key].is_break,
                    window_end=simul_breaks_ht[breaks_ht.key].window_end,
                    post_window_pos=simul_breaks_ht[breaks_ht.key].post_window_pos,
                )
            )
            logger.info("Adding transcript information to breaks HT globals...")
            breaks_ht = breaks_ht.annotate_globals(
                single_break=hl.struct(transcripts=one_break_transcripts),
                multiple_breaks=hl.struct(**filtered_transcripts),
                simul_breaks=hl.struct(transcripts=simul_break_transcripts),
            )

            logger.info("Creating no breaks HT...")
            no_breaks_ht = context_ht.filter(
                ~total_rmc_transcripts.contains(context_ht.transcript)
            )

            if CURRENT_VERSION == "2.1.1":
                logger.info("Reading in XG HT (one-off fix in v2.1.1)...")
                xg_ht = hl.read_table(f"{temp_path}/XG.ht").select(
                    "total_mu",
                    "total_exp",
                    "total_obs",
                    "cumulative_exp",
                    "cumulative_obs",
                    "overall_oe",
                )
                no_breaks_ht = no_breaks_ht.annotate(
                    xg_total_mu=xg_ht[no_breaks_ht.key].total_mu,
                    xg_total_exp=xg_ht[no_breaks_ht.key].total_exp,
                    xg_total_obs=xg_ht[no_breaks_ht.key].total_obs,
                    xg_cum_exp=xg_ht[no_breaks_ht.key].cumulative_exp,
                    xg_cum_obs=xg_ht[no_breaks_ht.key].cumulative_obs,
                    xg_oe=xg_ht[no_breaks_ht.key].overall_oe,
                )
                no_breaks_ht = no_breaks_ht.transmute(
                    total_mu=hl.if_else(
                        hl.is_defined(no_breaks_ht.xg_total_mu),
                        no_breaks_ht.xg_total_mu,
                        no_breaks_ht.total_mu,
                    ),
                    total_exp=hl.if_else(
                        hl.is_defined(no_breaks_ht.xg_total_exp),
                        no_breaks_ht.xg_total_exp,
                        no_breaks_ht.total_exp,
                    ),
                    total_obs=hl.if_else(
                        hl.is_defined(no_breaks_ht.xg_total_obs),
                        no_breaks_ht.xg_total_obs,
                        no_breaks_ht.total_obs,
                    ),
                    cumulative_exp=hl.if_else(
                        hl.is_defined(no_breaks_ht.xg_cum_exp),
                        no_breaks_ht.xg_cum_exp,
                        no_breaks_ht.cumulative_exp,
                    ),
                    cumulative_obs=hl.if_else(
                        hl.is_defined(no_breaks_ht.xg_cum_obs),
                        no_breaks_ht.xg_cum_obs,
                        no_breaks_ht.cumulative_obs,
                    ),
                    overall_oe=hl.if_else(
                        hl.is_defined(no_breaks_ht.xg_oe),
                        no_breaks_ht.xg_oe,
                        no_breaks_ht.overall_oe,
                    ),
                )

            # TODO: Add section chisq calculation here
            if (breaks_ht.count() + no_breaks_ht.count()) != context_ht.count():
                raise DataException(
                    "Row counts for breaks HT (one break, multiple breaks, simul breaks) and no breaks HT doesn't match context HT row count!"
                )

            logger.info("Checkpointing HTs...")
            breaks_ht = breaks_ht.checkpoint(f"{temp_path}/breaks.ht", overwrite=True)
            no_breaks_ht = no_breaks_ht.checkpoint(
                f"{temp_path}/no_breaks.ht", overwrite=True
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(LOGGING_PATH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script searches for regional missense constraint in gnomAD"
    )
    parser.add_argument(
        "--trimers", help="Use trimers instead of heptamers", action="store_true"
    )
    parser.add_argument(
        "--chisq-threshold",
        help="Chi-square significance threshold. Value should be 10.8 (single break) and 13.8 (two breaks) (values from ExAC RMC code).",
        type=float,
        default=10.8,
    )
    parser.add_argument(
        "--search-for-additional-breaks",
        help="Search for additional break in transcripts with one significant break",
        action="store_true",
    )
    parser.add_argument(
        "--finalize",
        help="Combine and reformat (finalize) RMC output",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite-transcript-ht",
        help="Overwrite the transcript HT (HT with start/end positions and transcript sizes), even if it already exists.",
        action="store_true",
    )
    parser.add_argument(
        "--remove-outlier-transcripts",
        help="Remove outlier transcripts (transcripts with too many/few LoF, synonymous, or missense variants)",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite", help="Overwrite existing data", action="store_true"
    )
    parser.add_argument(
        "--slack-channel",
        help="Send message to Slack channel/user",
        default="@kc (she/her)",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
