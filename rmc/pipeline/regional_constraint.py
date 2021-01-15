import argparse
import logging

import hail as hl

from gnomad.utils.slack import slack_notifications
from rmc.resources.basics import (
    constraint_prep,
    LOGGING_PATH,
    multiple_breaks,
    not_one_break,
    one_break,
    simul_break,
    temp_path,
)
from rmc.resources.grch37.exac import filtered_exac
from rmc.resources.grch37.gnomad import (
    constraint_ht,
    filtered_exomes,
    processed_exomes,
    prop_obs_coverage,
)
from rmc.resources.grch37.reference_data import processed_context
from rmc.resources.resource_utils import MISSENSE
from rmc.slack_creds import slack_token
from rmc.utils.constraint import (
    calculate_exp_per_transcript,
    calculate_observed,
    get_fwd_exprs,
    GROUPINGS,
    process_additional_breaks,
    process_transcripts,
    search_for_two_breaks,
)
from rmc.utils.generic import (
    filter_to_region_type,
    generate_models,
    get_coverage_correction_expr,
    keep_criteria,
    process_context_ht,
    process_vep,
)


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_missense_constraint")
logger.setLevel(logging.INFO)


def main(args):

    hl.init(log="/RMC.log")
    exac = args.exac

    try:
        if args.pre_process_data:
            # TODO: Add code to create annotations necessary for constraint_flag_expr and filter transcripts prior to running constraint
            logger.warning("Code currently only processes b37 data!")
            logger.info(
                "Filtering gnomAD exomes HT to missense variants in canonical transcripts only..."
            )
            exome_ht = processed_exomes.ht()
            exome_ht = process_vep(exome_ht, filter_csq=True, csq=MISSENSE)

            # Move nested annotations into top level annotations
            exome_ht = exome_ht.select(
                ac=exome_ht.freq[0].AC,
                af=exome_ht.freq[0].AF,
                pass_filters=exome_ht.pass_filters,
                exome_coverage=exome_ht.coverage.exomes.median,
                transcript_consequences=exome_ht.transcript_consequences,
            )

            logger.info("Preprocessing reference fasta (context) HT...")
            context_ht = process_context_ht("GRCh37", args.trimers)

            logger.info(
                "Filtering context HT to all sites not found in gnomAD exomes + all rare, covered sites in gnomAD"
            )
            exome_join = exome_ht[context_ht.key]
            context_ht = context_ht.filter(
                hl.is_missing(exome_join) | keep_criteria(exome_join)
            )
            # NOTE: need to repartition here to desired number of partitions!
            # NOTE: should use ~30k-40k partitions
            context_ht = context_ht.repartition(args.n_partitions)
            context_ht.write(processed_context.path, overwrite=args.overwrite)

            exome_ht = exome_ht.filter(keep_criteria(exome_ht))
            exome_ht.write(filtered_exomes.path, overwrite=args.overwrite)

            logger.info("Done preprocessing files")

        if args.prep_for_constraint:
            logger.info("Reading in exome HT...")
            if exac:
                exome_ht = filtered_exac.ht()

            else:
                exome_ht = filtered_exomes.ht()

            logger.info("Reading in context HT...")
            context_ht = processed_context.ht()

            logger.info("Building plateau and coverage models...")
            coverage_ht = prop_obs_coverage.ht()
            coverage_x_ht = hl.read_table(
                prop_obs_coverage.path.replace(".ht", "_x.ht")
            )
            coverage_y_ht = hl.read_table(
                prop_obs_coverage.path.replace(".ht", "_y.ht")
            )

            (
                coverage_model,
                plateau_models,
                plateau_x_models,
                plateau_y_models,
            ) = generate_models(
                coverage_ht, coverage_x_ht, coverage_y_ht, trimers=args.trimers
            )

            context_ht = context_ht.annotate_globals(
                plateau_models=plateau_models,
                plateau_x_models=plateau_x_models,
                plateau_y_models=plateau_y_models,
                coverage_model=coverage_model,
            )
            context_ht = context_ht.annotate(
                coverage_correction=get_coverage_correction_expr(
                    context_ht.exome_coverage, coverage_model, args.high_cov_cutoff,
                )
            )

            if not args.skip_calc_oe:
                logger.info(
                    "Adding coverage correction to mutation rate probabilities..."
                )
                context_ht = context_ht.annotate(
                    raw_mu_snp=context_ht.mu_snp,
                    mu_snp=context_ht.mu_snp
                    * get_coverage_correction_expr(
                        context_ht.exome_coverage, context_ht.coverage_model
                    ),
                )

                logger.info(
                    "Creating autosomes-only, chrX non-PAR-only, and chrY non-PAR-only HT versions..."
                )
                context_x_ht = filter_to_region_type(context_ht, "chrX")
                context_y_ht = filter_to_region_type(context_ht, "chrY")
                context_auto_ht = filter_to_region_type(context_ht, "autosomes")

                logger.info("Calculating expected values per transcript...")
                exp_ht = calculate_exp_per_transcript(
                    context_auto_ht, locus_type="autosomes", groupings=GROUPINGS,
                )
                exp_x_ht = calculate_exp_per_transcript(
                    context_x_ht, locus_type="X", groupings=GROUPINGS,
                )
                exp_y_ht = calculate_exp_per_transcript(
                    context_y_ht, locus_type="Y", groupings=GROUPINGS,
                )
                exp_ht = exp_ht.union(exp_x_ht).union(exp_y_ht)

                logger.info(
                    "Aggregating total observed variant counts per transcript..."
                )
                obs_ht = calculate_observed(exome_ht)

            else:
                logger.warning(
                    "Using observed and expected values calculated on gnomAD v2.1.1 exomes..."
                )
                exp_ht = (
                    constraint_ht.ht()
                    .key_by("transcript")
                    .select("obs_mis", "exp_mis", "oe_mis")
                )

                # Filter to canonical transcripts only and rename fields
                exp_ht = exp_ht.filter(exp_ht.canonical)
                exp_ht = exp_ht.transmute(
                    observed=exp_ht.obs_mis,
                    expected=exp_ht.exp_mis,
                    mu_agg=exp_ht.mu_mis,
                )
                obs_ht = exp_ht

            logger.info(
                "Annotating context HT with number of observed and expected variants per site..."
            )
            # Add observed variants to context HT
            context_ht = context_ht.annotate(_obs=exome_ht.index(context_ht.key))
            context_ht = context_ht.transmute(
                observed=hl.int(hl.is_defined(context_ht._obs))
            )

            logger.info(
                "Collecting by key to run constraint per base and not per base-allele..."
            )
            # Context HT is keyed by locus and allele, which means there is one row for every possible missense variant
            # This means that any locus could be present up to three times (once for each possible missense)
            # Collect by key here to ensure all loci are unique
            context_ht = context_ht.key_by("locus", "transcript").collect_by_key()
            context_ht = context_ht.annotate(
                # Collect the mutation rate probabilities at each locus
                mu_snp=hl.sum(context_ht.values.mu_snp),
                # Collect the observed counts for each locus
                # (this includes counts for each possible missense at the locus)
                observed=hl.sum(context_ht.values.observed),
                # Take just the first coverage value, since the locus should have the same coverage across the possible variants
                coverage=context_ht.values.exome_coverage[0],
            )

            logger.info(
                "Annotating total observed and expected values and overall observed/expected value "
                "(capped at 1) per transcript..."
            )
            context_ht = context_ht.annotate(
                total_exp=exp_ht[context_ht.transcript].expected,
                total_mu=exp_ht[context_ht.transcript].mu_agg,
                total_obs=obs_ht[context_ht.transcript].observed,
            )
            context_ht = context_ht.annotate(
                overall_oe=hl.min(context_ht.total_obs / context_ht.total_exp, 1)
            )

            context_ht = get_fwd_exprs(
                ht=context_ht,
                transcript_str="transcript",
                obs_str="observed",
                mu_str="mu_snp",
                total_mu_str="total_mu",
                total_exp_str="total_exp",
            )

            context_ht = context_ht.write(
                constraint_prep.path, overwrite=args.overwrite
            )

        if args.search_for_first_break:
            logger.info("Searching for transcripts with a significant break...")
            context_ht = constraint_prep.ht()
            context_ht = process_transcripts(context_ht, args.chisq_threshold)
            context_ht = context_ht.checkpoint(
                f"{temp_path}/first_break.ht", overwrite=True
            )

            logger.info(
                "Filtering HT to transcripts with one significant break and writing..."
            )
            is_break_ht = context_ht.filter(context_ht.is_break)
            transcripts = is_break_ht.aggregate(
                hl.agg.collect_as_set(is_break_ht.transcript), _localize=False
            )
            one_break_ht = context_ht.filter(
                transcripts.contains(context_ht.transcript)
            )
            one_break_ht = one_break_ht.annotate_globals(
                break_1_transcripts=transcripts
            )
            one_break_ht.repartition(args.n_partitions).write(
                one_break.path, overwrite=args.overwrite
            )

            logger.info(
                "Filtering HT to transcripts without a significant break and writing..."
            )
            not_one_break_ht = context_ht.anti_join(one_break_ht)
            not_one_break_ht.repartition(args.n_partitions).write(
                not_one_break.path, overwrite=args.overwrite
            )

        if args.search_for_additional_breaks:
            logger.info(
                "Searching for additional breaks in transcripts with at least one significant break..."
            )
            context_ht = one_break.ht()

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
                break_ht = break_ht.annotate(
                    break_list=break_ht.break_list.append(break_ht.is_break)
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
                    f"break_{break_num}_null": break_ht[context_ht.key].total_null,
                    f"break_{break_num}_alt": break_ht[context_ht.key].total_alt,
                }
                context_ht = context_ht.annotate(**annot_expr)

                break_ht = break_ht.filter(transcripts.contains(break_ht.transcript))
                break_num += 1

            context_ht.write(multiple_breaks.path, overwrite=args.overwrite)

        if args.search_for_simul_breaks:
            logger.info(
                "Searching for two simultaneous breaks in transcripts that didn't have \
                a single significant break..."
            )
            context_ht = not_one_break.ht()
            exome_ht = filtered_exomes.ht()

            logger.info("Getting start and end positions for each transcript...")
            transcript_ht = context_ht.group_by(context_ht.transcript).aggregate(
                end_pos=hl.agg.max(context_ht.locus.position),
                start_pos=hl.agg.min(context_ht.locus.position),
            )
            context_ht = context_ht.annotate(
                start_pos=transcript_ht[context_ht.transcript].start_pos,
                end_pos=transcript_ht[context_ht.transcript].end_pos,
            )

            logger.info("Searching for transcripts with simultaneous breaks...")
            all_transcripts = context_ht.aggregate(
                hl.agg.collect_as_set(context_ht.transcripts), _localize=False
            )
            ht = search_for_two_breaks(
                ht=context_ht,
                exome_ht=exome_ht,
                chisq_threshold=args.chisq_threshold,
                num_obs_var=args.num_obs_var,
            )
            is_break_ht = ht.filter(ht.is_break)
            transcripts = is_break_ht.aggregate(
                hl.agg.collect_as_set(is_break_ht.transcript), _localize=False
            )

            logger.info("Annotating globals and writing...")
            ht = ht.annotate_globals(
                no_break_transcripts=all_transcripts.difference(transcripts),
                simul_break_transcripts=transcripts,
            )
            ht.repartition(args.n_partitions).write(simul_break.path)

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
        "--exac", help="Use ExAC Table (not gnomAD Table)", action="store_true"
    )
    parser.add_argument(
        "--n_partitions",
        help="Desired number of partitions for output data",
        type=int,
        default=40000,
    )
    parser.add_argument(
        "--high_cov_cutoff",
        help="Coverage threshold for a site to be considered high coverage",
        type=int,
        default=40,
    )
    parser.add_argument(
        "--chisq_threshold",
        help="Chi-square significance threshold. Value should be 10.8 (single break) and 13.8 (two breaks) (values from ExAC RMC code).",
        type=float,
        default=10.8,
    )
    parser.add_argument(
        "--num_obs_var",
        help="Number of observed variants. Used when determining the window size for simultaneous breaks.",
        type=int,
        default=10,
    )
    parser.add_argument(
        "--pre_process_data", help="Pre-process data", action="store_true"
    )
    parser.add_argument(
        "--prep_for_constraint",
        help="Prepare tables for constraint calculations",
        action="store_true",
    )
    parser.add_argument(
        "--skip_calc_oe",
        help="Skip observed and expected variant calculations per transcript. Relevant only to gnomAD v2.1.1!",
        action="store_true",
    )
    parser.add_argument(
        "--search_for_first_break",
        help="Initial search for one break in all transcripts",
        action="store_true",
    )
    parser.add_argument(
        "--search_for_additional_breaks",
        help="Search for additional break in transcripts with one significant break",
        action="store_true",
    )
    parser.add_argument(
        "--search_for_simul_breaks",
        help="Search for two simultaneous breaks in transcripts without a single significant break",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite", help="Overwrite existing data", action="store_true"
    )
    parser.add_argument(
        "--slack_channel", help="Send message to Slack channel/user", default="@kc"
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
