import argparse
import logging

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists

from rmc.resources.basics import (
    not_one_break,
    temp_path,
)


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("search_for_two_breaks")
logger.setLevel(logging.INFO)


def search_for_two_breaks(
    ht: hl.Table,
    pos_ht: hl.Table,
    transcript: str,
    success_ht_path: str,
    scan_checkpoint_path: str,
    chisq_checkpoint_path: str,
    min_window_size: int,
    chisq_threshold: float = 13.8,
) -> None:
    """
    Search for windows of constraint in transcripts with simultaneous breaks.

    This function searches for breaks for all possible window sizes but only keeps break sizes >= `min_window_size`.
    `min_window_size` is the number of base pairs needed, on average, to see 10 missense variants (by default).

    For gnomAD v2.1, `min_window_size` is 100bp.

    :param str ht_path:Input Table filtered to contain only transcripts without one significant break.
    :param str pos_ht_path: Table containing all positions per transcript.
    :param str transcript: Transcript to search.
    :param str success_ht_path: Path to output Table.
    :param str scan_checkpoint_path: Path to checkpoint temporary HT (after completing scans).
    :param str chisq_checkpoint_path: Path to checkpoint temporary HT (after calculating null/alt distributions and associated chi square values).
    :param int min_window_size: Smallest possible window size to search for two simultaneous breaks.
    :param float chisq_threshold:  Chi-square significance threshold. Default is 13.8.
        Default is from ExAC RMC code and corresponds to a p-value of 0.999 with 2 degrees of freedom.
        (https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm)
    :return: None
    """
    ht = ht.filter(ht.transcript == transcript)

    # Reformat cumulative obs annotation from a dict ({transcript: int}) to just an int
    # Also reformat cumulative mu annotation from dict ({transcript: float}) to just float
    # TODO: Fix these annotations in the code above when they're written
    ht = ht.transmute(
        cumulative_obs=ht.cumulative_obs[ht.transcript],
        cumulative_mu=ht._mu_scan[ht.transcript],
    )

    # Create arrays of cumulative observed and expected missense values
    ht = ht.annotate(
        # NOTE: These scans do NOT need to be inclusive per row (do not need adjusting)
        # This is because these cumulative values are going to be used to calculate the
        # observed and expected missense counts for the section pre-window
        prev_obs=hl.scan.group_by(ht.transcript, hl.scan.collect(ht.cumulative_obs)),
        prev_mu=hl.scan.group_by(ht.transcript, hl.scan.collect(ht.cumulative_mu)),
    )
    # Select annotations needed for simultaneous breaks calculations
    ht = ht.select(
        "overall_oe",
        "prev_obs",
        "prev_mu",
        "total_obs",
        "total_exp",
        "total_mu",
        "reverse",
        "reverse_obs_exp",
    )
    # Checkpoint HT here because the scans are computationally intense
    ht = ht.checkpoint(scan_checkpoint_path, overwrite=True)

    # Translate mu to expected
    # ht = ht.annotate(
    #    prev_exp=hl.if_else(
    #        hl.len(ht.prev_mu) == 0,
    #        ht.prev_mu,
    #        hl.map(lambda x: (x / ht.total_mu) * ht.total_exp, ht.prev_mu),
    #    )
    # ).drop("prev_mu")
    ht = ht.annotate(
        prev_exp=hl.if_else(
            hl.is_missing(ht.prev_mu.get(ht.transcript)),
            hl.empty_array(hl.tfloat64),
            hl.map(
                lambda x: (x / ht.total_mu) * ht.total_exp, ht.prev_mu[ht.transcript]
            ),
        )
    ).drop("prev_mu")
    ht = ht.transmute(
        prev_obs=hl.if_else(
            hl.is_missing(ht.prev_obs.get(ht.transcript)),
            hl.empty_array(hl.tint64),
            ht.prev_obs[ht.transcript],
        )
    )

    # Run a quick validity check that each row has the same number of values in prev_obs and prev_exp
    mismatch_len_count = ht.aggregate(
        hl.agg.count_where(hl.len(ht.prev_obs) != hl.len(ht.prev_exp))
    )
    if mismatch_len_count != 0:
        ht = ht.filter(hl.len(ht.prev_obs) != hl.len(ht.prev_exp))
        ht.show()
        raise DataException(
            f"{mismatch_len_count} lines have different scan array lengths for prev obs and prev exp!"
        )

    # Calculate null and alt distributions
    ht = ht.annotate(
        nulls=hl.zip(ht.prev_obs, ht.prev_exp).starmap(
            lambda obs, exp: (
                hl.dpois(obs, exp * ht.overall_oe)
                * hl.dpois(
                    (ht.total_obs - ht.reverse.obs - obs),
                    (ht.total_exp - ht.reverse.exp - exp) * ht.overall_oe,
                )
                * hl.dpois(ht.reverse.obs, ht.reverse.exp * ht.overall_oe)
            )
        ),
        alts=hl.zip(ht.prev_obs, ht.prev_exp).starmap(
            lambda obs, exp: (
                hl.dpois(obs, exp * hl.min(obs / exp, 1))
                * hl.dpois(
                    (ht.total_obs - ht.reverse.obs - obs),
                    (ht.total_exp - ht.reverse.exp - exp)
                    * hl.min(
                        (ht.total_obs - ht.reverse.obs - obs)
                        / (ht.total_exp - ht.reverse.exp - exp),
                        1,
                    ),
                )
                * hl.dpois(ht.reverse.obs, ht.reverse.exp * ht.reverse_obs_exp)
            )
        ),
    )

    # Calculate chi square values
    ht = ht.annotate(
        chi_squares=hl.zip(ht.nulls, ht.alts).starmap(
            lambda null, alt: 2 * (hl.log10(alt) - hl.log10(null))
        )
    )

    # Get maximum chi square value for each position
    ht = ht.annotate(
        max_chisq_for_pos=hl.or_missing(
            hl.len(ht.chi_squares) > 0, hl.max(ht.chi_squares)
        )
    )

    # Select chi square annotations only (only required annotations) and checkpoint
    # Checkpoint here because HT branches and becomes both group_ht and max_chisq ht below
    ht = ht.select("chi_squares", "max_chisq_for_pos")
    ht = ht.checkpoint(chisq_checkpoint_path, overwrite=True)

    # Get maximum chi square value for each transcript
    group_ht = ht.group_by("transcript").aggregate(
        max_chisq_per_transcript=hl.agg.max(ht.max_chisq_for_pos)
    )
    ht = ht.annotate(
        max_chisq_per_transcript=group_ht[ht.transcript].max_chisq_per_transcript
    )

    # Find window start position associated with maximum chi square values
    max_chisq_ht = ht.filter(ht.max_chisq_for_pos == ht.max_chisq_per_transcript)
    max_chisq_ht = max_chisq_ht.annotate(
        chisq_index=max_chisq_ht.chi_squares.index(
            max_chisq_ht.max_chisq_per_transcript
        )
    )

    # Get all positions present in transcript to map chi square index back to window start position
    max_chisq_ht = max_chisq_ht.annotate(
        pos_per_transcript=pos_ht[max_chisq_ht.transcript].positions
    )
    max_chisq_ht = max_chisq_ht.annotate(
        window_start=max_chisq_ht.pos_per_transcript[max_chisq_ht.chisq_index]
    )

    # Make sure window start position is at least the same as minimum window size
    max_chisq_ht = max_chisq_ht.filter(
        (max_chisq_ht.max_chisq_for_pos >= chisq_threshold)
        & (max_chisq_ht.locus.position - max_chisq_ht.window_start >= min_window_size)
    )

    # Write output Table only if there is a significant break
    if max_chisq_ht.count() != 0:
        max_chisq_ht = max_chisq_ht.select("window_start", "max_chisq_per_transcript")
        max_chisq_ht.write(success_ht_path, overwrite=True)


def main(args):
    """Search for two simultaneous breaks in transcripts without evidence of a single significant break."""
    hl.init(log="/search_for_two_breaks.log")

    # Set hail flag to avoid method too large and out of memory errors
    hl._set_flags(no_whole_stage_codegen="1")

    # Set hail flag to speed up scans
    hl._set_flags(distributed_scan_comb_op="1")

    logger.info(
        "Searching for two simultaneous breaks in transcripts that didn't have \
        a single significant break..."
    )
    transcript_tsv_path = args.transcript_tsv

    if not file_exists(transcript_tsv_path):
        raise DataException(f"{transcript_tsv_path} doesn't exist!")

    transcripts = []
    with hl.hadoop_open(transcript_tsv_path) as i:
        for line in i:
            transcripts.append(line.strip())

    logger.info("Checking for positions per transcript HT...")
    # Table stored in temp bucket that is used to calculate constraint but can be deleted afterwards
    # This Table is grouped by transcript and has all positions in each transcript collected into a list
    pos_ht_path = f"{temp_path}/pos_per_transcript.ht"
    if not file_exists(pos_ht_path) or args.overwrite_pos_ht:
        ht = not_one_break.ht()
        pos_ht = ht.group_by("transcript").aggregate(
            positions=hl.sorted(hl.agg.collect(ht.locus.position)),
        )
        pos_ht.write(f"{temp_path}/pos_per_transcript.ht", overwrite=True)
    ht = not_one_break.ht()
    pos_ht = hl.read_table(pos_ht_path)

    transcript_success_map = {}
    transcript_ht_map = {}
    for transcript in transcripts:
        logger.info("Working on %s...", transcript)

        # Search for simultaneous breaks
        search_for_two_breaks(
            ht,
            pos_ht,
            transcript,
            f"{temp_path}/simul_breaks_{transcript}.ht",
            f"{temp_path}/simul_breaks_scan_collect_{transcript}.ht",
            f"{temp_path}/simul_breaks_chisq_{transcript}.ht",
            args.min_window_size,
            args.chisq_threshold,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script searches for regional missense constraint in gnomAD"
    )
    parser.add_argument(
        "--chisq-threshold",
        help="Chi-square significance threshold. Value should be 10.8 (single break) and 13.8 (two breaks) (values from ExAC RMC code).",
        type=float,
        default=10.8,
    )
    parser.add_argument(
        "--transcript-tsv",
        help="Path to store transcripts to search for two simultaneous breaks. Path should be to a file in Google cloud storage.",
        default=f"{temp_path}/no_break_transcripts.tsv",
    )
    parser.add_argument(
        "--overwrite-pos-ht",
        help="Overwrite the positions per transcript HT (HT keyed by transcript with a list of positiosn per transcript), even if it already exists.",
        action="store_true",
    )
    parser.add_argument(
        "--min-window-size",
        help="Smallest possible window size for simultaneous breaks. Determined by running --get-min-window-size.",
        type=int,
    )
    args = parser.parse_args()
    main(args)
