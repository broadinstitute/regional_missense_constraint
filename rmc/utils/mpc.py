import logging
from typing import Dict, List, Tuple

import hail as hl

from gnomad.resources.grch37.gnomad import public_release
from gnomad.resources.grch37.reference_data import vep_context

from rmc.resources.basics import (
    blosum,
    blosum_txt_path,
    grantham,
    grantham_txt_path,
    misbad,
    temp_path,
)
from rmc.resources.grch37.reference_data import clinvar_path_mis
from rmc.utils.generic import get_aa_map, process_vep
from rmc.utils.missense_badness import filter_codons, get_oe_annotation


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("mpc_utils")
logger.setLevel(logging.INFO)


def convert_score_list_to_ht(
    score_list: List[Dict[str, str]],
    schema: str = "array<struct{amino_acids: str, score: str}>",
    key_fields: Tuple[str] = ("ref", "alt"),
) -> hl.Table:
    """
    Convert list of amino acid changes/associated scores to Table format.

    :param List[Dict[str, str]] score_list: List of dictionaries containing amino acid changes (key) and associated scores (value).
    :param str schema: Schema of `score_list`. Default is 'array<struct{amino_acids: str, score: str}>'.
        Note that the dictionary keys must match values provided in this schema.
    :param str key_fields: Desired key fields for the new Table. Default is ("ref", "alt").
    """
    ht = hl.Table.parallelize(hl.literal(score_list, schema))
    if schema == "array<struct{amino_acids: str, score: str}>":
        ht = ht.transmute(
            ref=ht.amino_acids.split("_")[0], alt=ht.amino_acids.split("_")[1]
        )
    return ht.key_by(*key_fields)


def import_blosum():
    """
    Import BLOSUM score.

    Read in text file, convert to hail Table format, and write to resource path.

    :return: None; function writes HT to resource path.
    """
    # Get amino acid map (map 1 letter code to 3 letter code)
    aa_map = get_aa_map()

    # Create empty list to store BLOSUM scores
    # Will use this list later in the function to directly convert the scores into a Table format
    blosum_scores = []
    with hl.hadoop_open(blosum_txt_path) as b:
        for line in b:
            # Skip metadata header lines
            # e.g., # Matrix made by matblas from blosum62.iij
            if not line.startswith("#"):
                # Parse header line (starts with '.')
                if line.startswith("."):
                    header = line.strip().split("\t")
                    header_dict = {}

                    for counter, item in enumerate(header[1:]):
                        # Change asterisk to STOP
                        if item == "*":
                            item = "STOP"

                        # Store amino acid in header dict
                        header_dict[counter] = item

                else:
                    line = line.strip().split("\t")
                    # Get amino acid 1 letter code (and change asterisk to STOP)
                    aa = line[0]
                    if aa == "*":
                        aa = "STOP"

                    # Skip any amino acids that aren't in the amino acid map
                    # There are three in this file: B, X, Z
                    try:
                        aa = aa_map[aa]
                    except KeyError:
                        continue

                    for counter, item in enumerate(line[1:]):
                        alt_aa = header_dict[counter]
                        try:
                            alt_aa = aa_map[alt_aa]
                        except KeyError:
                            continue

                        # Add amino acid change and score to list
                        blosum_scores.append(
                            {"amino_acids": f"{aa}_{alt_aa}", "score": item}
                        )

    # Convert list of dictionaries to hail Table
    ht = convert_score_list_to_ht(blosum_scores)
    ht.write(blosum.path)


def import_grantham():
    """
    Import Grantham score.

    Read in text file, convert to hail Table format, and write to resource path.

    :return: None; function writes HT to resource path.
    """
    # Get amino acid map (map 1 letter code to 3 letter code)
    aa_map = get_aa_map()

    # Create empty list to store Grantham scores
    # Will use this list later in the function to directly convert the scores into a Table format
    grantham_scores = []
    with open(grantham_txt_path) as g:
        for line in g:
            # Grab header line (starts with '.')
            if line.startswith("."):
                header = line.strip().split("\t")
                header_dict = {}
                for counter, item in enumerate(header):
                    if item == ".":
                        header_dict[counter] = item
                    else:
                        header_dict[counter] = aa_map[item]
            else:
                line = line.strip().split("\t")
                aa = aa_map[line[0]]

                for counter, item in enumerate(line):
                    alt_aa = header_dict[counter]
                    if alt_aa != ".":
                        grantham_scores.append(
                            {"amino_acids": f"{aa}-{alt_aa}", "score": item}
                        )

    # Convert list of dictionaries to hail Table
    ht = convert_score_list_to_ht(grantham_scores)
    ht.write(grantham.path)


def prepare_pop_path_ht(
    gnomad_data_type: str = "exomes", af_threshold: float = 0.01
) -> None:
    """
    Prepare Table with 'population' (common gnomAD missense) and 'pathogenic' (ClinVar pathogenic/likely pathogenic missense) variants.

    .. note::
        This function reads in data from a requester-pays bucket and will fail if requester-pays
        is not enabled on the cluster.

    :param str gnomad_data_type: gnomAD data type. Used to retrieve public release Table.
        Must be one of "exomes" or "genomes" (check is done within `public_release`).
        Default is "exomes".
    :param float af_threshold: Allele frequency cutoff to filter gnomAD public dataset.
        Variants *above* this threshold will be kept.
        Default is 0.01.
    """
    logger.info("Reading in ClinVar P/LP missense variants in severe HI genes...")
    clinvar_ht = clinvar_path_mis.ht()
    clinvar_ht = clinvar_ht.annotate(pop_v_path="pop")

    logger.info("Importing gnomAD public data and filtering to common variants...")
    gnomad_ht = public_release(gnomad_data_type).ht()
    gnomad_ht = gnomad_ht.filter(gnomad_ht.freq[0].AF > af_threshold)
    gnomad_ht = gnomad_ht.annotate(pop_v_path="path")

    logger.info("Joining ClinVar and gnomAD HTs...")
    ht = clinvar_ht.select("pop_v_path").union(gnomad_ht.select("pop_v_path"))
    ht = ht.checkpoint(f"{temp_path}/joint_clinvar_gnomad.ht", overwrite=True)

    logger.info("Adding CADD, BLOSUM, Grantham, RMC annotations and checkpointing...")
    # CADD (not sure if it needs to be split)
    cadd = hl.experimental.load_dataset(
        name="cadd", version="1.6", reference_genome="GRCh37"
    )
    cadd = hl.split_multi(cadd)
    ht = ht.annotate(cadd=hl.struct(**cadd[ht.key]))
    # BLOSUM and Grantham
    blosum_ht = blosum.ht()
    grantham_ht = grantham.ht()
    ht = ht.annotate(blosum=blosum_ht[ht.key].score, grantham=grantham_ht[ht.key].score)
    # Missense observed/expected (OE) ratio
    ht = get_oe_annotation(ht)
    ht = ht.checkpoint(f"{temp_path}/joint_clinvar_gnomad_temp.ht", overwrite=True)

    logger.info("Getting PolyPhen-2 and codon annotations from VEP context HT...")
    context_ht = vep_context.ht().select_globals().select("vep", "was_split")
    context_ht = context_ht.filter(hl.is_defined(ht[context_ht].key))
    context_ht = process_vep(context_ht)
    context_ht = context_ht.annotate(
        polyphen=hl.struct(
            prediction=ht.transcript_consequences.polyphen_prediction,
            score=ht.transcript_consequences.polyphen_score,
        )
    )
    context_ht = context_ht.select(
        "polyphen", codons=context_ht.transcript_consequences.codons
    )
    context_ht = filter_codons(context_ht)
    context_ht = context_ht.checkpoint(f"{temp_path}/polyphen.ht", overwrite=True)

    logger.info("Adding PolyPhen-2 and codon annotations to joint ClinVar/gnomAD HT...")
    ht = ht.annotate(**context_ht[ht.key])

    logger.info("Getting missense badness annotation...")
    mb_ht = misbad.ht()
    ht = ht.annotate(misbad=mb_ht[ht.ref, ht.alt].misbad)

    logger.info("Filtering to rows with defined annotations and checkpointing...")
    ht = ht.filter(
        hl.is_defined(ht.cadd.phred)
        & hl.is_defined(ht.blosum)
        & hl.is_defined(ht.grantham)
        & hl.is_defined(ht.oe)
        & hl.is_defined(ht.polyphen.score)
    )
    ht = ht.checkpoint(f"{temp_path}/joint_clinvar_gnomad.ht", overwrite=True)
