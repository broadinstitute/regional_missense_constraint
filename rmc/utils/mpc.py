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
    joint_clinvar_gnomad,
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
    :return: None; function writes Table to resource path.
    """
    logger.info("Reading in ClinVar P/LP missense variants in severe HI genes...")
    clinvar_ht = clinvar_path_mis.ht()
    clinvar_ht = clinvar_ht.annotate(pop_v_path=0)

    logger.info("Importing gnomAD public data and filtering to common variants...")
    gnomad_ht = public_release(gnomad_data_type).ht()
    gnomad_ht = gnomad_ht.filter(gnomad_ht.freq[0].AF > af_threshold)
    gnomad_ht = gnomad_ht.annotate(pop_v_path=1)

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
        & hl.is_defined(ht.misbad)
    )
    ht.write(joint_clinvar_gnomad.path, overwrite=True)


def run_regressions(
    regression_type: str,
    output_fname: str,
    variables: List[str] = ["oe", "misbad", "polyphen"],
    additional_variables: List[str] = ["blosum", "grantham"],
) -> None:
    """
    Run single variable and joint regressions and pick best model.

    These regressions are used to determine the fitted score that is used to predict MPC scores.

    Fitted score from ExAC:
        fitted_score = 4.282793 + (4.359682*obs_exp) + (-3.654815*misbad) + (-3.512215*pph2)
                    + (2.585361*obs_exp*misbad) + (1.350056*obs_exp*pph2)

    Relationship between fitted score and MPC (from ExAC):
        mpc = -log10(n_less)/82932)
        n_less = number of ExAC variants with fitted_score < a given fitted_score

    :return: None; function writes model coefficients to local file.
    """
    import pandas as pd
    from patsy import dmatrices
    import statsmodels
    import statsmodels.api as sm

    assert regression_type in {
        "single",
        "additive",
        "all",
        "specific",
    }, "regression_type must be one of 'single', 'additive', 'all', 'specific'!"

    # Convert HT to pandas dataframe as logistic regression aggregations aren't currently possible in hail
    ht = joint_clinvar_gnomad.ht()
    df = ht.to_pandas()

    def _run_glm(
        formula: str,
    ) -> statsmodels.genmod.generalized_linear_model.GLMResultsWrapper:
        """
        Run logistic regression using input formula and return model results.

        MPC formula (from ExAC):
        `pop_v_path ~ obs_exp + mis_badness3 + obs_exp:mis_badness3 + polyphen2 + obs_exp:polyphen2`

        For formula reference, see: https://learn-scikit.oneoffcoder.com/patsy.html.

        :param str formula: String containing R-style formula defining model.
        :return: Logistic regression results.
        """
        # Create design matrices to fit model
        # NOTE: If run into a TypeError here, fix with df['field'].astype()
        # Example error: TypeError: Cannot interpret 'string[python]' as a data type
        # Example fix: df['misbad3'] = df['misbad3'].astype(float)
        y, X = dmatrices(formula, data=df, return_type="dataframe")
        model = sm.GLM(y, X, family=sm.families.Binomial()).fit()
        logger.info("%s summary: %s", formula, model.summary())
        logger.info("AIC: %i", model.aic)
        return model

    logger.info("Run single variable regressions...")
    # E.g., mod.misbad3 <- glm(pop_v_path ~ mis_badness3, data=cleaned_joint_exac_clinvar.scores, family=binomial)
    single_var_res = {}
    single_var_aic = []
    for var in [variables] + [additional_variables]:
        logger.info("Running single variable regression...")
        # Create design matrices
        formula = f"pop_v_path ~ {var}"
        model = _run_glm(formula)
        single_var_aic.append(model.aic)
        single_var_res[var] = model.params

    # Find lowest AIC for single variable regressions and corresponding model
    min_single_aic = min(single_var_aic)
    min_single_aic_var = variables[single_var_aic.index(min_single_aic)]
    logger.info(
        "Model with smallest AIC for single variable regressions used %i",
        min_single_aic_var,
    )

    logger.info("Run joint (additive interactions only) regression...")
    # Including
    add_formula = f"pop_v_path ~ {' + '.join(variables)}"
    add_model = _run_glm(add_formula)

    logger.info("Running joint regression with all interactions...")
    mult_formula = f"pop_v_path ~ {' * '.join(variables)}"
    mult_model = _run_glm(mult_formula)

    logger.info("Running joint regression with specific interactions...")
    # Currently hardcoded to be formula from ExAC
    spec_formula = "pop_v_path ~ oe + misbad + oe:misbad + polyphen + oe:polyphen"
    spec_model = _run_glm(spec_formula)

    single_var_aic.extend([add_model.aic, mult_model.aic, spec_model.aic])
    min_aic = min(single_var_aic)
    logger.info("Lowest model AIC: %i", min_aic)
    if min_aic == min_single_aic:
        logger.info(
            "Single variable regression using %s had the lowest AIC", min_single_aic_var
        )
        logger.info("Coefficients: %s", single_var_res[min_single_aic_var])
        single_var_res[min_single_aic_var].to_csv(output_fname)
    elif min_aic == add_model.aic:
        logger.info(
            "Joint regression using additive interactions (%s) had the lowest AIC",
            add_formula,
        )
        logger.info("Coefficients: %s", add_model.params)
        add_model.params.to_csv(output_fname)
    elif min_aic == mult_model.aic:
        logger.info(
            "Joint regression using all interactions (%s) had the lowest AIC",
            mult_formula,
        )
        logger.info("Coefficients: %s", mult_model.params)
        mult_model.params.to_csv(output_fname)
    else:
        logger.info(
            "Joint regression using specific interactions (%s) had the lowest AIC",
            spec_formula,
        )
        logger.info("Coefficients: %s", spec_model.params)
        spec_model.params.to_csv(output_fname)
