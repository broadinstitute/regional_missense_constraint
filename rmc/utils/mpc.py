import logging
import pandas as pd
import pickle
import statsmodels
import statsmodels.api as sm
from typing import Dict, List, Tuple, Union

import hail as hl

from gnomad.resources.grch37.gnomad import public_release
from gnomad.resources.grch37.reference_data import vep_context
from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists

from rmc.resources.basics import (
    blosum,
    blosum_txt_path,
    gnomad_fitted_score,
    grantham,
    grantham_txt_path,
    joint_clinvar_gnomad,
    misbad,
    mpc_model_pkl_path,
    temp_path,
)
from rmc.resources.grch37.reference_data import cadd, clinvar_path_mis
from rmc.utils.generic import get_aa_map, process_vep
from rmc.utils.missense_badness import annotate_and_filter_codons, get_oe_annotation


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("mpc_utils")
logger.setLevel(logging.INFO)


def convert_score_list_to_ht(
    score_list: List[Dict[str, Union[str, float]]],
    schema: str = "array<struct{amino_acids: str, score: float}>",
    key_fields: Tuple[str] = ("ref", "alt"),
) -> hl.Table:
    """
    Convert list of amino acid changes/associated scores to Table format.

    :param List[Dict[str, float]] score_list: List of dictionaries. Each dictionary contains two keys:
        amino_acids (value: ref and alt amino acids) and score (value: associated score).
    :param str schema: Schema of `score_list`. Default is 'array<struct{amino_acids: str, score: float}>'.
        Note that the dictionary keys must match the field names provided in this schema
        (amino_acids and score).
    :param str key_fields: Desired key fields for the new Table. Default is ("ref", "alt").
    """
    ht = hl.Table.parallelize(hl.literal(score_list, schema))
    if schema == "array<struct{amino_acids: str, score: float}>":
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
                            {"amino_acids": f"{aa}_{alt_aa}", "score": float(item)}
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
    with hl.hadoop_open(grantham_txt_path) as g:
        for line in g:
            # Grab header line (starts with '.')
            if line.startswith("."):
                header = line.strip().split("\t")
                header_dict = {}
                for counter, item in enumerate(header[1:]):
                    header_dict[counter] = aa_map[item]
            else:
                line = line.strip().split("\t")
                aa = aa_map[line[0]]

                for counter, item in enumerate(line[1:]):
                    alt_aa = header_dict[counter]
                    grantham_scores.append(
                        {"amino_acids": f"{aa}_{alt_aa}", "score": float(item)}
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

    logger.info("Adding CADD...")
    # TODO: Make sure future CADD HTs have already been split
    cadd_ht = cadd.ht()
    cadd_ht = cadd_ht.transmute(
        raw=cadd_ht.RawScore,
        phred=cadd_ht.PHRED,
    )
    ht = ht.annotate(cadd=hl.struct(**cadd_ht[ht.key]))

    logger.info("Getting PolyPhen-2 and codon annotations from VEP context HT...")
    context_ht = vep_context.ht().select_globals().select("vep", "was_split")
    context_ht = context_ht.filter(hl.is_defined(ht[context_ht.key]))
    context_ht = process_vep(context_ht)
    context_ht = context_ht.annotate(
        polyphen=hl.struct(
            prediction=context_ht.transcript_consequences.polyphen_prediction,
            score=context_ht.transcript_consequences.polyphen_score,
        )
    )
    context_ht = context_ht.select(
        "polyphen",
        codons=context_ht.transcript_consequences.codons,
        most_severe_consequence=context_ht.transcript_consequences.most_severe_consequence,
        transcript=context_ht.transcript_consequences.transcript_id,
    )
    context_ht = annotate_and_filter_codons(context_ht)
    context_ht = context_ht.checkpoint(f"{temp_path}/polyphen.ht", overwrite=True)

    logger.info(
        "Adding PolyPhen-2, codon, and transcript annotations to joint ClinVar/gnomAD HT..."
    )
    ht = ht.annotate(**context_ht[ht.key])

    logger.info("Getting regional missense constraint missense o/e annotation...")
    ht = get_oe_annotation(ht)

    logger.info("Getting missense badness annotation...")
    mb_ht = misbad.ht()
    ht = ht.annotate(misbad=mb_ht[ht.ref, ht.alt].misbad)

    logger.info("Adding BLOSUM and Grantham annotations...")
    if not file_exists(blosum.path):
        import_blosum()
    blosum_ht = blosum.ht()
    if not file_exists(grantham.path):
        import_grantham()
    grantham_ht = grantham.ht()
    ht = ht.annotate(
        blosum=blosum_ht[ht.ref, ht.alt].score,
        grantham=grantham_ht[ht.ref, ht.alt].score,
    )

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
    variables: List[str] = ["oe", "misbad", "polyphen"],
    additional_variables: List[str] = ["blosum", "grantham"],
) -> None:
    """
    Run single variable and joint regressions and pick best model.

    These regressions are used to determine the fitted score that is used to predict MPC scores.

    For a variant v:
    Fitted score (from ExAC):
        fitted_score(v) = 4.282793 + (4.359682*v[obs_exp]) + (-3.654815*v[misbad]) + (-3.512215*v[pph2])
                    + (2.585361*v[obs_exp]*v[misbad]) + (1.350056*v[obs_exp]*v[pph2])

    Relationship between fitted score and MPC (from ExAC):
        mpc(v) = -log10(n_less(v))/82932)
        n_less(v) = number of common (MAF > 0.01) ExAC variants with fitted_score < fitted_score(v)

    :param List[str] variables: Variables to include in all regressions (single, joint).
        Default is ["oe", "misbad", "polyphen"].
    :param List[str] additional_variables: Additional variables to include in single variable regressions only.
        Default is ["blosum", "grantham"].
    :return: None; function writes Table with gnomAD fitted scores
        and model coefficients as pickle to resource paths.
    """
    from patsy import dmatrices

    # Convert HT to pandas dataframe as logistic regression aggregations aren't currently possible in hail
    ht = joint_clinvar_gnomad.ht()
    ht = ht.transmute(polyphen=ht.polyphen.score)
    df = ht.to_pandas()

    def _run_glm(
        formula: str,
    ) -> Tuple[
        pd.core.frame.DataFrame,
        statsmodels.genmod.generalized_linear_model.GLMResultsWrapper,
    ]:
        """
        Run logistic regression using input formula and return model results.

        MPC formula (from ExAC):
        `pop_v_path ~ obs_exp + mis_badness3 + obs_exp:mis_badness3 + polyphen2 + obs_exp:polyphen2`
        Missense badness was calculated three times in ExAC.
        The third version (mis_badness3) is the version that was released.

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
        return (X, model)

    logger.info("Run single variable regressions...")
    # E.g., mod.misbad3 <- glm(pop_v_path ~ mis_badness3, data=cleaned_joint_exac_clinvar.scores, family=binomial)
    single_var_res = {}
    single_var_aic = []
    all_var = variables + additional_variables
    for var in all_var:
        logger.info("Running %s regression...", var)
        # Create design matrices
        formula = f"pop_v_path ~ {var}"
        X, model = _run_glm(formula)
        single_var_aic.append(model.aic)
        single_var_res[var] = model.params

    # Find lowest AIC for single variable regressions and corresponding model
    min_single_aic = min(single_var_aic)
    min_single_aic_var = all_var[single_var_aic.index(min_single_aic)]
    logger.info(
        "Model with smallest AIC for single variable regressions used %s",
        min_single_aic_var,
    )

    logger.info("Running joint (additive interactions only) regression...")
    add_formula = f"pop_v_path ~ {' + '.join(variables)}"
    add_X, add_model = _run_glm(add_formula)

    logger.info("Running joint regression with all interactions...")
    mult_formula = f"pop_v_path ~ {' * '.join(variables)}"
    mult_X, mult_model = _run_glm(mult_formula)

    logger.info("Running joint regression with specific interactions...")
    # Currently hardcoded to be formula from ExAC
    spec_formula = "pop_v_path ~ oe + misbad + oe:misbad + polyphen + oe:polyphen"
    spec_X, spec_model = _run_glm(spec_formula)

    all_model_aic = single_var_aic + [add_model.aic, mult_model.aic, spec_model.aic]
    min_aic = min(all_model_aic)
    min_aic = min(single_var_aic)
    logger.info("Lowest model AIC: %i", min_aic)
    if min_aic == min_single_aic:
        logger.info(
            "Single variable regression using %s had the lowest AIC", min_single_aic_var
        )
        logger.info("Coefficients: %s", single_var_res[min_single_aic_var])
        model = single_var_res
    elif min_aic == add_model.aic:
        logger.info(
            "Joint regression using additive interactions (%s) had the lowest AIC",
            add_formula,
        )
        logger.info("Coefficients: %s", add_model.params)
        model = add_model
        X = add_X
    elif min_aic == mult_model.aic:
        logger.info(
            "Joint regression using all interactions (%s) had the lowest AIC",
            mult_formula,
        )
        logger.info("Coefficients: %s", mult_model.params)
        model = mult_model
        X = mult_X
    else:
        logger.info(
            "Joint regression using specific interactions (%s) had the lowest AIC",
            spec_formula,
        )
        logger.info("Coefficients: %s", spec_model.params)
        model = spec_model
        X = spec_X

    logger.info("Annotating gnomAD variants with fitted score...")
    model = spec_model
    X = spec_X
    df["fitted_score"] = model.predict(X)

    logger.info("Converting gnomAD variants dataframe into Table and writing...")
    ht = hl.Table.from_pandas(df)
    ht = ht.filter(ht.pop_v_path == 0)
    ht.write(gnomad_fitted_score.path, overwrite=True)

    logger.info("Saving model as pickle...")
    with hl.hadoop_open(mpc_model_pkl_path, "wb") as p:
        pickle.dump(model, p, protocol=pickle.HIGHEST_PROTOCOL)


def annotate_mpc(
    ht: hl.Table,
    output_path: str,
    n_less_eq0_float: float = 0.83,
    interaction_char: str = ":",
    intercept_str: str = "Intercept",
) -> None:
    """
    Annotate Table with MPC component variables and calculate MPC using relationship defined in `mpc_rel_vars`.

    Relationship in `mpc_rel_vars` is the formula used to calculate a variant's fitted score.
    A variant of interest's fitted score is combined with the number of common (MAF > 0.01) variants in gnomAD with fitted scores < the fitted score for the variant of interest to determine the variant's MPC score.
    to determine a variant's MPC score.

    For more information on the fitted score and MPC calculation, see the docstring of `run_regressions`.

    .. note::
        Assume input Table is keyed by locus and alleles.

    :param hl.Table ht: Input Table to be annotated.
    :param str output_path: Where to write Table after adding MPC annotations.
    :param float n_less_eq0_float: Set `n_less` annotation to this float if `n_less` is 0.
        This avoids errors in the `hl.log10` call and ensures that MPC for variants with a fitted score more severe than any common gnomAD variant (`n_less` = 0) is more severe but by a controlled amount compared to MPC for variants with a fitted score more severe than one common gnomAD variant (`n_less` = 1).
        slightly more severe than variants seen only once in gnomAD.
    :param str interaction_char: Character representing interactions in MPC model. Must be one of "*", ":".
        Default is ":".
    :param str intercept_str: Name of intercept variable in MPC model pickle. Default is "Intercept".
    :return: None; function writes Table to specified output path.
    """
    assert interaction_char in {"*", ":"}, "interaction_char must be one of '*' or ':'!"

    logger.info("Extracting MPC model relationships from pickle...")
    with hl.hadoop_open(mpc_model_pkl_path, "rb") as p:
        model = pickle.load(p)
    mpc_rel_vars = model.params.to_dict()
    try:
        intercept = mpc_rel_vars.pop(intercept_str)
    except KeyError:
        raise DataException(
            f"{intercept_str} not in model parameters! Please double check and rerun."
        )

    logger.info("Annotating HT with MPC variables...")
    if "transcript" not in ht.row:
        # Get transcript annotation from polyphen HT
        # (context HT filtered to contain transcript, ref/alt amino acids, and polyphen annotation)
        polyphen_ht = hl.read_table(f"{temp_path}/polyphen.ht").select("transcript")
        ht = ht.annotate(transcript=polyphen_ht[ht.key].transcript)

    variables = mpc_rel_vars.keys()
    if "oe" in variables:
        logger.info("Getting regional missense constraint missense o/e annotation...")
        ht = get_oe_annotation(ht)

    if "misbad" in variables:
        logger.info("Getting missense badness annotation...")
        if "ref" not in ht.row:
            polyphen_ht = hl.read_table(f"{temp_path}/polyphen.ht").select("ref", "alt")
            ht = ht.annotate(ref=polyphen_ht[ht.key].ref, alt=polyphen_ht[ht.key].alt)
        mb_ht = misbad.ht()
        ht = ht.annotate(misbad=mb_ht[ht.ref, ht.alt].misbad)

    if "polyphen" in variables:
        logger.info("Annotating HT with Polyphen...")
        polyphen_ht = hl.read_table(f"{temp_path}/polyphen.ht").select("polyphen")
        ht = ht.annotate(polyphen=polyphen_ht[ht.key].polyphen.score)

    logger.info("Aggregating gnomAD fitted scores...")
    gnomad_ht = gnomad_fitted_score.ht()
    gnomad_var_count = gnomad_ht.count()
    gnomad_ht = gnomad_ht.group_by("fitted_score").aggregate(n_var=hl.agg.count())
    gnomad_ht = gnomad_ht.order_by("fitted_score")
    gnomad_ht = gnomad_ht.key_by("fitted_score")
    gnomad_ht = gnomad_ht.annotate(n_less=hl.scan.sum(gnomad_ht.n_var))
    #
    gnomad_ht = gnomad_ht.annotate(
        n_less=hl.if_else(
            gnomad_ht.n_less == 0,
            n_less_eq0_float,
            gnomad_ht.n_less,
        )
    )
    gnomad_ht = gnomad_ht.checkpoint(f"{temp_path}/gnomad_group.ht", overwrite=True)

    logger.info("Annotating fitted scores...")
    variable_dict = {
        variable: mpc_rel_vars[variable]
        for variable in mpc_rel_vars
        if variable.isalpha()
    }
    interactions_dict = {
        variable: mpc_rel_vars[variable]
        for variable in variables
        if not variable.isalpha()
    }
    annot_expr = [(ht[variable] * mpc_rel_vars[variable]) for variable in variable_dict.keys()]
    # NOTE: This assumes we won't have more than one variable interacting
    interaction_annot_expr = [
        (
            ht[variable.split(interaction_char)[0]]
            * ht[variable.split(interaction_char)[1]]
            * mpc_rel_vars[variable]
        )
        for variable in interactions_dict.keys()
    ]
    annot_expr.extend(interaction_annot_expr)
    combined_annot_expr = hl.fold(lambda i, j: i + j, 0, annot_expr)

    ht = ht.annotate(fitted_score=coefficient + combined_annot_expr)
    ht = ht.annotate(n_less=gnomad_ht[ht.fitted_score].n_less)
    ht = ht.annotate(mpc=-(hl.log10(ht.n_less / gnomad_var_count)))
    ht.write(output_path)
