import logging
import pickle
from typing import Dict, List, Tuple, Union

import hail as hl
import numpy as np
import pandas as pd
import statsmodels
import statsmodels.api as sm
from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists
from patsy import dmatrices
from itertools import combinations

from rmc.resources.basics import TEMP_PATH_WITH_FAST_DEL
from rmc.resources.reference_data import (
    blosum,
    blosum_txt_path,
    cadd,
    clinvar_plp_mis_haplo,
    FOLD_K,
    grantham,
    grantham_txt_path,
    train_val_test_transcripts_path,
)
from rmc.resources.rmc import (
    CURRENT_FREEZE,
    context_with_oe,
    context_with_oe_dedup,
    gnomad_fitted_score_path,
    joint_clinvar_gnomad_path,
    misbad_path,
    mpc_model_pkl_path,
    mpc_release,
)
from rmc.utils.generic import (
    get_aa_map,
    get_gnomad_public_release,
    keep_criteria,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("mpc_utils")
logger.setLevel(logging.INFO)


def convert_score_list_to_ht(
    score_list: List[Dict[str, Union[str, float]]],
    schema: str = "array<struct{amino_acids: str, score: float}>",
    key_fields: Tuple[str, str] = ("ref", "alt"),
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
    gnomad_data_type: str = "exomes",
    af_threshold: float = 0.001,
    overwrite_temp: bool = False,
    do_k_fold_training: bool = False,
    freeze: int = CURRENT_FREEZE,
    adj_freq_index: int = 0,
    cov_threshold: int = 0,
) -> None:
    """
    Prepare Table(s) with 'population' and 'pathogenic' variants in given set of transcripts.

    'Population' variants are defined as common gnomAD missenses.
    'Pathogenic' variants are defined as ClinVar pathogenic/likely pathogenic missenses in
    predicted-haploinsufficient genes.

    .. note::
        Assumes tables containing all variants in canonical transcripts and their
        missense O/E exist (both duplicated and dedup versions).

    :param str gnomad_data_type: gnomAD data type. Used to retrieve public release Table.
        Must be one of "exomes" or "genomes" (check is done within `public_release`).
        Default is "exomes".
    :param float af_threshold: Allele frequency cutoff to filter gnomAD public dataset.
        Variants at frequencies greater than or equal to this threshold will be kept.
        Default is 0.001.
    :param bool overwrite_temp: Whether to overwrite intermediate temporary data if it already exists.
        If False, will read existing intermediate temporary data rather than overwriting.
        Default is False.
    :param do_k_fold_training: Whether to generate k-fold Tables with the training transcripts.
        If False, will use all training transcripts to prepare a single Table for a single model.
        If True, will prepare Tables for k models corresponding to the k-folds of the training transcripts.
        Default is False.
    :param int freeze: RMC data freeze number. Default is CURRENT_FREEZE.
    :param adj_freq_index: Index of array that contains allele frequency information calculated on
        high quality (adj) genotypes across genetic ancestry groups. Default is 0.
    :param int cov_threshold: Coverage threshold used to filter context Table. Default is 0.
    :return: None; function writes Table(s) to resource path.
    """
    logger.info("Reading in ClinVar P/LP missense variants in severe HI genes...")
    clinvar_ht = clinvar_plp_mis_haplo.ht()
    # TODO: Change severe gene definition? Include TS genes?
    clinvar_ht = clinvar_ht.annotate(pop_v_path=0)

    logger.info(
        "Importing gnomAD public data and filtering to high quality, common variants..."
    )
    gnomad_ht = get_gnomad_public_release(gnomad_data_type, adj_freq_index)
    gnomad_ht = gnomad_ht.filter(
        keep_criteria(
            ac_expr=gnomad_ht.ac,
            af_expr=gnomad_ht.af,
            filters_expr=gnomad_ht.filters,
            cov_expr=gnomad_ht.gnomad_coverage,
            af_threshold=af_threshold,
            cov_threshold=cov_threshold,
            filter_to_rare=False,
        )
    )
    gnomad_ht = gnomad_ht.annotate(pop_v_path=1)

    logger.info("Joining ClinVar and gnomAD HTs and filtering out overlaps...")
    ht = clinvar_ht.select("pop_v_path").union(gnomad_ht.select("pop_v_path"))
    # Remove variants that are in both the benign and pathogenic set
    counts = ht.group_by(*ht.key).aggregate(n=hl.agg.count())
    counts = counts.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/clinvar_gnomad_counts.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )
    overlap = counts.filter(counts.n > 1)
    logger.info("%i ClinVar P/LP variants are common in gnomAD", overlap.count())
    ht = ht.anti_join(overlap)
    ht = ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/joint_clinvar_gnomad.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    logger.info("Adding CADD...")
    # TODO: Make sure future CADD HTs have already been split
    cadd_ht = cadd.ht()
    cadd_ht = cadd_ht.transmute(raw=cadd_ht.RawScore, phred=cadd_ht.PHRED)
    ht = ht.annotate(cadd=hl.struct(**cadd_ht[ht.key]))

    logger.info("Getting transcript annotations...")
    # If OE-annotated VEP context HT filtered to missense variants
    # in canonical transcripts does not exist, create it
    if not file_exists(context_with_oe_dedup.versions[freeze].path) or not file_exists(
        context_with_oe.versions[freeze].path
    ):
        raise DataException("OE-annotated context table does not exist!")

    # Get transcript annotation from deduplicated context HT resource
    context_ht = context_with_oe_dedup.versions[freeze].ht()
    ht = ht.annotate(transcript=context_ht[ht.key].transcript)
    ht = ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/joint_clinvar_gnomad_transcript.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    logger.info(
        "Getting PolyPhen-2, codon, and regional missense constraint missense o/e annotations..."
    )
    context_ht = context_with_oe.versions[freeze].ht()
    ht = ht.annotate(**context_ht[ht.locus, ht.alleles, ht.transcript])
    ht = ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/joint_clinvar_gnomad_context_annot.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

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

    logger.info("Filtering to rows with defined annotations...")
    ht = ht.filter(
        hl.is_defined(ht.cadd.phred)
        & hl.is_defined(ht.blosum)
        & hl.is_defined(ht.grantham)
        & hl.is_defined(ht.oe)
        & hl.is_defined(ht.polyphen.score)
    )
    ht = ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/joint_clinvar_gnomad_annot.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    def _add_misbad_and_filter_transcripts(
        ht: hl.Table,
        fold: int = None,
    ) -> None:
        """
        Add missense badness calculated on a set of transcripts to table, filter to those transcripts, and write out.

        :param hl.Table ht: Input table.
        :param int fold: Fold number in training set to select training transcripts from.
            If not None, the Table is generated from variants in only training transcripts
            from the specified fold of the overall training set. If None, the Table is generated from
            variants in all training transcripts. Default is None.
        :return: None; function writes HT to specified path.
        """
        mb_ht = hl.read_table(misbad_path(fold=fold, freeze=freeze))
        ht = ht.annotate(misbad=mb_ht[ht.ref, ht.alt].misbad)
        filter_transcripts = hl.experimental.read_expression(
            train_val_test_transcripts_path(fold=fold)
        )
        ht = ht.filter(
            hl.is_defined(ht.misbad) & filter_transcripts.contains(ht.transcript)
        )
        ht.write(joint_clinvar_gnomad_path(fold=fold, freeze=freeze), overwrite=True)

    if not do_k_fold_training:
        logger.info(
            "Adding misbad from training transcripts, filtering transcripts, writing out...",
        )
        _add_misbad_and_filter_transcripts(ht=ht)
    else:
        logger.info(
            "Adding misbad, filtering transcripts, writing out for the %i-fold training sets...",
            FOLD_K,
        )
        for i in range(1, FOLD_K + 1):
            _add_misbad_and_filter_transcripts(ht=ht, fold=i)


def run_glm(
    df: pd.core.frame.DataFrame,
    formula: str,
) -> statsmodels.genmod.generalized_linear_model.GLMResultsWrapper:
    """
    Run logistic regression using input formula and return model results.

    For formula reference, see: https://learn-scikit.oneoffcoder.com/patsy.html.

    :param pd.core.frame.DataFrame df: Table holding data to fit model with.
    :param str formula: R-style formula defining model.
    :return: Model object holding logistic regression results.
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


# TODO: Separate out comparing models with AIC and training regressions from all combinations
# of a set of variables
# After merging with model flexibility code changes
def get_min_aic_model(
    df: pd.core.frame.DataFrame,
    variables: List[str],
) -> statsmodels.genmod.generalized_linear_model.GLMResultsWrapper:
    """
    Run single variable and joint regressions using a given set of transcripts and pick best model.

    :param pd.core.frame.DataFrame df: Table holding data to fit model with.
    :param List[str] variables: Variables to include in all regressions (single, joint).
    :return: Logistic regression results.
    """
    model_types = ["single", "additive", "multiplicative", "spec"]
    # Model type labels for use in logging information
    log_model_types = {
        "single": "single variable",
        "additive": "joint (no interactions)",
        "multiplicative": "joint (with interactions)",
        "spec": "specific formula",
    }

    def _find_min_aic_model(
        var_combs: List[Tuple[str, ...]],
        model_type: str = "single",
    ) -> Tuple[statsmodels.genmod.generalized_linear_model.GLMResultsWrapper, str,]:
        """
        Run logistic regressions on different variable combinations and return results on model with lowest AIC.

        :param List[Tuple[str]] var_combs: Variable combinations to generate regressions with.
        :param str model_type: Regression formula type.
            One of "single", "additive", or "multiplicative".
            If "single", each variable combination must consist of only one variable.
            If "additive", no interactions will be included.
            If "multiplicative", all interaction terms will be included.
            Default is "single".
        :return: Tuple of logistic regression model results for variable combination with lowest AIC.
            Tuple is composed of the model object and the regression formula.
        """
        if model_type not in {"single", "additive", "multiplicative"}:
            raise DataException(
                "Model type must be one of 'single', 'additive', or 'multiplicative'!"
            )
        if (model_type == "single") and (not all([len(x) == 1 for x in var_combs])):
            raise DataException(
                "Please check that all single variable inputs consist of only one variable!"
            )

        models = {}
        model_aics = []
        model_formulas = []
        formula_interaction_char = "*" if model_type == "multiplicative" else "+"
        log_model_type = log_model_types[model_type]

        for var_comb in var_combs:
            logger.info("Running %s regression for %s...", log_model_type, var_comb)
            formula = f"pop_v_path ~ {formula_interaction_char.join(var_comb)}"
            model = run_glm(df, formula)
            models[var_comb] = model
            model_aics.append(model.aic)
            model_formulas.append(formula)

        # Find model with lowest AIC
        min_model_aic = min(model_aics)
        min_model_idx = model_aics.index(min_model_aic)
        min_model_var = var_combs[min_model_idx]
        min_model = models[min_model_var]
        min_model_formula = model_formulas[min_model_idx]
        if model_aics.count(min_model_aic) > 1:
            logger.warning(
                """
                There is a tie for minimum AIC in the %s regressions.
                This function will use the first variable combination it finds by default!
                """,
                log_model_type,
            )
        logger.info(
            "Model with smallest AIC for %s regressions used variables: %s",
            log_model_type,
            min_model_var,
        )

        return (min_model, min_model_formula)

    min_models = {}
    min_model_formulas = {}
    min_model_aics = {}

    # Find model with lowest AIC for single variable regressions
    min_models["single"], min_model_formulas["single"] = _find_min_aic_model(
        [(x,) for x in variables]
    )

    # List possible variable combinations for joint regressions
    var_combs = []
    for k in range(2, len(variables) + 1):
        var_combs += [x for x in combinations(variables, k)]
    # Find models with lowest AICs for joint variable regressions (with or without interactions)
    for model_type in ["additive", "multiplicative"]:
        min_models[model_type], min_model_formulas[model_type] = _find_min_aic_model(
            var_combs, model_type
        )

    logger.info("Running joint regression with specific interactions...")
    # Currently hardcoded to be formula from ExAC
    min_model_formulas[
        "spec"
    ] = "pop_v_path ~ oe + misbad + oe:misbad + polyphen + oe:polyphen"
    min_models["spec"] = run_glm(df, min_model_formulas["spec"])

    min_model_aics = {x: min_models[x].aic for x in model_types}
    overall_min_aic = min(min_model_aics.values())
    logger.info("Lowest model AIC: %f", overall_min_aic)
    if list(min_model_aics.values()).count(overall_min_aic) > 1:
        logger.warning(
            """
            There is a tie for minimum AIC over model types.
            This function will use the first model it finds by default
            (single variable -> no interactions -> all interactions -> specific interactions)!
            """
        )
    overall_min_model_type = list(min_model_aics.keys())[
        list(min_model_aics.values()).index(overall_min_aic)
    ]
    overall_min_model = min_models[overall_min_model_type]
    logger.info(
        "Model with overall lowest AIC is %s regression (%s)",
        overall_min_model_type,
        min_model_formulas[overall_min_model_type],
    )
    logger.info("Coefficients: %s", overall_min_model.params)

    return overall_min_model


def run_regressions(
    use_model_formula: bool = False,
    model_formula: str = "pop_v_path ~ oe + misbad + oe:misbad + polyphen + oe:polyphen",
    variables: List[str] = ["oe", "misbad", "polyphen"],
    additional_variables: List[str] = ["blosum", "grantham"],
    overwrite_temp: bool = False,
    do_k_fold_training: bool = False,
    freeze: int = CURRENT_FREEZE,
) -> None:
    """
    Train regression model(s) on given set of transcripts to determine fitted scores for MPC calculation.

    Depending on the `use_model_formula` parameter, this function can be used to:
    (1) Train a model with a specified formula, or
    (2) Train single variable and joint regressions on combinations of the input variables
    and pick the best model out of these based on AIC.

    Regression formula from ExAC:
        `pop_v_path ~ obs_exp + mis_badness3 + obs_exp:mis_badness3 + polyphen2 + obs_exp:polyphen2`
    Missense badness was calculated three times in ExAC.
    The third version (mis_badness3) is the version that was released.

    Fitted score for a variant v from ExAC:
        fitted_score(v) = 4.282793 + (4.359682*v[obs_exp]) + (-3.654815*v[misbad]) + (-3.512215*v[pph2])
                        + (2.585361*v[obs_exp]*v[misbad]) + (1.350056*v[obs_exp]*v[pph2])

    Relationship between fitted score and MPC from ExAC:
        mpc(v) = -log10(n_less(v))/82932)
        n_less(v) = number of common (AC > 121, AF > 0.001) ExAC variants with fitted_score < fitted_score(v)

    Note that higher MPC scores indicate increased predicted missense deleteriousness, and
    smaller n_less values and fitted scores correspond to higher MPC scores.

    :param bool use_model_formula: Whether to use a specified model formula.
        If True, only one model will be built using the formula specified in `model_formula`.
        If False, all possible combinations of variables in `variables` and `additional_variables`
        will be used to build models. The single best model based on AIC will be kept.
        Default is False.
    :param str model_formula: Model formula to use in regression.
        Default is 'pop_v_path ~ oe + misbad + oe:misbad + polyphen + oe:polyphen'.
        This was the formula selected in ExAC.
    :param List[str] variables: Variables to include in all regressions (single, joint).
        Default is ['oe', 'misbad', 'polyphen'].
    :param List[str] additional_variables: Additional variables to include in single variable regressions only.
        Default is ['blosum', 'grantham'].
    :param bool overwrite_temp: Whether to overwrite intermediate temporary data if it already exists.
        If False, will read existing intermediate temporary data rather than overwriting.
        Default is False.
    :param bool do_k_fold_training: Whether to generate k-fold models with the training transcripts.
        If False, will use all training transcripts in calculation of a single model.
        If True, will calculate k models corresponding to the k-folds of the training transcripts.
        Default is False.
    :param int freeze: RMC data freeze number. Default is CURRENT_FREEZE.
    :return: None; function writes Table(s) with gnomAD fitted scores and pickled model(s) to resource paths.
    """
    if not use_model_formula and do_k_fold_training:
        raise DataException(
            "Training separate models with k folds may result in different optimums!"
        )

    def _generate_regression_model(
        temp_label: str,
        fold: int = None,
    ):
        """
        Generate a regression model on the given set of transcripts.

        :param str temp_label: Suffix to add to temporary data paths to avoid conflicting names for different models.
        :param int fold: Fold number in training set to select training transcripts from.
            If not None, the Table is generated from variants in only training transcripts
            from the specified fold of the overall training set. If None, the Table is generated from
            variants in all training transcripts. Default is None.
        :return: None; function writes Table with gnomAD fitted scores and pickled model to resource paths.
        """
        # Convert HT to pandas dataframe as logistic regression aggregations aren't currently possible in hail
        ht = hl.read_table(joint_clinvar_gnomad_path(fold=fold, freeze=freeze))
        ht = ht.transmute(polyphen=ht.polyphen.score)
        df = ht.to_pandas()

        # NOTE: `Table.to_pandas()` outputs a DataFrame with nullable dtypes
        # for int and float fields, which patsy cannot handle
        # These must be converted to non-nullable standard numpy dtypes
        # prior to working with patsy
        # See: https://github.com/hail-is/hail/commit/da557655ef1da99ddd1887bd32b33f4b5adcec9f
        for column in df.columns:
            if df[column].dtype == "Int32":
                df[column] = df[column].astype(np.int32)
            if df[column].dtype == "Int64":
                df[column] = df[column].astype(np.int64)
            if df[column].dtype == "Float32":
                df[column] = df[column].astype(np.float32)
            if df[column].dtype == "Float64":
                df[column] = df[column].astype(np.float64)

        if use_model_formula:
            import re

            # Check input formula is formulated properly
            formula_components = re.split("[~+:*]", model_formula.replace(" ", ""))
            formula_y = formula_components.pop(0)
            formula_xs = set(formula_components)
            if (formula_y != "pop_v_path") or (
                not all([x in df.columns.values for x in formula_xs])
            ):
                raise DataException("Model formula is not formulated properly!")

            model = run_glm(df, model_formula)
        else:
            all_variables = variables + additional_variables
            model = get_min_aic_model(df, all_variables)

        logger.info("Saving model as pickle...")
        with hl.hadoop_open(mpc_model_pkl_path(fold=fold, freeze=freeze), "wb") as p:
            pickle.dump(model, p, protocol=pickle.HIGHEST_PROTOCOL)

        logger.info(
            "Annotating gnomAD variants with fitted score and writing to gnomad_fitted_score path..."
        )
        ht = ht.filter(ht.pop_v_path == 1)
        ht = calculate_fitted_scores(
            ht=ht,
            temp_label=temp_label,
            overwrite_temp=overwrite_temp,
            fold=fold,
            freeze=freeze,
        )
        ht.write(gnomad_fitted_score_path(fold=fold, freeze=freeze), overwrite=True)

    if not do_k_fold_training:
        logger.info(
            "Running regression(s) on training transcripts...",
        )
        _generate_regression_model(temp_label="_train")
    else:
        logger.info(
            "Running regressions on the %i-fold training sets...",
            FOLD_K,
        )
        for i in range(1, FOLD_K + 1):
            _generate_regression_model(temp_label=f"_train_fold{i}", fold=i)
    # TODO: Make this show/save a table of AICs and model formulas instead of printing every one


def calculate_fitted_scores(
    ht: hl.Table,
    temp_label: str,
    overwrite_temp: bool = False,
    fold: int = None,
    interaction_char: str = ":",
    intercept_str: str = "Intercept",
    freeze: int = CURRENT_FREEZE,
) -> hl.Table:
    """
    Use MPC model chosen in `run_regressions` to calculate fitted scores for input Table.

    .. note::
        - Function will remove any rows with undefined MPC feature annotations from input Table.
        - Input Table is assumed to be keyed by ['locus', 'alleles'].
        - If the model uses transcript-specific features (e.g. RMC O/E) and there are variants
            in the input HT in multiple transcripts, only the most severe fitted score
            (and the corresponding annotations it was calculated from) per variant
            will be retained.

    :param hl.Table ht: Input Table with variants to be annotated.
    :param str temp_label: Suffix to add to temporary data paths to avoid conflicting names for different models.
    :param bool overwrite_temp: Whether to overwrite intermediate temporary data if it already exists.
        If False, will read existing intermediate temporary data rather than overwriting.
        Default is False.
    :param int fold: Fold number in training set to select training transcripts from.
        If not None, the Table is generated from variants in only training transcripts
        from the specified fold of the overall training set. If None, the Table is generated from
        variants in all training transcripts. Default is None.
    :param str interaction_char: Character representing interactions in MPC model. Must be one of "*", ":".
        Default is ":".
    :param str intercept_str: Name of intercept variable in MPC model pickle. Default is "Intercept".
    :param int freeze: RMC data freeze number. Default is CURRENT_FREEZE.
    :return: Table annotated with MPC model features and fitted scores, filtered to variants with defined annotations.
    """
    assert interaction_char in {"*", ":"}, "interaction_char must be one of '*' or ':'!"

    logger.info("Extracting MPC model relationships from pickle...")
    with hl.hadoop_open(mpc_model_pkl_path(fold=fold, freeze=freeze), "rb") as p:
        model = pickle.load(p)
    mpc_rel_vars = model.params.to_dict()
    try:
        intercept = mpc_rel_vars.pop(intercept_str)
    except KeyError:
        raise DataException(
            f"{intercept_str} not in model parameters! Please double check and rerun."
        )

    variable_dict = {k: v for k, v in mpc_rel_vars.items() if not interaction_char in k}
    variables = variable_dict.keys()
    interactions_dict = {k: v for k, v in mpc_rel_vars.items() if interaction_char in k}

    # Select fields from context HT containing annotations to extract
    # NOTE: `ref` and `alt` fields are also extracted for annotation of missense badness
    context_ht_annots = {"oe", "polyphen"}
    annots = context_ht_annots.intersection(variables)
    logger.info("Annotating HT with MPC variables from context HT (%s)...", annots)
    context_ht = context_with_oe.versions[freeze].ht()
    # Re-key context HT by locus and alleles (original key contains transcript)
    # to enable addition of annotations from each transcript corresponding to a variant
    context_ht = context_ht.key_by("locus", "alleles").select(
        *({"transcript", "ref", "alt"} | annots)
    )
    if "polyphen" in annots:
        context_ht = context_ht.transmute(polyphen=context_ht.polyphen.score)
    # Add annotations from all transcripts per variant
    # This creates separate rows for each transcript a variant is in
    scores_ht = context_ht.join(ht.select())
    scores_ht = scores_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/fitted_score_annots_context{temp_label}.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    # Add annotations not in context HT
    if "misbad" in variables:
        logger.info("Annotating HT with missense badness...")
        mb_ht = hl.read_table(misbad_path(fold=fold, freeze=freeze))
        scores_ht = scores_ht.annotate(
            misbad=mb_ht[scores_ht.ref, scores_ht.alt].misbad
        )
        annots.add("misbad")

    if "blosum" in variables:
        logger.info("Annotating HT with BLOSUM...")
        if not file_exists(blosum.path):
            import_blosum()
        blosum_ht = blosum.ht()
        scores_ht = scores_ht.annotate(
            blosum=blosum_ht[scores_ht.ref, scores_ht.alt].score
        )
        annots.add("blosum")

    if "grantham" in variables:
        logger.info("Annotating HT with Grantham...")
        if not file_exists(grantham.path):
            import_grantham()
        grantham_ht = grantham.ht()
        scores_ht = scores_ht.annotate(
            grantham=grantham_ht[scores_ht.ref, scores_ht.alt].score
        )
        annots.add("grantham")

    if annots != variables:
        raise DataException("Check that all MPC model features are accounted for!")

    logger.info("Filtering to defined annotations...")
    filter_expr = True
    for annot in annots:
        filter_expr &= hl.is_defined(scores_ht[annot])
    scores_ht = scores_ht.filter(filter_expr)
    # Checkpoint here to force joins and filter to complete
    scores_ht = scores_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/fitted_score_annots_filtered{temp_label}.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    logger.info("Calculating fitted scores...")
    annot_expr = [
        (scores_ht[variable] * mpc_rel_vars[variable]) for variable in variables
    ]
    interaction_annot_expr = []
    for interaction in interactions_dict.keys():
        expr = mpc_rel_vars[interaction]
        for variable in interaction.split(interaction_char):
            expr *= scores_ht[variable]
        interaction_annot_expr.append(expr)

    annot_expr.extend(interaction_annot_expr)
    combined_annot_expr = hl.fold(lambda i, j: i + j, 0, annot_expr)

    logger.info("Computing fitted scores...")
    scores_ht = scores_ht.annotate(fitted_score=intercept + combined_annot_expr)
    scores_ht = scores_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/fitted_scores_dup{temp_label}.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )

    # Deduplicate variants in multiple transcripts by retaining only the most severe score
    # per variant and corresponding annotations for that transcript
    logger.info("Filtering to most severe fitted score for each variant...")
    min_scores_ht = scores_ht.select("transcript", "fitted_score")
    min_scores_ht = min_scores_ht.collect_by_key()
    min_scores_ht = min_scores_ht.annotate(
        fitted_score=hl.nanmin(min_scores_ht.values["fitted_score"])
    )
    min_scores_ht = min_scores_ht.annotate(
        transcript=min_scores_ht.values.find(
            lambda x: x["fitted_score"] == min_scores_ht["fitted_score"]
        ).transcript
    )
    min_scores_ht = min_scores_ht.drop("values")
    min_scores_ht = min_scores_ht.key_by("locus", "alleles", "transcript")
    scores_ht = scores_ht.key_by("locus", "alleles", "transcript").join(min_scores_ht)
    scores_ht = scores_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/fitted_scores_dedup{temp_label}.ht",
        _read_if_exists=not overwrite_temp,
        overwrite=overwrite_temp,
    )
    logger.info("Annotating fitted score and model features onto input HT...")
    return ht.join(scores_ht)


def aggregate_gnomad_fitted_scores(
    n_less_eq0_float: float = 0.83,
    fold: int = None,
    freeze: int = CURRENT_FREEZE,
) -> None:
    """
    Aggregate gnomAD fitted scores to count number of variants with a score less than a given score.

    :param float n_less_eq0_float: Set `n_less` annotation to this float if `n_less` is 0.
        This avoids errors in the `hl.log10` call and ensures that MPC for variants with a fitted score
        more severe than any common gnomAD variant score (`n_less` = 0) is more severe (by a controlled amount)
        compared to MPC for variants with a fitted score more severe than one common gnomAD variant (`n_less` = 1).
    :param int fold: Fold number in training set to select training transcripts from.
        If not None, the Table is generated from variants in only training transcripts
        from the specified fold of the overall training set. If None, the Table is generated from
        variants in all training transcripts. Default is None.
    :param int freeze: RMC data freeze number. Default is CURRENT_FREEZE.
    :return: None; function writes Table to resource path.
    """
    logger.info("Aggregating gnomAD fitted scores...")
    gnomad_ht = hl.read_table(gnomad_fitted_score_path(fold=fold, freeze=freeze))
    gnomad_ht = gnomad_ht.group_by("fitted_score").aggregate(n_var=hl.agg.count())
    gnomad_ht = gnomad_ht.order_by("fitted_score")
    gnomad_ht = gnomad_ht.key_by("fitted_score")
    gnomad_ht = gnomad_ht.annotate(n_less=hl.scan.sum(gnomad_ht.n_var))
    # Make n_less a non-zero value if it is zero
    gnomad_ht = gnomad_ht.annotate(
        n_less=hl.if_else(
            gnomad_ht.n_less == 0,
            n_less_eq0_float,
            gnomad_ht.n_less,
        )
    )
    # Add index annotation to table and convert from int64
    # (default value returned by `add_index`) to int32
    # This is necessary for code in `annotate_mpc` downstream
    # (`annotate_mpc will try to join this index field with an int32 field;
    # `hl.binary_search` returns an int32 by default)
    gnomad_ht = gnomad_ht.add_index()
    gnomad_ht = gnomad_ht.annotate(idx=hl.int(gnomad_ht.idx))
    gnomad_ht = gnomad_ht.key_by("idx")
    gnomad_ht.write(
        gnomad_fitted_score_path(is_grouped=True, fold=fold, freeze=freeze),
        overwrite=True,
    )


def annotate_mpc(
    ht: hl.Table,
    output_ht_path: str,
    temp_label: str,
    use_release: bool = True,
    overwrite_temp: bool = False,
    model_train_fold: int = None,
    freeze: int = CURRENT_FREEZE,
) -> None:
    # TODO: Should this be `fold` or `all_k_folds`?
    # TODO: Check how many variants overlap between training + test transcripts
    """
    Annotate input Table with MPC scores.

    Depending on the `use_release` parameter, this function can be used to:
    (1) Annotate with scores from the MPC release HT (VEP context HT with MPC scores
    calculated from all training data), or
    (2) Annotate with MPC scores calculated directly using a given MPC model
    (e.g., if saving scores for the whole VEP context HT is not desired).

    To calculate the MPC score, a variant of interest's fitted score from the MPC model
    is combined with the number of common (AF > 0.001) variants in gnomAD with fitted scores
    < the fitted score for the variant of interest.

    For more information on the fitted score and MPC calculation, see the docstring of `run_regressions`.

    :param hl.Table ht: Input Table.
    :param str output_ht_path: Path to write out annotated table (with duplicate variants removed.
    :param str temp_label: Suffix to add to temporary data paths to avoid conflicting names for different models.
    :param bool use_release: Whether to use the MPC HT release in annotation of MPC scores.
        If True, MPC scores will be taken from the MPC release HT. If False, MPC scores
        will be computed using the regression model saved for the given fold of the training set.
        Default is True.
    :param bool overwrite_temp: Whether to overwrite intermediate temporary data if it already exists.
        If False, will read existing intermediate temporary data rather than overwriting.
        Default is False.
    :param int model_train_fold: Fold number in training set used to generate MPC model.
        If not None, MPC model used to generate annotated scores is trained on variants
        in training transcripts from the specified fold of the overall training set.
        If None, model is trained on variants in all training transcripts.
        Default is None.
    :param int freeze: RMC data freeze number. Default is CURRENT_FREEZE.
    :return: None; function writes Table to resource path.
    """
    assert list(ht.key) == [
        "locus",
        "alleles",
    ], "HT must have key ['locus', 'alleles']!"
    # NOTE: hl.tlocus() will use the reference genome build from the hail initialization
    # if the user doesn't specify a reference build
    # This means that for v2, where hail is initialized with GRCh37, this will check for
    # a build 37 locus
    assert ht.key.locus.dtype == hl.tlocus() and ht.key.alleles.dtype == hl.tarray(
        hl.tstr
    ), "'locus' must be a LocusExpression, and 'alleles' must be an array of strings!"

    if use_release:
        if model_train_fold is not None:
            raise DataException(
                "MPC release is generated only from all training transcripts!"
            )
        if not file_exists(mpc_release.versions[freeze].path):
            raise DataException("MPC release has not yet been generated!")

        logger.info("Annotating MPC scores using MPC release HT and writing out...")
        mpc_ht = mpc_release.versions[freeze].ht()
        ht = ht.annotate(mpc=mpc_ht[ht.locus, ht.alleles].mpc)
        ht.write(output_ht_path, overwrite=True)
    else:
        logger.info("Calculating fitted scores on input variants...")
        scores_ht = calculate_fitted_scores(
            ht=ht.select(),
            temp_label=temp_label,
            overwrite_temp=overwrite_temp,
            freeze=freeze,
        )

        logger.info("Aggregating fitted scores for gnomAD common variants...")
        # TODO: Add model labels in paths - separate folder for each model
        fitted_group_path = gnomad_fitted_score_path(
            is_grouped=True, fold=model_train_fold, freeze=freeze
        )
        if not file_exists(fitted_group_path):
            aggregate_gnomad_fitted_scores(fold=model_train_fold, freeze=freeze)
        gnomad_ht = hl.read_table(fitted_group_path)
        gnomad_scores = gnomad_ht.aggregate(
            hl.sorted(hl.agg.collect(gnomad_ht.fitted_score))
        )
        gnomad_scores_len = len(gnomad_scores)

        # Get total number of gnomAD common variants
        gnomad_var_count = hl.read_table(
            gnomad_fitted_score_path(fold=model_train_fold, freeze=freeze)
        ).count()

        logger.info("Getting n_less values for input variants...")
        # Annotate HT with sorted array of gnomAD fitted scores
        scores_ht = scores_ht.annotate_globals(gnomad_scores=gnomad_scores)
        # Checkpoint here to force the gnomAD join to complete
        scores_ht = scores_ht.checkpoint(
            f"{TEMP_PATH_WITH_FAST_DEL}/mpc_temp_gnomad{temp_label}.ht",
            _read_if_exists=not overwrite_temp,
            overwrite=overwrite_temp,
        )

        # Search all gnomAD scores to find first score that is
        # greater than or equal to score to be annotated
        # `binary_search` will return the index of the first gnomAD fitted score that
        # is >= the score of interest
        # e.g., if the score of interest is 0.45, and gnomAD fitted scores are
        # [0.3, 0.4, 0.5], then `binary_search` will return an index of 2
        # the `n_less` of 0.5 will contain the counts of variants with gnomAD scores of
        # 0.3 and 0.4 due to the non-inclusive nature of scans
        # (n_less[0.5] = n_var[0.3] + n_var[0.4])
        scores_ht = scores_ht.annotate(
            idx=hl.binary_search(scores_ht.gnomad_scores, scores_ht.fitted_score)
        )
        scores_ht = scores_ht.annotate(
            n_less=hl.if_else(
                # Make n_less equal to total gnomAD common variant count if
                # index is equal to the length of the gnomAD scores array
                scores_ht.idx == gnomad_scores_len,
                gnomad_var_count,
                gnomad_ht[scores_ht.idx].n_less,
            )
        )
        # Checkpoint here to force the binary search to compute
        scores_ht = scores_ht.checkpoint(
            f"{TEMP_PATH_WITH_FAST_DEL}/mpc_temp_binary{temp_label}.ht",
            _read_if_exists=not overwrite_temp,
            overwrite=overwrite_temp,
        )

        logger.info("Calculating MPC scores...")
        scores_ht = scores_ht.annotate(
            mpc=-(hl.log10(scores_ht.n_less / gnomad_var_count))
        )
        scores_ht = scores_ht.checkpoint(
            f"{TEMP_PATH_WITH_FAST_DEL}/mpc{temp_label}.ht",
            _read_if_exists=not overwrite_temp,
            overwrite=overwrite_temp,
        )

        logger.info("Annotating MPC scores to input table and writing out...")
        ht = ht.annotate(mpc=scores_ht[ht.key].mpc)
        ht.write(output_ht_path, overwrite=True)


def create_mpc_release_ht(
    overwrite_temp: bool = True,
    freeze: int = CURRENT_FREEZE,
) -> None:
    """
    Annotate variants in VEP context Table with MPC scores calculated using all training transcripts.

    :param bool overwrite_temp: Whether to overwrite intermediate temporary data if it already exists.
        If False, will read existing intermediate temporary data rather than overwriting.
        Default is True.
    :param int freeze: RMC data freeze number. Default is CURRENT_FREEZE.
    :return: None; function writes Table to resource path.
    """
    annotate_mpc(
        ht=context_with_oe.versions[freeze].ht().select(),
        output_ht_path=mpc_release.versions[freeze].path,
        temp_label="_release",
        use_release=False,
        overwrite_temp=overwrite_temp,
        freeze=freeze,
    )
