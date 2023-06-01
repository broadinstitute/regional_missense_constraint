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
    overwrite_output: bool = True,
    do_k_fold_training: bool = False,
    freeze: int = CURRENT_FREEZE,
    adj_freq_index: int = 0,
    cov_threshold: int = 0,
) -> None:
    """
    Prepare Table with 'population' (common gnomAD missense) and 'pathogenic' (ClinVar pathogenic/likely pathogenic missense) variants.

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
    :param bool overwrite_output: Whether to entirely overwrite final output data if it already exists.
        If False, will read and modify existing output data by adding or modifying columns rather than overwriting entirely.
        If True, will clear existing output data and write new output data.
        The output Tables are the OE-annotated context Tables with duplicated or deduplicated sections
        and the population/pathogenic variant Table.
        Default is True.
    :param do_k_fold_training: Whether to generate k-fold models with the training transcripts.
        If False, will use all training transcripts in calculation of a single model.
        If True, will calculate k models corresponding to the k-folds of the training transcripts.
        Default is False.
    :param int freeze: RMC data freeze number. Default is CURRENT_FREEZE.
    :param adj_freq_index: Index of array that contains allele frequency information calculated on
        high quality (adj) genotypes across genetic ancestry groups. Default is 0.
    :param cov_threshold: Coverage threshold used to filter context Table. Default is 0.
    :return: None; function writes Table to resource path.
    """
    logger.info("Reading in ClinVar P/LP missense variants in severe HI genes...")
    clinvar_ht = clinvar_plp_mis_haplo.ht()
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
    cadd_ht = cadd_ht.transmute(
        raw=cadd_ht.RawScore,
        phred=cadd_ht.PHRED,
    )
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
        mb_ht = hl.read_table(
            misbad_path(
                fold=fold,
                freeze=freeze,
            )
        )
        ht = ht.annotate(misbad=mb_ht[ht.ref, ht.alt].misbad)
        filter_transcripts = hl.experimental.read_expression(
            train_val_test_transcripts_path(fold=fold)
        )
        ht = ht.filter(
            hl.is_defined(ht.misbad) & filter_transcripts.contains(ht.transcript)
        )
        ht.write(
            joint_clinvar_gnomad_path(
                fold=fold,
                freeze=freeze,
            ),
            overwrite=overwrite_output,
        )

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


def run_regressions(
    variables: List[str] = ["oe", "misbad", "polyphen"],
    additional_variables: List[str] = ["blosum", "grantham"],
    overwrite: bool = True,
    do_k_fold_training: bool = False,
    freeze: int = CURRENT_FREEZE,
) -> None:
    """
    Run single variable and joint regressions and pick best model using training data.

    These regressions are used to determine the fitted score that is used to predict MPC scores.

    For a variant v:
    Fitted score (from ExAC):
        fitted_score(v) = 4.282793 + (4.359682*v[obs_exp]) + (-3.654815*v[misbad]) + (-3.512215*v[pph2])
                    + (2.585361*v[obs_exp]*v[misbad]) + (1.350056*v[obs_exp]*v[pph2])

    Relationship between fitted score and MPC (from ExAC):
        mpc(v) = -log10(n_less(v))/82932)
        n_less(v) = number of common (AC > 121, AF > 0.001) ExAC variants with fitted_score < fitted_score(v)

    Note that higher MPC scores predict increased missense deleteriousness, and
    smaller n_less values and fitted scores will lead to higher MPC scores.

    :param List[str] variables: Variables to include in all regressions (single, joint).
        Default is ["oe", "misbad", "polyphen"].
    :param List[str] additional_variables: Additional variables to include in single variable regressions only.
        Default is ["blosum", "grantham"].
    :param bool overwrite: Whether to overwrite gnomAD fitted score table if it already exists. Default is True.
    :param do_k_fold_training: Whether to generate k-fold models with the training transcripts.
        If False, will use all training transcripts in calculation of a single model.
        If True, will calculate k models corresponding to the k-folds of the training transcripts.
        Default is False.
    :param int freeze: RMC data freeze number. Default is CURRENT_FREEZE.
    :return: None; function writes Table with gnomAD fitted scores
        and model coefficients as pickle to resource paths.
    """

    def _run_glm(
        df: pd.core.frame.DataFrame,
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

        :param pd.core.frame.DataFrame df: Table holding data to fit model with.
        :param str formula: R-style formula defining model.
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

    def _run_regressions_on_transcripts(
        temp_label: str,
        fold: int = None,
    ) -> None:
        """
            Run single variable and joint regressions and pick best model using a given set of transcripts.

            :param str temp_label: Suffix to add to temporary data paths to avoid conflicting names for different models.
        :param int fold: Fold number in training set to select training transcripts from.
            If not None, the Table is generated from variants in only training transcripts
            from the specified fold of the overall training set. If None, the Table is generated from
            variants in all training transcripts. Default is None.
            :return: None; function writes HT to specified path.
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

        logger.info("Run single variable regressions...")
        # E.g., mod.misbad3 <- glm(pop_v_path ~ mis_badness3, data=cleaned_joint_exac_clinvar.scores, family=binomial)
        single_var_res = {}
        single_var_aic = []
        single_var_X = []
        all_var = variables + additional_variables
        for var in all_var:
            logger.info("Running %s regression...", var)
            # Create design matrices
            formula = f"pop_v_path ~ {var}"
            X, model = _run_glm(df=df, formula=formula)
            single_var_X.append(X)
            single_var_aic.append(model.aic)
            single_var_res[var] = model.params

        # Find lowest AIC for single variable regressions and corresponding model
        min_single_aic = min(single_var_aic)
        single_var_idx = single_var_aic.index(min_single_aic)
        min_single_aic_var = all_var[single_var_idx]
        min_single_X = single_var_X[single_var_idx]
        if single_var_aic.count(min_single_aic) > 1:
            logger.warning(
                """
                There is a tie for minimum AIC.
                This function will use the first variable it finds by default!
                """
            )
        logger.info(
            "Model with smallest AIC for single variable regressions used %s",
            min_single_aic_var,
        )

        logger.info("Running joint (additive interactions only) regression...")
        add_formula = f"pop_v_path ~ {' + '.join(variables)}"
        add_X, add_model = _run_glm(df=df, formula=add_formula)

        logger.info("Running joint regression with all interactions...")
        mult_formula = f"pop_v_path ~ {' * '.join(variables)}"
        mult_X, mult_model = _run_glm(df=df, formula=mult_formula)

        logger.info("Running joint regression with specific interactions...")
        # Currently hardcoded to be formula from ExAC
        spec_formula = "pop_v_path ~ oe + misbad + oe:misbad + polyphen + oe:polyphen"
        spec_X, spec_model = _run_glm(df=df, formula=spec_formula)

        all_model_aic = single_var_aic + [add_model.aic, mult_model.aic, spec_model.aic]
        min_aic = min(all_model_aic)
        logger.info("Lowest model AIC: %f", min_aic)
        if all_model_aic.count(min_aic) > 1:
            logger.warning(
                """
                There is a tie for minimum AIC.
                This function will use the first model it finds by default
                (single variable -> additive interactions -> all interactions -> specific interactions)!
                """
            )
        if min_aic == min_single_aic:
            logger.info(
                "Single variable regression using %s had the lowest AIC",
                min_single_aic_var,
            )
            logger.info("Coefficients: %s", single_var_res[min_single_aic_var])
            model = single_var_res
            X = min_single_X
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

        logger.info("Saving model as pickle...")
        with hl.hadoop_open(mpc_model_pkl_path(fold=fold, freeze=freeze), "wb") as p:
            pickle.dump(model, p, protocol=pickle.HIGHEST_PROTOCOL)

        logger.info(
            "Annotating gnomAD variants with fitted score and writing to gnomad_fitted_score path..."
        )
        ht = ht.filter(ht.pop_v_path == 1)
        ht = calculate_fitted_scores(
            ht=ht, temp_label=temp_label, fold=fold, freeze=freeze
        )
        ht.write(
            gnomad_fitted_score_path(fold=fold, freeze=freeze),
            overwrite=overwrite,
        )

    if not do_k_fold_training:
        logger.info(
            "Running regressions on training transcripts...",
        )
        _run_regressions_on_transcripts(temp_label="_train")
    else:
        logger.info(
            "Adding misbad, filtering transcripts, writing out for the %i-fold training sets...",
            FOLD_K,
        )
        for i in range(1, FOLD_K + 1):
            _run_regressions_on_transcripts(
                temp_label=f"_train_fold{i}",
                fold=i,
            )


def calculate_fitted_scores(
    ht: hl.Table,
    temp_label: str = "",
    overwrite_temp: bool = True,
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
        Default is the empty string.
    :param overwrite_temp: Whether to overwrite intermediate temporary data if it already exists.
        If False, will read existing intermediate temporary data rather than overwriting. Default is True.
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
        f"{TEMP_PATH_WITH_FAST_DEL}/fitted_score_annots_context{temp_label}.ht"
    )

    # Add annotations not in context HT
    if "misbad" in variables:
        logger.info("Annotating HT with missense badness...")
        mb_ht = hl.read_table(misbad_path(fold=fold, freeze=freeze))
        scores_ht = scores_ht.annotate(
            misbad=mb_ht[scores_ht.ref, scores_ht.alt].misbad
        )
        annots.add("misbad")

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
        overwrite=overwrite_temp,
    )

    logger.info("Calculating fitted scores...")
    annot_expr = [
        (scores_ht[variable] * mpc_rel_vars[variable]) for variable in variables
    ]
    interaction_annot_expr = []
    for variable in interactions_dict.keys():
        if len(variable.split(interaction_char)) > 2:
            # NOTE: This assumes that variable is in the format of x:y:z or x*y*z
            # (Assumes the number of variables is maximum 3)
            if len(variable.split(interaction_char)) > 2:
                interaction_annot_expr.append(
                    scores_ht[variable.split(interaction_char)[0]]
                    * scores_ht[variable.split(interaction_char)[1]]
                    * scores_ht[variable.split(interaction_char)[2]]
                    * mpc_rel_vars[variable]
                )
            else:
                interaction_annot_expr.append(
                    scores_ht[variable.split(interaction_char)[0]]
                    * scores_ht[variable.split(interaction_char)[1]]
                    * mpc_rel_vars[variable]
                )

    annot_expr.extend(interaction_annot_expr)
    combined_annot_expr = hl.fold(lambda i, j: i + j, 0, annot_expr)

    logger.info("Computing fitted scores...")
    scores_ht = scores_ht.annotate(fitted_score=intercept + combined_annot_expr)
    scores_ht = scores_ht.checkpoint(
        f"{TEMP_PATH_WITH_FAST_DEL}/fitted_scores_dup{temp_label}.ht"
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
        f"{TEMP_PATH_WITH_FAST_DEL}/fitted_scores_dedup{temp_label}.ht"
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
    use_release: bool = True,
    overwrite_temp: bool = True,
    overwrite_output: bool = True,
    fold: int = None,
    freeze: int = CURRENT_FREEZE,
) -> None:
    """
    Annotate input Table with MPC scores.

    This function can annotate with scores from the MPC release HT (VEP context HT with MPC scores
    calculated from all training data) or it can calculate MPC scores directly using
    a given MPC model (e.g., if saving scores for the whole VEP context HT is not needed).

    To calculate the MPC score, a variant of interest's fitted score from the MPC model
    is combined with the number of common (AF > 0.001) variants in gnomAD with fitted scores
    < the fitted score for the variant of interest.

    For more information on the fitted score and MPC calculation, see the docstring of `run_regressions`.

    :param hl.Table ht: Input Table.
    :param str output_ht_path: Path to write out annotated table (with duplicate variants removed.
    :param bool overwrite_temp: Whether to overwrite intermediate temporary data if it already exists.
        If False, will read existing intermediate temporary data rather than overwriting.
        Default is True.
    # TODO: Update `overwrite_output` documentation
    :param bool overwrite_output: Whether to entirely overwrite final output data if it already exists.
        If False, will read and modify existing output data by adding or modifying columns rather than overwriting entirely.
        If True, will clear existing output data and write new output data.
        The output Tables are the MPC score Tables with duplicated or deduplicated loci.
        Default is True.
    :param int fold: Fold number in training set used to generate MPC model.
        If not None, MPC scores are generated from variants in only training transcripts
        from the specified fold of the overall training set. If None, scores are generated from
        variants in all training transcripts. Default is None.
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
        if fold is not None:
            raise DataException(
                "MPC release is generated only from all training transcripts!"
            )
        if not file_exists(mpc_release.versions[freeze].path):
            raise DataException("MPC release has not yet been generated!")

        logger.info("Annotating MPC scores using MPC release HT and writing out...")
        mpc_ht = mpc_release.versions[freeze].ht()
        ht = ht.annotate(mpc=mpc_ht[ht.locus, ht.alleles].mpc)
        ht.write(output_ht_path, overwrite=overwrite_output)
    else:
        logger.info("Calculating fitted scores on input variants...")
        scores_ht = calculate_fitted_scores(ht=ht.select(), freeze=freeze)

        logger.info("Aggregating fitted scores for gnomAD common variants...")
        fitted_group_path = gnomad_fitted_score_path(
            is_grouped=True, fold=fold, freeze=freeze
        )
        if not file_exists(fitted_group_path):
            aggregate_gnomad_fitted_scores(fold=fold, freeze=freeze)
        gnomad_ht = hl.read_table(fitted_group_path)
        scores = gnomad_ht.aggregate(hl.sorted(hl.agg.collect(gnomad_ht.fitted_score)))
        scores_len = len(scores)

        # Get total number of gnomAD common variants
        gnomad_var_count = hl.read_table(
            gnomad_fitted_score_path(fold=fold, freeze=freeze)
        ).count()

        logger.info("Getting n_less values for input variants...")
        # Annotate HT with sorted array of gnomAD fitted scores
        scores_ht = scores_ht.annotate_globals(gnomad_scores=scores)
        # Checkpoint here to force the gnomAD join to complete
        scores_ht = scores_ht.checkpoint(
            f"{TEMP_PATH_WITH_FAST_DEL}/mpc_temp_gnomad.ht",
            _read_if_exists=not overwrite_temp,
            overwrite=overwrite_temp,
        )
        # TODO: Add temp label here and other places as needed

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
                scores_ht.idx == scores_len,
                gnomad_var_count,
                gnomad_ht[scores_ht.idx].n_less,
            )
        )
        # Checkpoint here to force the binary search to compute
        scores_ht = scores_ht.checkpoint(
            f"{TEMP_PATH_WITH_FAST_DEL}/mpc_temp_binary.ht",
            _read_if_exists=not overwrite_temp,
            overwrite=overwrite_temp,
        )

        logger.info("Calculating MPC scores...")
        scores_ht = scores_ht.annotate(
            mpc=-(hl.log10(scores_ht.n_less / gnomad_var_count))
        )
        scores_ht = scores_ht.checkpoint(
            f"{TEMP_PATH_WITH_FAST_DEL}/mpc.ht",
            overwrite=overwrite_temp,
            _read_if_exists=not overwrite_temp,
        )

        logger.info("Annotating MPC scores to input table and writing out...")
        ht = ht.annotate(mpc=scores_ht[ht.key].mpc)
        ht.write(output_ht_path, overwrite=overwrite_output)


def create_mpc_release_ht(
    overwrite_temp: bool = True,
    overwrite_output: bool = True,
    freeze: int = CURRENT_FREEZE,
) -> None:
    """
    Annotate VEP context Table with MPC scores calculated using training transcripts.

    :param bool overwrite_temp: Whether to overwrite intermediate temporary data if it already exists.
        If False, will read existing intermediate temporary data rather than overwriting.
        Default is True.
    :param bool overwrite_output: Whether to entirely overwrite final output data if it already exists.
        If False, will read and modify existing output data by adding or modifying columns rather than overwriting entirely.
        If True, will clear existing output data and write new output data.
        The output Tables are the MPC score Tables with duplicated or deduplicated loci.
        Default is True.
    :param int freeze: RMC data freeze number. Default is CURRENT_FREEZE.
    :return: None; function writes Table to resource path.
    """
    annotate_mpc(
        ht=context_with_oe.versions[freeze].ht().select(),
        output_ht_path=mpc_release.versions[freeze].path,
        use_release=False,
        overwrite_temp=overwrite_temp,
        overwrite_output=overwrite_output,
        freeze=freeze,
    )
