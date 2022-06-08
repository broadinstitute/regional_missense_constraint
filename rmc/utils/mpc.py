import logging
import pandas as pd
from patsy import dmatrices
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
    context_with_oe,
    context_with_oe_dedup,
    gnomad_fitted_score,
    gnomad_fitted_score_group,
    grantham,
    grantham_txt_path,
    joint_clinvar_gnomad,
    misbad,
    mpc_model_pkl_path,
    mpc_release,
    mpc_release_dedup,
    temp_path,
)
from rmc.resources.grch37.reference_data import cadd, clinvar_path_mis
from rmc.resources.resource_utils import MISSENSE
from rmc.utils.generic import get_aa_map, get_constraint_transcripts, process_vep
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


def get_annotations_from_context_ht_vep(
    ht: hl.Table, annotations: List[str] = ["polyphen", "sift"]
) -> hl.Table:
    """
    Get desired annotations from context Table's VEP annotation and annotate/filter to informative amino acids.

    Function will get 'codons', 'most_severe_consequence', and 'transcript' annotations by default.
    Additional annotations should be specified using `annotations`.

    :param hl.Table ht: Input VEP context Table.
    :param List[str] annotations: Additional annotations to extract from context HT's VEP field.
    :return: VEP context Table with rows filtered to remove loci that are non-coding and retain only those with informative amino acids
        and columns filtered to keep only specified annotations.
    """
    annot_expr = {
        "codons": ht.transcript_consequences.codons,
        "most_severe_consequence": ht.transcript_consequences.most_severe_consequence,
        "transcript": ht.transcript_consequences.transcript_id,
    }
    if "polyphen" in annotations:
        annot_expr["polyphen"] = hl.struct(
            prediction=ht.transcript_consequences.polyphen_prediction,
            score=ht.transcript_consequences.polyphen_score,
        )
    if "sift" in annotations:
        annot_expr["sift"] = hl.struct(
            prediction=ht.transcript_consequences.sift_prediction,
            score=ht.transcript_consequences.sift_score,
        )
    ht = ht.annotate(**annot_expr)
    annotations.extend(["codons", "most_severe_consequence", "transcript"])
    ht = ht.select(*annotations)
    return annotate_and_filter_codons(ht)


def create_context_with_oe(
    missense_str: str = MISSENSE,
    temp_path_with_del: str = "gs://gnomad-tmp/mpc",
    n_partitions: int = 30000,
    overwrite: bool = False,
) -> None:
    """
    Filter VEP context Table to missense variants in canonical transcripts, and add missense observed/expected.

    Function writes two Tables:
        - `context_with_oe`: Table of all missense variants in canonical transcripts + all annotations (includes duplicate loci).
            Annotations are: transcript, most severe consequence, codons, reference and alternate amino acids,
            Polyphen-2, and SIFT.
        - `context_with_oe_dedup`: Deduplicated version of `context_with_oe` that only contains missense o/e and transcript annotations.

    :param str missense_str: String representing missense variant consequence. Default is MISSENSE.
    :param str temp_path_with_del: Path to bucket to store temporary data with automatic deletion policy.
        Default is 'gs://gnomad-tmp/mpc'.
        TODO: Update this to `temp_path` (and set automatic deletion policy.)
    :param int n_partitions: Number of desired partitions for the VEP context Table.
        Repartition VEP context Table to this number on read.
        Default is 30000.
    :param bool overwrite: Whether to overwrite temporary data if it already exists.
        If False, will read existing temporary data rather than overwriting.
        Default is False.
    :return: None; function writes Table to resource path.
    """
    logger.info("Importing set of transcripts to keep...")
    transcripts = get_constraint_transcripts(outlier=False)

    logger.info(
        "Reading in VEP context HT and filtering to missense variants in canonical transcripts..."
    )
    # Using vep_context.path to read in table with fewer partitions
    # VEP context resource has 62164 partitions
    ht = (
        hl.read_table(vep_context.path, _n_partitions=n_partitions)
        .select_globals()
        .select("vep", "was_split")
    )
    ht = process_vep(ht, filter_csq=True, csq=missense_str)
    ht = ht.filter(transcripts.contains(ht.transcript_consequences.transcript_id))
    ht = get_annotations_from_context_ht_vep(ht)
    # Save context Table to temporary path with specified deletion policy because this is a very large file
    # and relevant information will be saved at `context_with_oe` (written below)
    # Checkpointing here to force this expensive computation to compute
    # before joining with RMC tables in `get_oe_annotation`
    # This computation is expensive, so overwrite only if specified (otherwise, read existing file)
    ht = ht.checkpoint(
        f"{temp_path_with_del}/vep_context_mis_only_annot.ht",
        _read_if_exists=not overwrite,
        overwrite=overwrite,
    )

    logger.info(
        "Adding regional missense constraint missense o/e annotation and writing to resource path..."
    )
    ht = get_oe_annotation(ht)
    ht = ht.checkpoint(
        context_with_oe.path, _read_if_exists=not overwrite, overwrite=overwrite
    )

    logger.info(
        "Creating dedup context with oe (with oe and transcript annotations only)..."
    )
    ht = ht.select("transcript", "oe")
    ht = ht.collect_by_key()
    ht = ht.annotate(oe=hl.min(ht.values.oe))
    ht = ht.annotate(transcript=ht.values.find(lambda x: x.oe == ht.oe).transcript)
    ht.write(
        context_with_oe_dedup.path, _read_if_exists=not overwrite, overwrite=overwrite
    )


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
    if not file_exists(context_with_oe.path):
        create_context_with_oe(overwrite=True)
    context_ht = context_with_oe.ht()

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

    logger.info("Filtering to rows with defined annotations and writing out...")
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
        n_less(v) = number of common (AF > 0.01) ExAC variants with fitted_score < fitted_score(v)

    Note that higher MPC scores predict increased missense deleteriousness, and
    smaller n_less values and fitted scores will lead to higher MPC scores.

    :param List[str] variables: Variables to include in all regressions (single, joint).
        Default is ["oe", "misbad", "polyphen"].
    :param List[str] additional_variables: Additional variables to include in single variable regressions only.
        Default is ["blosum", "grantham"].
    :return: None; function writes Table with gnomAD fitted scores
        and model coefficients as pickle to resource paths.
    """
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
    single_var_X = []
    all_var = variables + additional_variables
    for var in all_var:
        logger.info("Running %s regression...", var)
        # Create design matrices
        formula = f"pop_v_path ~ {var}"
        X, model = _run_glm(formula)
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
    logger.info("Lowest model AIC: %f", min_aic)
    '''if all_model_aic.count(min_aic) > 1:
        logger.warning(
            """
            There is a tie for minimum AIC.
            This function will use the first model it finds by default
            (single variable -> additive interactions -> all interactions -> specific interactions)!
            """
        )
    if min_aic == min_single_aic:
        logger.info(
            "Single variable regression using %s had the lowest AIC", min_single_aic_var
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
        X = spec_X'''

    model = spec_model
    X = spec_X
    logger.info("Saving model as pickle...")
    with hl.hadoop_open(mpc_model_pkl_path, "wb") as p:
        pickle.dump(model, p, protocol=pickle.HIGHEST_PROTOCOL)

    logger.info(
        "Annotating gnomAD variants with fitted score and writing to gnomad_fitted_score path..."
    )
    ht = ht.filter(ht.pop_v_path == 1)
    ht = calculate_fitted_scores(ht)
    ht.write(gnomad_fitted_score.path, overwrite=True)


def calculate_fitted_scores(
    ht: hl.Table, interaction_char: str = ":", intercept_str: str = "Intercept"
) -> hl.Table:
    """
    Use model and relationship determined in `run_regressions` to calculate fitted scores for input Table.

    .. note::
        - Assumes model has been written as pickle to `mpc_model_pkl_path`.

    :param hl.Table ht: Input Table with variants to be annotated.
    :param str interaction_char: Character representing interactions in MPC model. Must be one of "*", ":".
        Default is ":".
    :param str intercept_str: Name of intercept variable in MPC model pickle. Default is "Intercept".
    :return: Table annotated with fitted score.
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
    variables = mpc_rel_vars.keys()
    annotations = []
    if "misbad" in variables:
        logger.info("Getting missense badness annotation...")
        mb_ht = misbad.ht()
        ht = ht.annotate(misbad=mb_ht[ht.ref, ht.alt].misbad)
        annotations.append("misbad")

    if "oe" in variables:
        annotations.append("oe")

    if "polyphen" in variables:
        ht = ht.transmute(polyphen=ht.polyphen.score)
        annotations.append("polyphen")

    logger.info("Filtering to defined annotations and checkpointing...")
    filter_expr = True
    for field in annotations:
        filter_expr &= hl.is_defined(ht[field])
    ht = ht.filter(filter_expr)
    # Checkpoint here to force missense badness join and filter to complete
    ht = ht.checkpoint(f"{temp_path}/fitted_score_temp_join.ht", overwrite=True)

    logger.info("Annotating fitted scores...")
    variable_dict = {
        variable: mpc_rel_vars[variable]
        for variable in variables
        if not interaction_char in variable
    }
    interactions_dict = {
        variable: mpc_rel_vars[variable]
        for variable in variables
        if interaction_char in variable
    }
    for variable in interactions_dict.keys():
        if len(variable.split(interaction_char)) > 2:
            raise DataException(
                "Code currently doesn't handle interactions between more than two variables!"
            )
    annot_expr = [
        (ht[variable] * mpc_rel_vars[variable]) for variable in variable_dict.keys()
    ]
    # NOTE: This assumes that variable is in the format of x:y or x*y
    # and won't handle cases where variable is x:y:z or x*y*z correctly
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

    logger.info("Computing fitted score and returning...")
    return ht.annotate(fitted_score=intercept + combined_annot_expr)


def aggregate_gnomad_fitted_scores(n_less_eq0_float: float = 0.83) -> None:
    """
    Aggregate gnomAD fitted scores to count number of variants with a score less than a given score.

    :param float n_less_eq0_float: Set `n_less` annotation to this float if `n_less` is 0.
        This avoids errors in the `hl.log10` call and ensures that MPC for variants with a fitted score
        more severe than any common gnomAD variant score (`n_less` = 0) is more severe (by a controlled amount)
        compared to MPC for variants with a fitted score more severe than one common gnomAD variant (`n_less` = 1).
    :return: None; function writes Table to resource path.
    """
    logger.info("Aggregating gnomAD fitted scores...")
    gnomad_ht = gnomad_fitted_score.ht()
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
    gnomad_ht.write(gnomad_fitted_score_group.path, overwrite=True)


def create_mpc_release_ht(
    overwrite: bool = True,
) -> None:
    """
    Annotate VEP context Table with MPC component variables and calculate MPC using relationship defined in `mpc_rel_vars`.

    Relationship in `mpc_rel_vars` is the formula used to calculate a variant's fitted score.
    A variant of interest's fitted score is combined with the number of common (AF > 0.01)
    variants in gnomAD with fitted scores < the fitted score for the variant of interest
    to determine a variant's MPC score.

    For more information on the fitted score and MPC calculation, see the docstring of `run_regressions`.

    :param bool overwrite: Whether to overwrite output table if it already exists. Default is True.
    :return: None; function writes Table to resource path.
    """
    logger.info("Calculating fitted scores...")
    ht = context_with_oe.ht()
    ht = calculate_fitted_scores(ht)

    logger.info("Aggregating gnomAD fitted scores...")
    if not file_exists(gnomad_fitted_score_group.path):
        aggregate_gnomad_fitted_scores()
    gnomad_ht = gnomad_fitted_score_group.ht()
    scores = gnomad_ht.aggregate(hl.sorted(hl.agg.collect(gnomad_ht.fitted_score)))
    scores_len = len(scores)

    # Get total number of gnomAD common variants
    common_gnomad_ht = joint_clinvar_gnomad.ht()
    common_gnomad_ht = common_gnomad_ht.filter(common_gnomad_ht.pop_v_path == 1)
    gnomad_var_count = common_gnomad_ht.count()

    logger.info("Getting n_less annotation...")
    # Annotate HT with sorted array of gnomAD fitted scores
    ht = ht.annotate_globals(gnomad_scores=scores)
    # Checkpoint here to force the gnomAD join to complete
    ht = ht.checkpoint(f"{temp_path}/mpc_temp_gnomad.ht", overwrite=True)

    # Search all gnomAD scores to find first score that is
    # greater than or equal to score to be annotated
    # `binary_search` will return the index of the first gnomAD fitted score that
    # is >= the score of interest
    # e.g., if the score of interest is 0.45, and gnomAD fitted scores are
    # [0.3, 0.4, 0.5], then `binary_search` will return an index of 2
    # the `n_less` of 0.5 will contain the counts of variants with gnomAD scores of
    # 0.3 and 0.4 due to the non-inclusive nature of scans
    # (n_less[0.5] = n_var[0.3] + n_var[0.4])
    ht = ht.annotate(idx=hl.binary_search(ht.gnomad_scores, ht.fitted_score))
    ht = ht.annotate(
        n_less=hl.if_else(
            # Make n_less equal to total gnomAD common variant count if
            # index is equal to the length of the gnomAD scores array
            ht.idx == scores_len,
            gnomad_var_count,
            gnomad_ht[ht.idx].n_less,
        )
    )
    # Checkpoint here to force the binary search to compute
    ht = ht.checkpoint(f"{temp_path}/mpc_temp_binary.ht", overwrite=True)
    ht = ht.annotate(mpc=-(hl.log10(ht.n_less / gnomad_var_count)))
    ht.write(mpc_release.path, overwrite=overwrite)


def annotate_mpc(
    ht: hl.Table,
    output_path: str,
    overwrite: bool = True,
) -> None:
    """
    Annotate Table with MPC score using MPC release Table (`mpc_release`).

    .. note::
        - Assume input Table is keyed by locus and alleles.

    :param hl.Table ht: Input Table to be annotated.
    :param str output_path: Where to write Table after adding MPC annotations.
    :param bool overwrite: Whether to overwrite specified output path if it already exists.
        Default is True.
    :return: None; function writes Table to specified output path.
    """
    assert (
        len(ht.key) == 2 and "locus" in ht.key and "alleles" in ht.key
    ), "HT key must be keyed by 'locus' and 'alleles'!"
    # NOTE: hl.tlocus() will use the reference genome build from the hail initialization
    # if the user doesn't specify a reference build
    # This means that for v2, where hail is initialized with GRCh37, this will check for
    # a build 37 locus
    assert ht.key.locus.dtype == hl.tlocus() and ht.key.alleles.dtype == hl.tarray(
        hl.tstr
    ), "'locus' must be a LocusExpression, and 'alleles' must be an array of strings!"

    if not file_exists(mpc_release_dedup.path):
        mpc_ht = mpc_release.ht()
        mpc_ht = mpc_ht.select("transcript", "mpc")
        mpc_ht = mpc_ht.collect_by_key()
        mpc_ht.write(mpc_release_dedup.path, overwrite=True)
    mpc_ht = mpc_release_dedup.ht()
    ht = ht.annotate(mpc_list=mpc_ht[ht.key].values)
    ht.write(output_path, overwrite=overwrite)
