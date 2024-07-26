"""Utilities for loading resource data."""
import logging

import hail as hl
from gnomad.utils.file_utils import check_file_exists_raise_error, file_exists
from gnomad.utils.filtering import filter_to_clinvar_pathogenic
from gnomad.utils.liftover import default_lift_data

from rmc.resources.basics import TEMP_PATH_WITH_FAST_DEL
from rmc.resources.reference_data import (
    autism_de_novo_2022_tsv_path,
    clinvar,
    clinvar_plp_mis_haplo,
    clinvar_plp_mis_triplo,
    dosage_ht,
    dosage_tsv_path,
    haplo_genes_path,
    ndd_de_novo,
    ndd_de_novo_2020_tsv_path,
    triplo_genes_path,
)
from rmc.resources.resource_utils import CURRENT_BUILD, MISSENSE

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("data_loading_utils")
logger.setLevel(logging.INFO)


####################################################################################
## Assessment utils
####################################################################################


def import_dosage(
    overwrite: bool,
    haplo_threshold: float = 0.86,
    triplo_threshold: float = 0.94,
) -> None:
    """
    Import gene dosage sensitivity information.

    Also create HailExpressions with haploinsufficient
    and triplosensitive genes.

    :param overwrite: Whether to overwrite output data.
    :param haplo_threshold: pHaplo score threshold for determining whether a gene is predicted haploinsufficient.
        Default is 0.86 (from Collins et al. paper).
    :param triplo_threshold: pTriplo score threshold for determining whether a gene is predicted triplosensitive.
        Default is 0.94 (from Collins et al. paper).
    :return: None; function writes data to resource paths.
    """
    ht = hl.import_table(dosage_tsv_path, impute=True, force=True)
    ht = ht.transmute(gene=ht["#gene"])
    ht = ht.key_by("gene")
    ht = ht.annotate_globals(
        haplo_cutoff=haplo_threshold,
        triplo_cutoff=triplo_threshold,
    )
    ht = ht.checkpoint(
        dosage_ht.path, _read_if_exists=not overwrite, overwrite=overwrite
    )

    haplo_genes = ht.filter(ht.pHaplo >= haplo_threshold)
    haplo_genes = haplo_genes.aggregate(hl.agg.collect_as_set(haplo_genes.gene))
    hl.experimental.write_expression(
        haplo_genes,
        haplo_genes_path,
        overwrite=overwrite,
    )
    triplo_genes = ht.filter(ht.pTriplo >= triplo_threshold)
    triplo_genes = triplo_genes.aggregate(hl.agg.collect_as_set(triplo_genes.gene))
    hl.experimental.write_expression(
        triplo_genes,
        triplo_genes_path,
        overwrite=overwrite,
    )


def import_clinvar(
    overwrite: bool,
    missense_str: str = MISSENSE,
    build: str = CURRENT_BUILD,
    min_partitions: int = 1000,
) -> None:
    """
    Import ClinVar HT and pathogenic/likely pathogenic missense variants.

    Also filter P/LP missense variants to variants in haploinsufficient
    and triplosensitive genes.

    .. note::
        This function assumes that the ClinVar VCF has been added to the path
        f"{TEMP_PATH_WITH_FAST_DEL}/clinvar.tsv.gz".

    :param bool overwrite: Whether to overwrite output data.
    :param missense_str: String that corresponds to missense variant consequence.
        Default is MISSENSE.
    :param build: Reference genome build. Default is CURRENT_BUILD.
    :param min_partitions: Number of minimum partitions desired for ClinVar HT. Only used if importing ClinVar data from VCF.
        Default is 1000.
    :return: None; writes HTs and HEs to resource paths.
    """
    if (
        not file_exists(clinvar_plp_mis_haplo.path)
        or not file_exists(clinvar_plp_mis_triplo.path)
        or overwrite
    ):
        logger.info("Importing ClinVar VCF and filtering to P/LP missense...")
        clinvar_vcf_path = f"{TEMP_PATH_WITH_FAST_DEL}/clinvar.tsv.gz"
        check_file_exists_raise_error(
            clinvar_vcf_path,
            error_if_not_exists=True,
            error_if_not_exists_msg=(
                f"ClinVar VCF not found at {clinvar_vcf_path}. Please import and try"
                " again!"
            ),
        )
        ht = hl.import_vcf(
            clinvar_vcf_path,
            min_partitions=min_partitions,
            drop_samples=True,
            force_bgz=True,
            reference_genome=build,
            skip_invalid_loci=True,
        )
        ht = ht.checkpoint(clinvar.versions[build].path, overwrite=True)
        ht = filter_to_clinvar_pathogenic(ht)
        ht = ht.annotate(mc=ht.info.MC)
        ht = ht.filter(ht.mc.any(lambda x: x.contains(missense_str)))

        logger.info("Getting gene information from ClinVar HT...")
        ht = ht.annotate(gene=ht.info.GENEINFO.split(":")[0])
        ht = ht.checkpoint(f"{TEMP_PATH_WITH_FAST_DEL}/clinvar.ht", overwrite=True)
        logger.info(
            "Number of variants after filtering to P/LP missense: %i", ht.count()
        )

        logger.info("Filtering to variants in haploinsufficient genes...")
        if not file_exists(haplo_genes_path):
            import_dosage(overwrite)
        hi_genes = hl.experimental.read_expression(haplo_genes_path)
        haplo_ht = ht.filter(hi_genes.contains(ht.gene))
        haplo_ht = haplo_ht.checkpoint(
            clinvar_plp_mis_haplo.path,
            _read_if_exists=not overwrite,
            overwrite=overwrite,
        )
        logger.info(
            "Number of variants after filtering to HI genes: %i", haplo_ht.count()
        )

        logger.info("Filtering to variants in triplosensitive genes...")
        triplo_genes = hl.experimental.read_expression(triplo_genes_path)
        triplo_ht = ht.filter(triplo_genes.contains(ht.gene))
        triplo_ht = triplo_ht.checkpoint(
            clinvar_plp_mis_triplo.path,
            _read_if_exists=not overwrite,
            overwrite=overwrite,
        )
        logger.info(
            "Number of variants after filtering to TS genes: %i",
            triplo_ht.count(),
        )


def import_fu_data(overwrite: bool, liftover: bool = False) -> None:
    """
    Import de novo variants from Fu et al. (2022) paper.

    Function imports variants from TSV into HT, removes malformed rows,
    and optionally lifts data from GRCh38 to GRCh37.

    :param overwrite: Whether to overwrite Table if it exists.
    :param liftover: Whether to lift data from GRCh38 to GRCh37. Default is False.
    :return: None; Function writes Table to temporary path.
    """
    fu_ht = hl.import_table(
        autism_de_novo_2022_tsv_path,
        impute=True,
        # Skip blank lines at the bottom of this TSV
        missing="",
        skip_blank_lines=True,
    )
    # Remove lines from bottom of TSV that are parsed incorrectly upon import
    # These lines contain metadata about the TSV, e.g.:
    # "Supplementary Table 20. The de novo SNV/indel variants used in TADA
    # association analyses from assembled ASD cohorts"
    fu_ht = fu_ht.filter(~hl.is_missing(fu_ht.Role))
    fu_ht = fu_ht.annotate(
        locus=hl.parse_locus(
            hl.format(
                "chr%s:%s",
                fu_ht.Variant.split(":")[0],
                fu_ht.Variant.split(":")[1],
            ),
            reference_genome="GRCh38",
        ),
        alleles=[fu_ht.Variant.split(":")[2], fu_ht.Variant.split(":")[3]],
    )

    if liftover:
        logger.info("Lifting data from b38 to b37...")
        fu_ht = default_lift_data(fu_ht)

    # Rename 'Proband' > 'ASD' and 'Sibling' > 'control'
    fu_ht = fu_ht.transmute(role=hl.if_else(fu_ht.Role == "Proband", "ASD", "control"))
    fu_ht = fu_ht.group_by("locus", "alleles").aggregate(
        role=hl.agg.collect(fu_ht.role),
    )
    fu_ht.write(f"{TEMP_PATH_WITH_FAST_DEL}/fu_dn.ht", overwrite=overwrite)


def import_kaplanis_data(overwrite: bool, liftover: bool = True) -> None:
    """
    Import de novo variants from Kaplanis et al. (2020) paper.

    Function imports variants from TSV into HT, and filters to cases ascertained for
    developmental delay/intellectual disability only.
    Input TSV also contains autistic individuals, but these individuals overlap with
    data from Fu et al. paper and are therefore not retained.

    :param overwrite: Whether to overwrite Table if it exists.
    :param liftover: Whether to liftover data from GRCh37 to GRCh38. Default is True.
    :return: None; Function writes Table to temporary path.
    """
    kap_ht = hl.import_table(ndd_de_novo_2020_tsv_path, impute=True)
    kap_ht = kap_ht.transmute(
        locus=hl.locus(kap_ht.chrom, kap_ht.pos),
        alleles=[kap_ht.ref, kap_ht.alt],
    )
    kap_ht = kap_ht.filter(hl.is_snp(kap_ht.alleles[0], kap_ht.alleles[1]))
    kap_ht = kap_ht.filter(kap_ht.case_control == "DD")
    kap_ht = kap_ht.group_by("locus", "alleles").aggregate(
        case_control=hl.agg.collect(kap_ht.case_control)
    )
    if liftover:
        logger.info("Lifting data from b37 to b38...")
        kap_ht = default_lift_data(kap_ht)
    kap_ht.write(f"{TEMP_PATH_WITH_FAST_DEL}/kaplanis_dn.ht", overwrite=overwrite)


def import_de_novo_variants(
    overwrite: bool, n_partitions: int = 5000, liftover_b38: bool = True
) -> None:
    """
    Import de novo missense variants.

    .. note::
        These files currently only exist for build GRCh37.

    :param bool overwrite: Whether to overwrite de novo Table.
    :paran n_partitions: Number of partitions for input Tables.
        Used to repartition Tables on read.
        Will also help determine number of partitions in final Table.
    :param liftover_b38: Whether to lift data from GRCh37 to GRCh38. Default is True.
    :return: None; writes HT to resource path.
    """
    fu_ht_path = f"{TEMP_PATH_WITH_FAST_DEL}/fu_dn.ht"
    kaplanis_ht_path = f"{TEMP_PATH_WITH_FAST_DEL}/kaplanis_dn.ht"
    if not file_exists(fu_ht_path) or overwrite:
        # NOTE: Fu data is in GRCh38
        import_fu_data(overwrite=overwrite, liftover=not liftover_b38)
    if not file_exists(kaplanis_ht_path) or overwrite:
        # NOTE: Kaplanis data is in GRCh37
        import_kaplanis_data(overwrite=overwrite, liftover=liftover_b38)

    fu_ht = hl.read_table(fu_ht_path, _n_partitions=n_partitions)
    kap_ht = hl.read_table(kaplanis_ht_path, _n_partitions=n_partitions)
    ht = kap_ht.join(fu_ht, how="outer")

    # Union sample types (DD, ASD, control)
    ht = ht.transmute(
        sample_set=hl.case()
        .when(
            hl.is_defined(ht.case_control) & hl.is_defined(ht.role),
            ht.case_control.extend(ht.role),
        )
        .when(hl.is_defined(ht.case_control), ht.case_control)
        .when(hl.is_defined(ht.role), ht.role)
        .or_missing()
    )
    ht.write(ndd_de_novo.path, overwrite=overwrite)
