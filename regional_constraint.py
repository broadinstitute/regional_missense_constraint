from gnomad_hail.utils.generic import get_reference_genome
from gnomad_hail.utils.slack import *
from regional_missense_constraint.resources.basics import *
from regional_missense_constraint.utils.generic import *


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("regional_missense_constraint")
logger.setLevel(logging.INFO)


def calculate_expected(context_ht: hl.Table, exac: bool) -> hl.Table:
    """
    Annotates context ht with expected variants count (adjusted by mutation rate, divergence score, and region type

    :param Table context_ht: Context ht
    :param bool exac: Whether the data is ExAC data
    :return: Context ht with expected variant annotations
    :rtype: Table
    """
    logger.info('Adjusting mutation rate with divergence scores')
    # p_mut2 = p_mut1*(1 + (0.31898*div_score))
    context_ht = context_ht.annotate(mu=(context_ht.mu_snp * (1 + (0.31898*context_ht.div))))

    logger.info('Processing context ht to calculate number of expected variants')
    # NOTE: count variants from Konrad's constraint code
    # https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/generic.py
    '''exp_counts = count_variants(context_ht).variant_count # i think this pickle is actually used for mu calculations
    outfile = exac_exp_var_pickle if exac else exp_var_pickle
    with hl.hadoop_open(outfile, 'wb') as o:
        pickle.dump(exp_counts, o)'''
    logger.info('Counting variants in context ht')
    exp_ht = count_variants(context_ht)
    exp_ht = exp_ht.transmute(possible_variants=exp_ht.variant_count) # rename variant_counts (from count_variants function) to possible variants
    if exac:
        exp_ht = exp_ht.annotate(expected_variants=hl.case()
                                                        .when(exp_ht.region_type == 'x_nonpar', (0.3715167 + 7796945 * exp_ht.mu))
                                                        .when(exp_ht.region_type == 'y', (0.05330181 + 2457366 * exp_ht.mu))
                                                        .default(0.4190964 + 11330208 * exp_ht.mu))

    # TODO: get expected variant calculation for gnomAD
    return exp_ht


def calculate_observed(exome_ht: hl.Table, exac: bool) -> hl.Table:
    """
    Something about observed variants

    :param Table exome_ht: Input exome ht with observed variants
    :param bool exac: Whether the data is ExAC data
    :return: Exome ht something with observed variants
    :rtype: Table
    """
    if exac:
        # keep criteria from Kaitlin (manuscript, also in her code): adjusted AC <= 123 and VQSLOD >= -2.632
        obs_ht = exome_ht.filter((exome_ht.info.AC_Adj <= 123) & (exome_ht.info.VQSLOD >= -2.632))
        obs_ht = count_variants(obs_ht)
        obs_ht = obs_ht.transmute(observed_variants=obs_ht.variant_count)

    return obs_ht

    # TODO: get observed counts for gnomAD


def main(args):

    hl.init(log='/RMC.log')
    exac = args.exac

    if args.pre_process_data:
        logger.info('Preprocessing reference fasta and gencode files')
        process_context_ht('GRCh37', args.trimers)

        logger.info('Preprocessing gencode gtf information')
        gencode_ht = process_gencode_ht('GRCh37')

        logger.info('Filtering gnomAD exomes ht to missense variants only')
        exome_ht = hl.read_table(processed_exomes_ht_path)
        exome_ht = filter_to_missense(exome_ht)
        exome_ht.write(filtered_exomes_ht_path, overwrite=args.overwrite)

        logger.info('Done preprocessing files')

        # NOTE: from April 2019
        # total rows in context_ht:  4014886371 (full_context_ht_path)
        # total rows in context_ht after filtering to missenses in protein coding, canonical transcripts: 92683899 # this is processed context ht
        # total rows scored with average divergence score (=not in div_scores): 49443864
    
    logger.info('Reading in exome ht') 
    if exac:
        exome_ht = hl.read_table(filt_exac_cov_ht)
    else:   
        exome_ht = prepare_ht(hl.read_table(filtered_exomes_ht_path), args.trimers)  
    
    logger.info('Inferring build of exome ht')
    rg = get_reference_genome(exome_ht.locus) 

    logger.info('Reading in context ht and gencode ht') 
    context_ht = hl.read_table(get_processed_context_ht_path(rg.name))
    gencode_ht = get_processed_gencode_ht_path(rg.name)

    if args.test:
        contigs = rg.contigs[21]
        logger.info('Filtering to chr22 for testing')
        context_ht = hl.filter_intervals(context_ht, [hl.parse_locus_interval(rg.contigs[21], reference_genome=rg)])
        exome_ht = hl.filter_intervals(exome_ht, [hl.parse_locus_interval(rg.contigs[21], reference_genome=rg)])
    else:
        logger.info('Filtering to autosomes and PAR regions')
        context_ht = filter_to_autosome_and_par(context_ht)
        exome_ht = filter_to_autosome_and_par(exome_ht)

    if args.calc_exp:
        exp_ht = calculate_expected(context_ht, exac)

    if args.calc_obs:
        pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser('This script searches for regional missense constraint in gnomAD')

    parser.add_argument('--pre_process_data', help='Pre-process data', action='store_true')
    parser.add_argument('--calc_exp', help='Calculate expected variant counts', action='store_true')  
    parser.add_argument('--calc_obs', help='Calculated observed variant counts', action='store_true')
 
    parser.add_argument('--trimers', help='Use trimers instead of heptamers', action='store_false')  
    parser.add_argument('--exac', help='Use ExAC ht (not gnomAD ht)', action='store_true')
    parser.add_argument('--test', help='Filter to chr22 (for code testing purposes)', action='store_true')
    parser.add_argument('--pre_process_data', help='Pre-process data', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@kc')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
