from gnomad_hail.utils.generic import get_reference_genome
from gnomad_hail.utils.slack import *
from regional_missense_constraint.resources.basics import *


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("regional_missense_constraint")
logger.setLevel(logging.INFO)


def main(args):

    hl.init(log='/RMC.log')

    if args.pre_process_data:
        logger.info('Preprocessing reference fasta and gencode files')
        process_context_ht(args.build, args.trimers)

        logger.info('Preprocessing gencode gtf information')
        gencode_ht = process_gencode_ht(args.build)

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
    if args.exac:
        exome_ht = prepare_ht(hl.read_table(filtered_exac_ht), args.trimers)
    else:   
        exome_ht = prepare_ht(hl.read_table(filtered_exomes_ht_path), args.trimers)  
    
    logger.info('Inferring build of exome ht')
    rg = get_reference_genome(exome_ht.locus).name 

    logger.info('Reading in context ht and gencode ht') 
    context_ht = hl.read_table(get_processed_context_ht_path(args.build))
    gencode_ht = get_processed_gencode_ht(args.build)
 
    if args.test:
        contigs = rg.contigs[21]
        logger.info('Filtering to chr22 for testing')
        context_ht = hl.filter_intervals(context_ht, [hl.parse_locus_interval(rg.contigs[21], reference_genome=rg)])
        exome_ht = hl.filter_intervals(exome_ht, [hl.parse_locus_interval(rg.contigs[21], reference_genome=rg)])
    else:
        logger.info('Filtering to autosomes and PAR regions')
        context_ht = filter_to_autosome_and_par(context_ht)
        exome_ht = filter_to_autosome_and_par(exome_ht)

    # NOTE: average divergence score from Kaitlin
    # https://github.com/ksamocha/regional_constraint/blob/master/per_base_regional.py#L1547
    logger.info('Annotating transcripts in context ht with divergence scores')
    context_ht = context_ht.annotate(div=(hl.cond(
                                                context_ht.div_scores.contains(context_ht.vep.transcript_consequences.transcript_id),
                                                context_ht.div_scores[context_ht.vep.transcript_consequences.transcript_id],
                                                0.0564635)))
    
    #logger.info('Processing full context')


if __name__ == '__main__':
    parser = argparse.ArgumentParser('This script searches for regional missense constraint in gnomAD')

    parser.add_argument('--exac', help='Use ExAC ht (not gnomAD ht)', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite existing data', action='store_true')
    parser.add_argument('--pre_process_data', help='Pre-process data', action='store_true')
    parser.add_argument('--test', help='Filter to chr22 (for code testing purposes)', action='store_true')
    parser.add_argument('--trimers', help='Use trimers instead of heptamers', action='store_false')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@kc')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
    main(args)
