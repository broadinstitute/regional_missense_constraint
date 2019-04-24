import hail as hl
from regional_missense_constraint.resources.basics import *


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("regional_missense_constraint")
logger.setLevel(logging.INFO)


def main(args):

    if args.pre_process_data:
        logger.info('Preprocessing reference fasta and gnomAD files')
        pre_process_data(args.build)

    #full_genome_ht = prepare_ht(get_processed_genomes_ht_path(), args.trimers) 
    #full_exome_ht = prepare_ht(get_processed_exomes_ht_path(), args.trimers)
    #full_genome_ht = prepare_ht(f'{FLAGSHIP_LOF}/model/genomes_processed.ht', args.trimers)
    full_exome_ht = prepare_ht(hl.read_table(f'{FLAGSHIP_LOF}/model/exomes_processed.ht'), args.trimers)

    if args.test:
        logger.info('Filtering to chr21 for testing')
        if args.build == 'GRCh37':
            contig = '21'
        else:
            contig = 'chr21'
        exome_ht = full_exome_ht.filter(full_exome_ht.locus.contig == contig)
    else:
        logger.info('Filtering to autosomes and PAR regions')
        exome_ht = full_exome_ht.filter(full_exome_ht.locus.in_autosome_or_par())

    logger.info('Import codon translation table and amino acid names and annotate as globals')
    exome_ht = exome_ht.annotate_globals(codon_translation=get_codon_lookup(), acid_names=get_acid_names())
    exome_ht.describe()
            

if __name__ == '__main__':
    parser = argparse.ArgumentParser('This script searches for regional missense constraint in gnomAD')

    pparser.add_argument('--build', help='Reference genome of dataset (37 or 38)', type=int, choices=[37, 38], default=37)
    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--pre_process_data', help='Pre-process data', action='store_true')
    parser.add_argument('--test', help='Filter to chr21 (for code testing purposes)', action='store_true')
    parser.add_argument('--trimers', help='Use trimers instead of heptamers', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@kc')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel.split(','), main, args)
    else:
        main(args)
