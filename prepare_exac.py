from gnomad_hail.utils import *
from gnomad_hail.utils.constants import CSQ_ORDER
from regional_missense_constraint.resources.basics import *


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("prepare_exac")
logger.setLevel(logging.INFO)


def filter_to_missense(ht: hl.Table) -> hl.Table:
    """
    Filter ExAC ht to missense variants

    :param Table ht: ExAC ht
    :return: ExAC ht filtered to missense variants
    :rtype: Table
    """
    logger.info(f'ht count before filtration {ht.count()}')
    csqs = hl.literal(CSQ_ORDER)

    ht = ht.explode(ht.info.CSQ)
    ht = ht.annotate(most_severe_consequence=csqs.find(lambda c: ht.info.CSQ.contains(c)))
    logger.info(f'Consequence counts:{ht.aggregate(hl.agg.counter(ht.most_severe_consequence))}')

    logger.info('Filtering to missense variants')
    missense = ['stop_lost', 'initiator_codon_variant', 'start_lost', 'protein_altering_variant', 'missense_variant']
    ht = ht.filter(hl.literal(MISSENSE).contains(ht.most_severe_consequence))
    logger.info(f'ht count after filtration: {ht.count()}')
    return ht

 
def main(args):

    hl.init(log='/prepare_exac.log')

    if args.import_vcf:
        logger.info('Importing ExAC VCF')
        mt = hl.import_vcf(exac_vcf, force_bgz=True)
        ht = mt.rows()
        ht = ht.naive_coalesce(1000).write(exac_ht, overwrite=args.overwrite)

    if args.filter_ht:
        logger.info('Filtering ExAC ht to only missense variants')
        ht = hl.read_table(exac_ht)
        ht = filter_to_missense(ht)
        ht.naive_coalesce(500).write(filtered_exac_ht)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('This script prepares the ExAC sites vcf for RMC testing')

    parser.add_argument('--import_vcf', help='Import ExAC VCF and write to ht', action='store_true')
    parser.add_argument('--filter_ht', help='Filter ExAC ht to missense variants', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite existing data', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@kc')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
