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
    logger.info('Deduplicating keys')
    ht = ht.distinct()
    logger.info(f'ht count after dedup: {ht.count()}')
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
        ht.naive_coalesce(500).write(filtered_exac_ht, overwrite=args.overwrite)

    #chroms = ['X', 'Y']
    #for i in range(1, 23):
    #    chroms.append[i]
    chroms = ['22']

    # NOTE: only calculated for chr22
    if args.import_cov:
        for c in chroms:
            tsv = f'{exac_tsv_path}/Panel.chr{c}.coverage.txt.gz'
            out = f'{exac_cov_path}/{c}_coverage.ht'
            ht = hl.import_table(tsv, min_partitions=100, impute=True, force_bgz=True)
            ht = ht.transmute(locus=hl.parse_locus(hl.format('%s:%s', ht['#chrom'], ht.pos)))
            ht = ht.rename({'1': 'over_1', '5': 'over_5', '10': 'over_10', '15': 'over_15', '20': 'over_20', '25': 'over_25', 
                            '30': 'over_30', '50': 'over_50', '100': 'over_100'})
            ht = ht.key_by('locus')
            ht.write(out, overwrite=args.overwrite)

    if args.join_cov:
        ht = hl.read_table(filtered_exac_ht)
        for c in chroms:
            cov_path = f'{exac_cov_path}/{c}_coverage.ht'
            cov_ht = hl.read_table(cov_path)
            ht = ht.annotate(coverage=cov_ht[ht.locus])
        ht.write(filt_exac_cov_ht, overwrite=args.overwrite)



if __name__ == '__main__':
    parser = argparse.ArgumentParser('This script prepares the ExAC sites vcf for RMC testing')

    parser.add_argument('--import_vcf', help='Import ExAC VCF and write to ht', action='store_true')
    parser.add_argument('--filter_ht', help='Filter ExAC ht to missense variants', action='store_true')
    parser.add_argument('--import_cov', help='Import coverage files', action='store_true')
    parser.add_argument('--join_cov', help='Annotate ExAC ht with coverage', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite existing data', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@kc')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
