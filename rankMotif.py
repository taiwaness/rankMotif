#!/usr/bin/env python
#
# rankMotif
#
# Given a motif list, rankMotif will output the top-k motifs
# without redundancy.
#
# Author: Jian-Long Huang <jianlong@ntu.edu.tw>

__version__ = '1.5'

import os
import sys
import logging
import argparse
from rankmotif.basic import Pattern, PatternSet, merge_patterns
from rankmotif.scoring import PatternScoring
from rankmotif.ranking import Cluster, pfm
from rankmotif.seqio import parse_fasta


def main():
    parser = argparse.ArgumentParser(prog='rankMotif',
                                     description='Given a motif list, '
                                     'rankMotif will output the top-k motifs'
                                     'without redundancy.')
    parser.add_argument('-pset', required=True, metavar='<file>',
                        help='positive set of Chip-Chip sequences in FASTA format')
    parser.add_argument('-nset', required=True, metavar='<file>',
                        help='negative set of Chip-Chip sequences in FASTA format')
    parser.add_argument('-plist', required=True, metavar='<file>',
                        help='pattern list file')
    parser.add_argument('-seqtype', required=True, choices=['dna', 'rna'],
                        help='sequence type of the patterns')
    parser.add_argument('-out', required=True, metavar='<file>',
                        help='output directory')
    # parser.add_argument('-cpu', type=int, default=1,
    #                     help='Number of CPUs to perform the analysis (default: 1)')
    parser.add_argument('-oc', metavar='<file>',
                        help='support of nucleosome occupancy scores')
    parser.add_argument('-cs', metavar='<file>',
                        help='support of conservation scores')
    parser.add_argument('-sp', type=int, default=1, metavar='<int>',
                        help='weight of position scoring (default: 1)')
    parser.add_argument('-sn', type=int, default=1, metavar='<int>',
                        help='weight of nucleosome occupancy scoring. '
                        'This option is valid only when -oc is specified. (default: 1)')
    parser.add_argument('-sc', type=int, default=1, metavar='<int>',
                        help='weight of conservation scoring (default: 1)')
    parser.add_argument('-nc', type=int, default=5, metavar='<int>',
                        help='maximum number of clusters in the output (default: 5)')
    parser.add_argument('-np', type=int, default=5, metavar='<int>',
                        help='maximum number of patterns per cluster (default: 5)')
    parser.add_argument('-ws', type=int, metavar='<int>',
                        help='maximum window size of simpfm (default: auto)')
    parser.add_argument('-gc', type=float, metavar='<float>',
                        help='GC contents (default: auto)')
    parser.add_argument('-seqmask', choices=['yes', 'no'], default='no',
                        help='applying sequence mask (default: no)')
    parser.add_argument('-log', metavar='<file>',
                        help='log file (default: stdout)')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {0}'.format(__version__))
    args = parser.parse_args()

    log_config = {
        'level': logging.INFO,
        'format': '%(asctime)s [%(levelname)s] %(name)s: %(message)s',
        'datefmt': '%Y-%m-%d %H:%M:%S',
    }
    if args.log:
        log_config.update({'filename': args.log})
    else:
        log_config.update({'stream': sys.stdout})

    logging.basicConfig(**log_config)

    logger = logging.getLogger('main')

    if not os.path.exists(args.out):
        os.makedirs(args.out)

    if args.seqtype == 'dna':
        reverse_complement = True
    else:
        reverse_complement = False

    if args.seqmask == 'yes':
        seqmask = True
    else:
        seqmask = False

    pattern_set = PatternSet(reverse_complement)

    logger.info('building match tables of patterns')
    with open(args.plist, 'r') as fi:
        for line in fi:
            pattern = Pattern(line.strip())
            pattern.build_matchtable_pset(parse_fasta(args.pset), reverse_complement)
            pattern.build_matchtable_nset(parse_fasta(args.nset), reverse_complement)
            pattern_set.add(pattern)

    pattern_scoring = PatternScoring(sp_weight=args.sp, sn_weight=args.sn,
                                     sc_weight=args.sc)
    pattern_scoring.build(pattern_set, seqmask=seqmask, nuclocc=args.oc, consv=args.cs)

    cluster = Cluster(args.nc, 0.8, args.np, args.ws, reverse_complement)
    cluster.run(pattern_scoring, gc=args.gc, pset=args.pset)

    cluster_pfm = {}
    logger.info('merging patterns and calculating PFMs')

    with open(os.path.join(args.out, 'clustered_patterns.txt'), 'w') as fo_clu, \
            open(os.path.join(args.out, 'merged_patterns.txt'), 'w') as fo_mer, \
            open(os.path.join(args.out, 'match_sequences.txt'), 'w') as fo_mth:
        fo_clu.write('\t'.join(['cluster_no', 'pattern', 'pset_support\n']))
        fo_mer.write('\t'.join(['cluster_no', 'strand', 'pattern', 'pset_support\n']))
        fo_mth.write('\t'.join(['cluster_no', 'gene_name', 'start', 'sequence', 'strand\n']))
        # fo_mth.write('cluster: sequence\n')

        for i, j in cluster.results.iteritems():
            for p in j:
                pset_support = float(p.matchtable_pset.n_hitseqs) / p.matchtable_pset.n_seqs
                fo_clu.write('\t'.join([
                    str(i),
                    p.sequence.upper(),
                    str(round(pset_support, 2)),
                ]))

                fo_clu.write('\n')
                fo_clu.flush()

            merged = merge_patterns(j, reverse_complement)
            match_sequences = merged.extract_match_info(args.pset)

            for p in merged.patterns:
                strand = merged._strands.get(p)
                if strand == 1:
                    strand = '+'
                else:
                    strand = '-'
                pset_support = float(p.matchtable_pset.n_hitseqs) / p.matchtable_pset.n_seqs
                fo_mer.write('\t'.join([str(i), strand, p.sequence.upper(), str(round(pset_support, 2))]))
                fo_mer.write('\n')
                fo_mer.flush()

            for gene_name, strand, start, seq in match_sequences:
                if strand == 1:
                    strand = '+'
                else:
                    strand = '-'
                fo_mth.write('\t'.join([str(i), gene_name, str(start), seq.upper(), strand]))
                fo_mth.write('\n')
                fo_mth.flush()

            cluster_pfm.update({i: pfm([x[3] for x in match_sequences])})

    for i, j in cluster_pfm.iteritems():
        with open(os.path.join(args.out, 'cluster_{0}.pfm.txt'.format(i)), 'w') as fo:
            fo.write('cluster_{0}\t{1}\n'.format(i, str(len(j.get('a')))))
            for base in ['a', 't', 'c', 'g']:
                fo.write('\t'.join([str(i) for i in j.get(base)]))
                fo.write('\n')
                fo.flush()

    logger.info('Job has finished.')


if __name__ == '__main__':
    if sys.hexversion > 0x03000000:
        sys.exit('Unsupported python version: %s' % sys.version)
    main()
