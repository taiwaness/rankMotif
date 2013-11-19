#!/usr/bin/env python
#
# seetfbs - Discovering transcription factor binding sites from noisy data
#
# Created: 2013.11.19
# Version: 1.0

import os
import sys
import logging
import argparse
from seetfbs.basic import Pattern, PatternSet, merge_patterns
from seetfbs.scoring import PatternScoring
from seetfbs.ranking import Cluster, pfm
from seetfbs.seqio import parse_fasta_noheader


def main():
    parser = argparse.ArgumentParser(
        description='Discovering transcription factor binding sites from noisy data')
    parser.add_argument('-pset', required=True,
                        help='positive set of Chip-Chip sequences in FASTA format')
    parser.add_argument('-nset', required=True,
                        help='negative set of Chip-Chip sequences in FASTA format')
    parser.add_argument('-plist', required=True,
                        help='pattern list file')
    parser.add_argument('-seqtype', required=True, choices=['dna', 'rna'],
                        help='the sequence type of the patterns [dna|rna]')
    parser.add_argument('-out', required=True,
                        help='output directory')
    # parser.add_argument('-cpu', type=int, default=1,
    #                     help='Number of CPUs to perform the analysis (default: 1)')
    parser.add_argument('-sp', type=int, default=1,
                        help='the weight of position scoring (default: 1)')
    parser.add_argument('-cluster', type=int, default=5,
                        help='maximum number of clusters in the output (default: 5)')
    parser.add_argument('-n', type=int, default=5,
                        help='maximum number of patterns per cluster (default: 5)')
    parser.add_argument('-log',
                        help='log file (default: stdout)')
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

    pattern_set = PatternSet(reverse_complement)

    logger.info('building match tables of patterns')
    with open(args.plist, 'r') as fi:
        for line in fi:
            pattern = Pattern(line.strip())
            pattern.build_matchtable_pset(parse_fasta_noheader(args.pset), reverse_complement)
            pattern.build_matchtable_nset(parse_fasta_noheader(args.nset), reverse_complement)
            pattern_set.add(pattern)

    pattern_scoring = PatternScoring(args.sp)
    pattern_scoring.build(pattern_set)

    cluster = Cluster(args.cluster, 0.8, args.n, reverse_complement)
    cluster.run(pattern_scoring)

    cluster_pfm = {}
    logger.info('merging patterns and calculating PFMs')

    with open(os.path.join(args.out, 'clustered_patterns.txt'), 'w') as fo, \
            open(os.path.join(args.out, 'merged_patterns.txt'), 'w') as fo2, \
            open(os.path.join(args.out, 'match_sequences.txt'), 'w') as fo3:
        fo.write('cluster: pattern\n')
        fo2.write('cluster: strand pattern\n')
        fo3.write('cluster: sequence\n')

        for i, j in cluster.results.iteritems():
            for p in j:
                fo.write('{0}: {1}\n'.format(i, p.sequence))
                fo.flush()

            merged = merge_patterns(j, reverse_complement)
            match_sequences = merged.extract_match_sequences(args.pset)

            for p in merged.patterns:
                fo2.write('{0}: {1} {2}\n'.format(i, merged._strands.get(p), p.sequence))
                fo2.flush()

            for s in match_sequences:
                fo3.write('{0}: {1}\n'.format(i, s))
                fo3.flush()

            cluster_pfm.update({i: pfm(match_sequences)})

    for i, j in cluster_pfm.iteritems():
        with open(os.path.join(args.out, 'cluster_{0}.pfm.txt'.format(i)), 'w') as fo:
            fo.write('cluster_{0}\t{1}\n'.format(i, str(len(j.get('a')))))
            for base in ['a', 't', 'c', 'g']:
                fo.write('\t'.join([str(i) for i in j.get(base)]))
                fo.write('\n')
                fo.flush()

    logger.info('Job has finished.')


if __name__ == '__main__':
    main()
