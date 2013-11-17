#!/usr/bin/env python

import os
import argparse
from seetfbs.basic import Pattern, PatternSet
from seetfbs.scoring import PatternScoring
from seetfbs.ranking import Cluster
from seetfbs.seqio import parse_fasta_noheader


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-pset', required=True,
                        help='Positive set of Chip-Chip sequences in FASTA format')
    parser.add_argument('-nset', required=True,
                        help='Negative set of Chip-Chip sequences in FASTA format')
    parser.add_argument('-plist', required=True,
                        help='Pattern list file')
    parser.add_argument('-seqtype', required=True, choices=['dna', 'rna'],
                        help='The sequence type of the patterns [dna|rna]')
    parser.add_argument('-out', required=True,
                        help='Output directory')
    parser.add_argument('-cpu', type=int, default=1,
                        help='Number of CPUs to perform the analysis (default: 1)')
    parser.add_argument('-sp', type=int, default=1,
                        help='The weight of position scoring (default: 1)')
    parser.add_argument('-cluster', type=int, default=5,
                        help='Maximum number of clusters in the output (default: 5)')
    args = parser.parse_args()

    if not os.path.exists(args.out):
        os.makedirs(args.out)

    if args.seqtype == 'dna':
        reverse_complement = True
    else:
        reverse_complement = False

    pattern_set = PatternSet(reverse_complement)

    with open(args.plist, 'r') as fi:
        for line in fi:
            pattern = Pattern(line.strip())
            pattern.build_matchtable_pset(parse_fasta_noheader(args.pset), reverse_complement)
            pattern.build_matchtable_nset(parse_fasta_noheader(args.nset), reverse_complement)
            pattern_set.add(pattern)

    pattern_scoring = PatternScoring(args.sp)
    pattern_scoring.build(pattern_set)

    cluster = Cluster(args.cluster, similarity=0.8, reverse_complement=reverse_complement)
    cluster.run(pattern_scoring)

    with open(os.path.join(args.out, 'ranked_patterns.txt'), 'w') as fo:
        for i, j in cluster.clustered_patterns.iteritems():
            for p in j:
                fo.write('{0}: {1}'.format(i, p.sequence))
                fo.write('\n')
            fo.flush()


if __name__ == '__main__':
    main()
