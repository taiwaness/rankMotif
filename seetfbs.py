#!/usr/bin/env python

import os
import argparse
from seetfbs import ranking


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-ps', required=True, type=argparse.FileType('r'),
                        help='Positive set of Chip-Chip sequences in FASTA format')
    parser.add_argument('-ns', required=True, type=argparse.FileType('r'),
                        help='Negative set of Chip-Chip sequences in FASTA format')
    parser.add_argument('-plist', required=True, type=argparse.FileType('r'),
                        help='Pattern list file')
    parser.add_argument('-seqtype', required=True, choices=['dna', 'rna'],
                        help='The sequence type of the patterns [dna|rna]')
    parser.add_argument('-out', required=True,
                        help='Output directory')
    parser.add_argument('-cpu', type=int, default=1,
                        help='Number of CPUs to perform the analysis (default: 1)')
    parser.add_argument('-nc', type=int, default=5,
                        help='Number of clusters in the output (default: 5)')
    args = parser.parse_args()

    if not os.path.exists(args.o):
        os.mkdir(args.o)

if __name__ == '__main__':
    main()
