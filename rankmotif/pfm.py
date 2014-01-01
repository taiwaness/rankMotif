import argparse
from math import sqrt
from basic import Pattern
from seqio import revcomp


def pfm(pattern_sequence):
    """Claculate position frequency matrix (PFM) of matching sequences"""
    if isinstance(pattern_sequence, Pattern):
        sequences = []
        for match_sequences in pattern_sequence.matchtable_pset.match_sequences.itervalues():
            for strand, sequence in match_sequences:
                if strand == 2:
                    sequences.append(revcomp(sequence))
                else:
                    sequences.append(sequence)
    else:
        sequences = pattern_sequence

    ncol = len(sequences[0])
    matrix = {
        'a': [0] * ncol,
        't': [0] * ncol,
        'c': [0] * ncol,
        'g': [0] * ncol,
    }
    total = [0] * ncol

    for s in sequences:
        for i, j in enumerate(s):
            matrix.get(j)[i] += 1
            total[i] += 1

    # Normalization
    for i in xrange(ncol):
        matrix.get('a')[i] = float(matrix.get('a')[i]) / total[i]
        matrix.get('t')[i] = float(matrix.get('t')[i]) / total[i]
        matrix.get('c')[i] = float(matrix.get('c')[i]) / total[i]
        matrix.get('g')[i] = float(matrix.get('g')[i]) / total[i]

    return matrix


def simpfm(pfm_1, pfm_2, gc_content, max_wsize=7, reverse_complement=False):
    """Calculate the similarity scores of two PFMs and return the top one"""
    len_pfm_1 = len(pfm_1.get('a'))
    len_pfm_2 = len(pfm_2.get('a'))

    assert max_wsize > 0

    min_pfm_len = min(len_pfm_1, len_pfm_2)
    if max_wsize > min_pfm_len:
        max_wsize = min_pfm_len - 2
    alnlen = len_pfm_1 + len_pfm_2 - max_wsize
    max_score = 0

    for i in xrange(alnlen - len_pfm_1 + 1):
        for j in xrange(alnlen - len_pfm_2 + 1):
            if j > i:
                left_1 = 0
                left_2 = j - i
            elif i > j:
                left_1 = i - j
                left_2 = 0
            else:
                left_1 = 0
                left_2 = 0

            if len_pfm_2 + j - 1 > len_pfm_1 + i - 1:
                right_1 = (len_pfm_2 + j - 1) - (len_pfm_1 + i - 1)
                right_2 = 0
            elif len_pfm_1 + i - 1 > len_pfm_2 + j - 1:
                right_1 = 0
                right_2 = (len_pfm_1 + i - 1) - (len_pfm_2 + j - 1)
            else:
                right_1 = 0
                right_2 = 0

            p1 = expand_pfm(pfm_1, left_1, right_1, gc_content)
            p2 = expand_pfm(pfm_2, left_2, right_2, gc_content)
            score = simpfm_scoring(p1, p2)
            if score > max_score:
                max_score = score

    if reverse_complement:
        pfm_2 = reverse_pfm(pfm_2)

        for i in xrange(alnlen - len_pfm_1 + 1):
            for j in xrange(alnlen - len_pfm_2 + 1):
                if j > i:
                    left_1 = 0
                    left_2 = j - i
                elif i > j:
                    left_1 = i - j
                    left_2 = 0
                else:
                    left_1 = 0
                    left_2 = 0

                if len_pfm_2 + j - 1 > len_pfm_1 + i - 1:
                    right_1 = (len_pfm_2 + j - 1) - (len_pfm_1 + i - 1)
                    right_2 = 0
                elif len_pfm_1 + i - 1 > len_pfm_2 + j - 1:
                    right_1 = 0
                    right_2 = (len_pfm_1 + i - 1) - (len_pfm_2 + j - 1)
                else:
                    right_1 = 0
                    right_2 = 0

                p1 = expand_pfm(pfm_1, left_1, right_1, gc_content)
                p2 = expand_pfm(pfm_2, left_2, right_2, gc_content)
                score = simpfm_scoring(p1, p2)
                if score > max_score:
                    max_score = score

        return [(1, 'pfm_1'), (2, 'pfm_2'), max_score]
    else:
        return [(1, 'pfm_1'), (1, 'pfm_2'), max_score]


def expand_pfm(pfm, left_length, right_length, gc_content):
    at = (1 - gc_content) / 2
    cg = gc_content / 2
    expanded_pfm = {
        'a': ([at] * left_length) + pfm.get('a') + ([at] * right_length),
        't': ([at] * left_length) + pfm.get('t') + ([at] * right_length),
        'c': ([cg] * left_length) + pfm.get('c') + ([cg] * right_length),
        'g': ([cg] * left_length) + pfm.get('g') + ([cg] * right_length),

    }

    return expanded_pfm


def simpfm_scoring(pfm_1, pfm_2):
    scores = []
    alnlen = len(pfm_1.get('a'))
    for i in xrange(alnlen):
        s = 0
        for j in ['a', 't', 'c', 'g']:
            s += (pfm_1.get(j)[i] - pfm_2.get(j)[i]) ** 2
        scores.append(1 - sqrt(s) / sqrt(2))

    return sum(scores) / alnlen


def reverse_pfm(pfm):
    rv_pfm = {}
    rv_pfm.update({'a': pfm.get('t')[::-1]})
    rv_pfm.update({'t': pfm.get('a')[::-1]})
    rv_pfm.update({'c': pfm.get('g')[::-1]})
    rv_pfm.update({'g': pfm.get('c')[::-1]})

    return rv_pfm


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pfm1')
    parser.add_argument('pfm2')
    parser.add_argument('-gc', type=float, required=True, metavar='<float>',
                        help='GC contents')
    parser.add_argument('-r', action='store_true')
    parser.add_argument('-ws', type=int, default=7, metavar='<int>',
                        help='maximum window size of simpfm (default: 7)')
    args = parser.parse_args()

    pfm1 = {}
    with open(args.pfm1, 'r') as fi:
        gene, alnlen = fi.readline().split('\t')
        alnlen = int(alnlen)
        for i in ['a', 't', 'c', 'g']:
            pfm1.update({i: [float(x) for x in fi.readline().split('\t')]})

    pfm2 = {}
    with open(args.pfm2, 'r') as fi:
        gene, alnlen = fi.readline().split('\t')
        alnlen = int(alnlen)
        for i in ['a', 't', 'c', 'g']:
            pfm2.update({i: [float(x) for x in fi.readline().split('\t')]})

    simpfm_score = simpfm(pfm1, pfm2, args.gc, max_wsize=args.ws, reverse_complement=args.r)
    print(simpfm_score)


if __name__ == '__main__':
    main()
