import os
import logging
import subprocess
from tempfile import NamedTemporaryFile
from .scoring import PatternScoring
from .seqio import revcomp
from .basic import Pattern


class Cluster(object):

    def __init__(self, max_cluster=5, similarity=0.8,
                 max_patterns_per_cluster=5, reverse_complement=False):
        self._logger = logging.getLogger(self.__class__.__name__)
        self.max_cluster = max_cluster
        self.similarity = similarity
        self.max_patterns_per_cluster = max_patterns_per_cluster
        self.reverse_complement = reverse_complement
        self.results = {}

    def run(self, pattern_scoring):
        isinstance(pattern_scoring, PatternScoring)
        self._logger.info('begin running pattern clustering')

        self.results = {}
        ranked_patterns = sorted(
            pattern_scoring.results, key=lambda x: pattern_scoring.results.get(x), reverse=True)
        clusters = [0] * len(pattern_scoring.results)
        n_clusters = 0
        n_patterns = {}

        for i in range(len(ranked_patterns)):
            if n_clusters == self.max_cluster:
                break
            if clusters[i] == 0:
                n_clusters += 1
                clusters[i] = n_clusters
                n_patterns.update({n_clusters: 0})
                self.results.update({n_clusters: [ranked_patterns[i]]})
                for j in range(i + 1, len(ranked_patterns)):
                    if n_patterns.get(n_clusters) == self.max_patterns_per_cluster - 1:
                        break
                    score = sim_pfm(
                        pfm(ranked_patterns[i]), pfm(ranked_patterns[j]), self.reverse_complement)[2]
                    if clusters[j] == 0 and score >= self.similarity:
                        clusters[j] = n_clusters
                        self.results.get(n_clusters).append(ranked_patterns[j])
                        n_patterns[n_clusters] += 1

        return self


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
    for i in range(ncol):
        matrix.get('a')[i] = float(matrix.get('a')[i]) / total[i]
        matrix.get('t')[i] = float(matrix.get('t')[i]) / total[i]
        matrix.get('c')[i] = float(matrix.get('c')[i]) / total[i]
        matrix.get('g')[i] = float(matrix.get('g')[i]) / total[i]

    return matrix


def sim_pfm(pfm_1, pfm_2, reverse_complement=False):
    """Calculate the similarities of two PFMs"""

    f_pfm_1 = NamedTemporaryFile()
    f_pfm_1.write('\t'.join(['pfm_1', str(len(pfm_1.get('a')))]))
    f_pfm_1.write('\n')
    for i in ['a', 't', 'c', 'g']:
        f_pfm_1.write('\t'.join([str(s) for s in pfm_1.get(i)]))
        f_pfm_1.write('\n')
    f_pfm_1.flush()

    f_pfm_2 = NamedTemporaryFile()
    f_pfm_2.write('\t'.join(['pfm_2', str(len(pfm_2.get('a')))]))
    f_pfm_2.write('\n')
    for i in ['a', 't', 'c', 'g']:
        f_pfm_2.write('\t'.join([str(s) for s in pfm_2.get(i)]))
        f_pfm_2.write('\n')
    f_pfm_2.flush()

    cmd = [
        'perl',
        os.path.join(os.path.dirname(__file__), 'thirdparty/calcPwm.pl'),
        f_pfm_1.name,
        f_pfm_2.name,
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()

    max_score = 0
    best_pfm1 = None
    best_pfm2 = None

    for i in stdout.split('\n'):
        if i == '':
            continue
        data = i.split('\t')
        if not reverse_complement and ('rv' in data[0] or 'rv' in data[1]):
            continue
        similarity = float(data[2])
        if similarity > max_score:
            max_score = similarity
            best_pfm1, best_pfm2 = data[0], data[1]

    if reverse_complement and 'rv' in best_pfm1:
        pfm1_result = (2, 'pfm_1')
    else:
        pfm1_result = (1, 'pfm_1')

    if reverse_complement and 'rv' in best_pfm2:
        pfm2_result = (2, 'pfm_2')
    else:
        pfm2_result = (1, 'pfm_2')

    result = [pfm1_result, pfm2_result, max_score]

    return result
