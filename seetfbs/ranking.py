import os
import subprocess
from tempfile import NamedTemporaryFile
from scoring import PatternScoring
from seqio import revcomp
from basic import Pattern


class Cluster(object):

    def __init__(self, max_cluster=5, similarity=0.8):
        self.max_cluster = max_cluster
        self.similarity = similarity
        self.clustered_patterns = {}

    def run(self, pattern_scoring):
        isinstance(pattern_scoring, PatternScoring)

        self.clustered_patterns = {}
        ranked_patterns = sorted(
            pattern_scoring.results, key=lambda x: pattern_scoring.results.get(x), reverse=True)
        clusters = [0] * len(pattern_scoring.results)
        n_clusters = 0

        for i in range(len(ranked_patterns)):
            if n_clusters == self.max_cluster:
                break
            if clusters[i] == 0:
                n_clusters += 1
                clusters[i] = n_clusters
                self.clustered_patterns.update({n_clusters: [ranked_patterns[i]]})
                for j in range(i + 1, len(ranked_patterns)):
                    if clusters[j] == 0 and sim_pfm(pfm(ranked_patterns[i]), pfm(ranked_patterns[j])) >= self.similarity:
                        clusters[j] = n_clusters
                        self.clustered_patterns.get(i).append(ranked_patterns[j])

        return self


def pfm(pattern):
    """Claculate position frequency matrix (PFM)
    of matching sequences of a pattern"""
    assert isinstance(pattern, Pattern)

    sequences = []
    for strand, sequence in pattern.matchtable_pset.match_sequences.itervalues():
        if strand == 2:
            sequences.append(revcomp(sequence))
        else:
            sequences.append(sequence)

    ncol = max([len(i) for i in sequences])
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


def sim_pfm(pfm_1, pfm_2):
    """Calculate the similarities of two PFMs"""

    f_pfm_1 = NamedTemporaryFile()
    f_pfm_1.write('\t'.join('pfm_1', str(len(pfm_1.get('a')))))
    f_pfm_1.write('\n')
    for i in ['a', 't', 'c', 'g']:
        f_pfm_1.write('\t'.join(pfm_1.get(i)))
        f_pfm_1.write('\n')
    f_pfm_1.flush()

    f_pfm_2 = NamedTemporaryFile()
    f_pfm_2.write('\t'.join('pfm_2', str(len(pfm_2.get('a')))))
    f_pfm_2.write('\n')
    for i in ['a', 't', 'c', 'g']:
        f_pfm_2.write('\t'.join(pfm_2.get(i)))
        f_pfm_2.write('\n')
    f_pfm_2.flush()

    cmd = [
        'perl',
        os.path.join(__file__, 'thirdparty/calcPwm.pl'),
        f_pfm_1.name,
        f_pfm_2.name,
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()

    max_score = 0
    best_pfm1 = best_pfm2 = None
    for i in stdout.split('\n'):
        data = i.split('\t')
        similarity = float(data[2])
        if similarity > max_score:
            max_score = similarity
            best_pfm1, best_pfm2 = data[0], data[1]

    if 'rv' in best_pfm1:
        pfm1_result = (2, 'pfm_1')
    else:
        pfm1_result = (1, 'pfm_1')

    if 'rv' in best_pfm2:
        pfm2_result = (2, 'pfm_2')
    else:
        pfm2_result = (1, 'pfm_2')

    result = [pfm1_result, pfm2_result, max_score]

    return result
