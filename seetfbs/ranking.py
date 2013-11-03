import re
from math import sqrt
from seqio import revcomp


class Pattern(object):

    def __init__(self, sequence, reverse=True):
        self._reverse = reverse
        self.sequence = sequence
        self.num_wildcard = sequence.count('n')
        self.num_nonwildcard = len(sequence) - self.num_wildcard
        if reverse:
            self._re = re.compile('%s|%s' % (self.sequence.replace('n', '[atcg]'),
                                             revcomp(self.sequence).replace('n', '[atcg]')),
                                  re.IGNORECASE)
        else:
            self._re = re.compile(self.sequence.replace('n', '[atcg]'), re.IGNORECASE)
        self.ranking_score = None
        self.sd = None
        self.sp = None
        self._sp_match = 0
        self._sp_match_index = []
        self._sp_match_score = 0
        self.sc = None

    def finditer(self, sequence):
        return self._re.finditer(sequence)

    def preferential_occurrence(self, pset, nset):
        num_ps = 0
        num_ns = 0
        total = 0

        for i in pset:
            total += 1
            if self.re.search(i):
                num_ps += 1

        for i in nset:
            total += 1
            if self.re.search(i):
                num_ns += 1

        if num_ps == 0 or num_ns == 0:
            return None

        fp = float(num_ns) / total
        fn = float(num_ps) / total
        f = float(num_ps * fp + num_ns * fn) / (num_ps + num_ns)
        z_score = (fp - fn) / sqrt(f * (1 - f) * (1.0 / num_ps + 1.0 / num_ns))

        self.sd = z_score

    def pfm(self, pset):
        """Claculate position frequency matrix (PFM)"""
        matrix = {
            'A': [0] * len(self.sequence),
            'T': [0] * len(self.sequence),
            'C': [0] * len(self.sequence),
            'G': [0] * len(self.sequence),
        }
        for i in pset:
            for match in self.finditer(i):
                subseq = i[match.start(): match.end()].upper()
                for i, j in enumerate(subseq):
                    matrix.get(j)[i] += 1

        return matrix


class PatternCollect(object):
    """A collection of unique pattern objects for
    pattern ranking"""

    def __init__(self):
        self._collect = {}

    def __iter__(self):
        return self._collect.itervalues()

    def __len__(self):
        return len(self._collect)

    def add(self, pattern):
        """Add a pattern object"""
        assert isinstance(pattern, Pattern)

        if self._reverse and revcomp(pattern.sequence) in self._collect:
            self._collect.update({revcomp(pattern.sequence): pattern})
        else:
            self._collect.update({pattern.sequence: pattern})

    def preferential_occurrence(self, pset, nset):
        for i in self:
            i.preferential_occurrence(pset, nset)

    def position_scoring(self, pset):
        for sequence in pset:
            ntscore = {i: 0 for i in range(len(sequence))}
            for pattern in self:
                has_match = False
                for match in pattern.finditer(sequence):
                    has_match = True
                    for i, j in enumerate(pattern.sequence):
                        if j != 'n':
                            ntscore[match.start() + i] += 1
                            pattern._sp_match_index.append(match.start() + i)
                if has_match:
                    pattern._sp_match += 1
            for pattern in self:
                for i in pattern._sp_match_index:
                    pattern._sp_match_score += ntscore.get(i)
                pattern._sp_match_index = []
        for pattern in self:
            pattern.sp = float(pattern._sp_match_score) / (pattern._sp_match * pattern.num_nonwildcard)

    def pattern_scoring(self, sp_weight=1):
        for pattern in self:
            pattern.score = pattern.sd * (pattern.sp ** sp_weight)

    def clustering(self, max_cluster=5, similarity=0.8):
        ranked_patterns = sorted(self, key=lambda x: x.score, reverse=True)
        clusters = [0] * len(self)
        clustered_patterns = {}
        num_cluster = 0

        for i in len(self):
            if num_cluster == max_cluster:
                break
            if clusters[i] == 0:
                num_cluster += 1
                clusters[i] = num_cluster
                clustered_patterns.update({num_cluster: [ranked_patterns[i]]})
                for j in range(i + 1, len(self)):
                    if position_weight_matrix(ranked_patterns[i], ranked_patterns[j]) >= similarity:
                        clusters[j] = num_cluster
                        clustered_patterns.get(i).append(ranked_patterns[j])


def position_weight_matrix(pfm1, pfm2):
    pass
