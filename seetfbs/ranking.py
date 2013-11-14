from math import sqrt
from basic import PatternSet


class PreferentialOccurrence(object):
    """Calculate z-score of the pattern(s) in the positive set
    and the negative set"""

    def __init__(self):
        self.results = {}

    def build(self, pattern_set, append=False):
        assert isinstance(pattern_set, PatternSet)

        if not append:
            self.results = {}

        for pattern in PatternSet:
            mt_p = pattern.matchtable_pset
            mt_n = pattern.matchtable_nset

            fp = float(mt_p.n_hitseqs) / mt_p.n_seqs
            fn = float(mt_n.n_hitseqs) / mt_n.n_seqs
            f = float(mt_p.n_hitseqs * fp + mt_n.n_hitseqs * fn) / (mt_p.n_hitseqs + mt_n.n_hitseqs)
            z_score = (fp - fn) / sqrt(f * (1 - f) * (1.0 / mt_p.n_hitseqs + 1.0 / mt_n.n_hitseqs))

            self.results.update({pattern: z_score})

        return self


class PositionScoring(object):

    class PositionScoreMatrix(object):

        def __init__(self, weight=1):
            self.matrix = {}
            self.weight = weight

        def add(self, n_seqs, index):
            if n_seqs in self.matrix:
                if index in self.matrix.get(n_seqs):
                    self.matrix.get(n_seqs)[index] += self.weight
                else:
                    self.matrix.get(n_seqs).update({index: self.weight})
            else:
                self.matrix.update({n_seqs: {index: self.weight}})

        def get(self, n_seqs, index):
            if n_seqs not in self.matrix or index not in self.matrix.get(n_seqs):
                return 0
            else:
                return self.matrix.get(n_seqs).get(index)

    def __init__(self, psm_weight=1):
        self.results = {}
        self.psm_weight = psm_weight
        self._ntscore = self.PositionScoreMatrix(weight=psm_weight)

    def build(self, pattern_set, append=False):
        assert isinstance(pattern_set, PatternSet)

        if not append:
            self._ntscore = self.PositionScoreMatrix(weight=self.psm_weight)

        for pattern in pattern_set:
            for n_pset, indices in pattern.matchtable_pset.pos_matches.iteritems():
                for index in indices:
                    self._ntscore.add(n_pset, index)

        for pattern in pattern_set:
            match_score = 0
            for n_pset, indices in pattern.matchtable_pset.pos_matches.iteritems():
                for index in indices:
                    match_score += self._ntscore.get(n_pset, index)

            score = float(match_score) / (pattern.matchtable_pset.n_seqs * pattern.n_nonwildcards)
            self.results.update({pattern: score})

        return self


class PatternScoring(object):

    def __init__(self, sp_weight=1):
        self.results = {}
        self.sp_weight = sp_weight
        self._poccur = PreferentialOccurrence()
        self._pscore = PositionScoring()

    def build(self, pattern_set, append=False):
        assert isinstance(pattern_set, PatternSet)

        if not append:
            self.results = {}
            self._poccur = PreferentialOccurrence()
            self._pscore = PositionScoring()

        self._poccur.build(pattern_set, append)
        self._pscore.build(pattern_set, append)

        for pattern, score in self._poccur.results.iteritems():
            pattern_score = score * (self._pscore.get(pattern) ** self.sp_weight)
            self.results.update({pattern: pattern_score})

        return self
