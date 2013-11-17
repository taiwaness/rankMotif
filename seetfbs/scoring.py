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
            f = float(mt_p.n_seqs * fp + mt_n.n_seqs * fn) / (mt_p.n_seqs + mt_n.n_seqs)
            z_score = (fp - fn) / sqrt(f * (1 - f) * (1.0 / mt_p.n_seqs + 1.0 / mt_n.n_seqs))

            self.results.update({pattern: z_score})

        return self


class PositionScoring(object):

    class _PositionScoreMatrix(object):

        def __init__(self):
            self.matrix = {}

        def add(self, seqid, index):
            if seqid in self.matrix:
                if index in self.matrix.get(seqid):
                    self.matrix.get(seqid)[index] += 1
                else:
                    self.matrix.get(seqid).update({index: 1})
            else:
                self.matrix.update({seqid: {index: 1}})

        def get(self, seqid, index):
            if seqid not in self.matrix or index not in self.matrix.get(seqid):
                return 0
            else:
                return self.matrix.get(seqid).get(index)

    def __init__(self):
        self.results = {}
        self._ntscore = self._PositionScoreMatrix()

    def build(self, pattern_set, append=False):
        assert isinstance(pattern_set, PatternSet)

        if not append:
            self._ntscore = self._PositionScoreMatrix()

        # Calculate accumulative scores
        for pattern in pattern_set:
            for seqid, indices in pattern.matchtable_pset.pos_nonwildcards.iteritems():
                for index in indices:
                    self._ntscore.add(seqid, index)

        # Assign scores to each pattern
        for pattern in pattern_set:
            match_score = 0
            for seqid, indices in pattern.matchtable_pset.pos_nonwildcards.iteritems():
                for index in indices:
                    match_score += self._ntscore.get(seqid, index)

            score = float(match_score) / (pattern.matchtable_pset.n_hitseqs * pattern.n_nonwildcards)
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
