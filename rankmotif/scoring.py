import logging
from math import sqrt
from .basic import PatternSet


class PreferentialOccurrence(object):
    """Calculate z-score of the pattern(s) in the positive set
    and the negative set"""

    def __init__(self):
        self._logger = logging.getLogger(self.__class__.__name__)
        self.results = {}

    def build(self, pattern_set, append=False):
        assert isinstance(pattern_set, PatternSet)

        self._logger.info('building scores')

        if not append:
            self.results = {}

        for pattern in pattern_set:
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

        def add(self, seqid, indices):
            for index in indices:
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
        self._logger = logging.getLogger(self.__class__.__name__)
        self.results = {}
        self._ntscore = self._PositionScoreMatrix()

    def build(self, pattern_set, append=False, seqmask=False):
        assert isinstance(pattern_set, PatternSet)

        self._logger.info('building scores')

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
                if seqmask and seqid > pattern.matchtable_pset.n_hitseqs:
                    continue
                for index in indices:
                    for i in index:
                        match_score += self._ntscore.get(seqid, i)

            score = float(match_score) / (pattern.matchtable_pset.n_hitseqs * pattern.n_nonwildcards)
            self.results.update({pattern: score})

        return self


class NucleosomeOccupancyScoring(object):
    """Calculate the nucleosome occupancy scores of the pattern(s)"""

    def __init__(self):
        self._logger = logging.getLogger(self.__class__.__name__)
        self.results = {}

    def build(self, pattern_set, score_file=None, append=False):
        assert isinstance(pattern_set, PatternSet)

        self._logger.info('building scores')

        if not append:
            self.results = {}

        scores = parse_base_score(score_file, scale=True)
        for pattern in pattern_set:
            s = []
            for seqid, indices in pattern.matchtable_pset.pos_matches.iteritems():
                if scores.get(seqid):
                    for index in indices:
                        for i in index:
                            s.append(1 - scores.get(seqid).get(i))
            self.results.update({pattern: float(sum(s)) / len(s)})

        return self


class ConservationScoring(object):
    """Calculate the conservation scores of the pattern(s)"""

    def __init__(self):
        self._logger = logging.getLogger(self.__class__.__name__)
        self.results = {}

    def build(self, pattern_set, score_file=None, append=False):
        assert isinstance(pattern_set, PatternSet)

        self._logger.info('building scores')

        if not append:
            self.results = {}

        scores = parse_base_score(score_file, scale=False)
        for pattern in pattern_set:
            s = []
            for seqid, indices in pattern.matchtable_pset.pos_matches.iteritems():
                if scores.get(seqid):
                    for index in indices:
                        for i in index:
                            if scores.get(seqid).get(i):
                                s.append(scores.get(seqid).get(i))
            self.results.update({pattern: float(sum(s)) / len(s)})

        return self


class PatternScoring(object):

    def __init__(self, sp_weight=1, sn_weight=1, sc_weight=1):
        self.results = {}
        self.sp_weight = sp_weight
        self.sn_weight = sn_weight
        self.sc_weight = sc_weight
        self._poccur = PreferentialOccurrence()
        self._pscore = PositionScoring()
        self._noscore = NucleosomeOccupancyScoring()

    def build(self, pattern_set, append=False, seqmask=False, nuclocc=None,
              consv=None):
        assert isinstance(pattern_set, PatternSet)

        if not append:
            self.results = {}
            self._poccur = PreferentialOccurrence()
            self._pscore = PositionScoring()
            self._noscore = NucleosomeOccupancyScoring()
            self._csscore = ConservationScoring()

        self._poccur.build(pattern_set, append)
        self._pscore.build(pattern_set, append, seqmask)
        if nuclocc:
            self._noscore.build(pattern_set, score_file=nuclocc, append=append)
        if consv:
            self._csscore.build(pattern_set, score_file=consv, append=append)

        for pattern, score in self._poccur.results.iteritems():
            pattern_score = score
            pattern_score *= self._pscore.results.get(pattern) ** self.sp_weight
            if nuclocc:
                pattern_score *= self._noscore.results.get(pattern) ** self.sn_weight
            if consv:
                pattern_score *= self._csscore.results.get(pattern) ** self.sc_weight

            self.results.update({pattern: pattern_score})

        return self


def parse_base_score(fpath, scale=False):
    """Parse the base-based score data"""
    if isinstance(fpath, str):
        fi = open(fpath, 'r')
    else:
        return None

    scores = {}
    for line in fi:
        data = line.strip().split('\t')
        seqid = int(data[0])
        index = int(data[1]) - 1
        score = float(data[4])
        if seqid in scores:
            scores.get(seqid).update({index: score})
        else:
            scores.update({seqid: {index: score}})

    fi.close()

    # Scale score range to [0, 1]
    if scale:
        for seqid, score in scores.iteritems():
            values = score.values()
            min_score = min(values)
            max_score = max(values)
            interval = max_score - min_score
            for i, j in score.iteritems():
                score[i] = (j - min_score) / interval

    return scores
