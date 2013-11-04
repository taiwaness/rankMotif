from match import sqrt
from basic import Pattern, MatchTable, PatternCollection


class PreferentialOccurrence(object):
    """Calculate z-score of the pattern(s) in the positive set
    and the negative set"""

    def __init__(self):
        self.results = {}

    def run(self, pattern, append=False):
        assert isinstance(pattern, Pattern)
        assert isinstance(pattern.matchtable, MatchTable)

        if not append:
            self.results = {}

        mtable = pattern.matchtable
        fp = float(mtable.pmatch) / mtable.pnum
        fn = float(mtable.nmatch) / mtable.nnum
        f = float(mtable.pmatch * fp + mtable.nmatch * fn) / (mtable.pmatch + mtable.nmatch)
        z_score = (fp - fn) / sqrt(f * (1 - f) * (1.0 / mtable.pmatch + 1.0 / mtable.nmatch))

        self.results.update({pattern: z_score})

    def batch_run(self, pattern_collection, append=False):
        assert isinstance(pattern_collection, PatternCollection)

        if not append:
            self.results = {}

        for i in pattern_collection:
            self.run(i)


class PositionScoring(object):

    def __init__(self, psm_weight=1):
        self.results = {}
        self.psm_weight = psm_weight
        self._ntscore = PositionScoreMatrix(weight=psm_weight)

    def run(self, pattern_collection, append=False):
        assert isinstance(pattern_collection, PatternCollection)

        if not append:
            self._ntscore = PositionScoreMatrix(weight=self.psm_weight)
        for i in pattern_collection:
            for pnum, index in i.matchtable.pindex.match.iteritems():
                for j in index:
                    self._ntscore.add(pnum, j)

        for i in pattern_collection:
            match_score = 0
            for pnum, index in i.matchtable.pindex.match.iteritems():
                for j in index:
                    match_score += self._ntscore.get(pnum, index)

            score = float(match_score) / (i.matchtable.pnum * i.num_nonwildcard)
            self.results.update({i: score})


class PositionScoreMatrix(object):

    def __init__(self, weight=1):
        self.matrix = {}
        self.weight = weight

    def add(self, pnum, index):
        if pnum in self.matrix:
            if index in self.matrix.get(pnum):
                self.matrix.get(pnum)[index] += self.weight
            else:
                self.matrix.get(pnum).update({index: self.weight})
        else:
            self.matrix.update({pnum: {index: self.weight}})

    def get(self, pnum, index):
        if pnum not in self.matrix or index not in self.matrix.get(pnum):
            return 0
        else:
            return self.matrix.get(pnum).get(index)


class PatternScoring(object):

    def __init__(self, sp_weight=1):
        self.results = {}
        self.sp_weight = sp_weight

    def run(self, pattern_collection, append=False):
        assert isinstance(pattern_collection, PatternCollection)
        if not append:
            self.results = {}
        prefoccur = PreferentialOccurrence()
        pscoring = PositionScoring()
        prefoccur.batch_run(pattern_collection, append)
        pscoring.run(pattern_collection, append)

        for i, j in prefoccur.results.iteritems():
            score = j * (pscoring ** self.sp_weight)
            self.results.update({i: score})
