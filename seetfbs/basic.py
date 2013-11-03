import re
from seqio import revcomp
from matrix import position_weight_matrix


class Pattern(object):

    def __init__(self, sequence):
        self.sequence = sequence
        self.num_wildcard = sequence.count('n')
        self.num_nonwildcard = len(sequence) - self.num_wildcard
        self.matchtable = None

    def build(self, pset, nset, reverse_complement=False, append=False):
        self.matchtable = MatchTable(self, reverse_complement)
        self.matchtable.index_pset(pset, append)
        self.matchtable.index_nset(nset, append)

    # def pfm(self, pset):
    #     """Claculate position frequency matrix (PFM)"""
    #     matrix = {
    #         'A': [0] * len(self.sequence),
    #         'T': [0] * len(self.sequence),
    #         'C': [0] * len(self.sequence),
    #         'G': [0] * len(self.sequence),
    #     }
    #     for i in pset:
    #         for match in self.finditer(i):
    #             subseq = i[match.start(): match.end()].upper()
    #             for i, j in enumerate(subseq):
    #                 matrix.get(j)[i] += 1

    #     return matrix


class MatchTable(object):
    """Index of a pattern against matching sequences."""

    def __init__(self, pattern, reverse_complement=False):
        assert isinstance(pattern, Pattern)
        self.pattern = pattern
        self.reverse_complement = reverse_complement
        self.pmatch = 0
        self.nmatch = 0
        self.pindex = PatternPositionIndex()
        self.nindex = PatternPositionIndex()
        self.pnum = 0
        self.nnum = 0

    def index_pset(self, pset, append=False):
        if not append:
            self.pmatch = 0
            self.pindex = PatternPositionIndex()
            self.pnum = 0
        p = re.compile(self.pattern.sequence.replace('n', '[atcg]'),
                       re.IGNORECASE)
        prc = re.compile('%s' % (revcomp(self.pattern.sequence).replace('n', '[atcg]')),
                         re.IGNORECASE)

        for i in pset:
            self.pnum += 1
            has_match = False
            if self.reverse_complement:
                for j in prc.finditer(i):
                    has_match = True
                    self.pindex.append(self.pnum,
                                       self.pattern.sequence,
                                       i,
                                       j.start(),
                                       j.end() - 1,
                                       self.reverse_complement)
            else:
                for j in p.finditer(i):
                    self.pindex.append(self.pnum,
                                       self.pattern.sequence,
                                       i,
                                       j.start(),
                                       j.end() - 1,
                                       self.reverse_complement)

            if has_match:
                self.pmatch += 1

    def index_nset(self, nset, append=False):
        if not append:
            self.nmatch = 0
            self.nindex = PatternPositionIndex()
            self.nnum = 0
        p = re.compile(self.pattern.sequence.replace('n', '[atcg]'),
                       re.IGNORECASE)
        prc = re.compile('%s' % (revcomp(self.pattern.sequence).replace('n', '[atcg]')),
                         re.IGNORECASE)

        for i in nset:
            self.nnum += 1
            has_match = False
            if self.reverse_complement:
                for j in prc.finditer(i):
                    has_match = True
                    self.nindex.append(self.nnum,
                                       self.pattern.sequence,
                                       i,
                                       j.start(),
                                       j.end() - 1,
                                       self.reverse_complement)
            else:
                for j in p.finditer(i):
                    self.nindex.append(self.nnum,
                                       self.pattern.sequence,
                                       i,
                                       j.start(),
                                       j.end() - 1,
                                       self.reverse_complement)

            if has_match:
                self.nmatch += 1


class PatternPositionIndex(object):

    def __init__(self):
        self.seqid = set()
        self.match = {}
        self.match_wildcard = {}
        self.match_nonwilcard = {}
        self.match_sequence = {}

    def append(self, seqid, query, hit, hit_start, hit_end, reverse_complement=False):
        if seqid not in self.seqid:
            self.seqid.add(seqid)
            self.match.update({seqid: []})
            self.match_wildcard.update({seqid: []})
            self.match_nonwildcard.update({seqid: []})
            self.match_sequence.update({seqid: []})
        self.match.get(seqid).append(range(hit_start, hit_end + 1))
        wildcard = []
        nonwildcard = []
        for i, j in enumerate(query):
            if j == 'n':
                wildcard.append(hit_start + i)
            else:
                nonwildcard.append(hit_start + i)
        self.match_wildcard.get(seqid).append(wildcard)
        self.match_nonwildcard.get(seqid).append(nonwildcard)
        if reverse_complement:
            self.match_sequence.get(seqid).append(revcomp(hit[hit_start: hit_end + 1]))
        else:
            self.match_sequence.get(seqid).append(hit[hit_start: hit_end + 1])


class PatternCollection(object):
    """A collection of unique pattern objects for
    pattern ranking"""

    def __init__(self):
        self._collect = {}

    def __iter__(self):
        return self._collect.itervalues()

    def __len__(self):
        return len(self._collect)

    def add(self, pattern, reverse_complement=False):
        """Add a pattern object"""
        assert isinstance(pattern, Pattern)
        if reverse_complement and revcomp(pattern.sequence) in self._collect:
            self._collect.update({revcomp(pattern.sequence): pattern})
        else:
            self._collect.update({pattern.sequence: pattern})

    def build(self, pset, nset, reverse_complement=False, append=False):
        for i in self:
            i.build(pset, nset, reverse_complement, append)

    # def pattern_scoring(self, sp_weight=1):
    #     for pattern in self:
    #         pattern.score = pattern.sd * (pattern.sp ** sp_weight)

    # def clustering(self, max_cluster=5, similarity=0.8):
    #     ranked_patterns = sorted(self, key=lambda x: x.score, reverse=True)
    #     clusters = [0] * len(self)
    #     clustered_patterns = {}
    #     num_cluster = 0

    #     for i in len(self):
    #         if num_cluster == max_cluster:
    #             break
    #         if clusters[i] == 0:
    #             num_cluster += 1
    #             clusters[i] = num_cluster
    #             clustered_patterns.update({num_cluster: [ranked_patterns[i]]})
    #             for j in range(i + 1, len(self)):
    #                 if position_weight_matrix(ranked_patterns[i], ranked_patterns[j]) >= similarity:
    #                     clusters[j] = num_cluster
    #                     clustered_patterns.get(i).append(ranked_patterns[j])
