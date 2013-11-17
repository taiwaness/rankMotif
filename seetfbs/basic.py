import re
from .seqio import revcomp


class Pattern(object):

    def __init__(self, sequence):
        self.sequence = sequence.lower()
        self.n_wildcards = sequence.count('n')
        self.n_nonwildcards = len(sequence) - self.n_wildcards
        self.matchtable_pset = None
        self.matchtable_nset = None

    def __len__(self):
        return len(self.sequence)

    def build_matchtable_pset(self, seqset, reverse_complement=False, append=False):
        if not self.matchtable_pset or not append:
            self.matchtable_pset = MatchTable(reverse_complement)
        self.matchtable_pset.index(self.sequence, seqset, append)

        return self

    def build_matchtable_nset(self, seqset, reverse_complement=False, append=False):
        if not self.matchtable_nset or not append:
            self.matchtable_nset = MatchTable(reverse_complement)
        self.matchtable_nset.index(self.sequence, seqset, append)

        return self


class MatchTable(object):
    """Index of a pattern sequence against the matching sequences"""

    class _MatchPosition(object):

        def __init__(self):
            self.seqid = set()
            self.pos_matches = {}
            self.pos_wildcards = {}
            self.pos_nonwildcards = {}
            self.match_sequences = {}

        def add(self, seqid, query, hit, hit_start, hit_end, is_rc_match):
            if seqid not in self.seqid:
                self.seqid.add(seqid)
                self.pos_matches.update({seqid: []})
                self.pos_wildcards.update({seqid: []})
                self.pos_nonwildcards.update({seqid: []})
                self.match_sequences.update({seqid: []})
            self.pos_matches.get(seqid).append(range(hit_start, hit_end))
            pos_wildcards = []
            pos_nonwildcards = []
            for i, j in enumerate(query):
                if j == 'n':
                    pos_wildcards.append(hit_start + i)
                else:
                    pos_nonwildcards.append(hit_start + i)
            self.pos_wildcards.get(seqid).append(pos_wildcards)
            self.pos_nonwildcards.get(seqid).append(pos_nonwildcards)
            if is_rc_match:
                mseq = (2, hit[hit_start: hit_end])
                self.match_sequences.get(seqid).append(mseq)
            else:
                mseq = (1, hit[hit_start: hit_end])
                self.match_sequences.get(seqid).append(mseq)

    def __init__(self, reverse_complement=False):
        self.reverse_complement = reverse_complement
        self._match_position = self._MatchPosition()
        self.n_hitseqs = 0
        self.n_hitsites = 0
        self.n_seqs = 0

    def index(self, sequence, seqset, append=False):
        if not append:
            self._match_position = self._MatchPosition()
            self.n_hitseqs = 0
            self.n_hitsites = 0
            self.n_seqs = 0

        rc_sequence = revcomp(sequence)
        p = re.compile('({0})|({1})'.format(
            sequence.replace('n', '[atcg]'), rc_sequence.replace('n', '[atcg]')),
            re.IGNORECASE)

        for i in seqset:
            i = i.lower()
            self.n_seqs += 1
            has_match = False
            for j in p.finditer(i):
                has_match = True
                self.n_hitsites += 1
                if self.reverse_complement:
                    if j.group(1):
                        self._match_position.add(self.n_seqs, sequence, i, j.start(), j.end(), False)
                    elif j.group(2):
                        self._match_position.add(self.n_seqs, rc_sequence, i, j.start(), j.end(), True)
                elif j.group(1):
                    self._match_position.add(self.n_seqs, sequence, i, j.start(), j.end(), False)
            if has_match:
                self.n_hitseqs += 1

    @property
    def seqid(self):
        return self._match_position.seqid

    @property
    def pos_matches(self):
        return self._match_position.pos_matches

    @property
    def pos_wildcards(self):
        return self._match_position.pos_wildcards

    @property
    def pos_nonwildcards(self):
        return self._match_position.pos_nonwildcards

    @property
    def match_sequences(self):
        return self._match_position.match_sequences


class PatternSet(object):
    """A set of pattern objects without redundant pattern sequences"""

    def __init__(self, reverse_complement=False):
        self.reverse_complement = reverse_complement
        self._collect = {}

    def __iter__(self):
        return self._collect.itervalues()

    def __len__(self):
        return len(self._collect)

    def add(self, pattern):
        assert isinstance(pattern, Pattern)
        if self.reverse_complement and revcomp(pattern.sequence) in self._collect:
            self._collect.update({revcomp(pattern.sequence): pattern})
        else:
            self._collect.update({pattern.sequence: pattern})

    def iterseqs(self):
        return self._collect.iterkeys()

    def remove(self, sequence):
        self._collect.pop(sequence)


def merge_pattern(pattern_1, pattern_2):
    seq_1 = pattern_1.sequence
    seq_2 = pattern_2.sequence
    arrlen = len(seq_1) + len(seq_2) - 1

    max_score = -1
    best_seq1 = None
    best_seq2 = None
    best_strands = None

    strands_1 = [1, 1]
    for i in range(arrlen):
        s1 = '%s%s%s' % ('n' * i, seq_1, 'n' * (arrlen - i - len(seq_1) + 1))
        s2 = '%s%s%s' % ('n' * (arrlen - i - len(seq_2) + 1), seq_2, 'n' * i)
        score = full_alignment_scoring(s1, s2)
        if score > max_score:
            max_score = score
            best_seq1 = s1
            best_seq2 = s2
            best_strands = strands_1

    strands_2 = [2, 1]
    rc_seq_1 = revcomp(seq_1)
    for i in range(arrlen):
        s1 = '%s%s%s' % ('n' * i, rc_seq_1, 'n' * (arrlen - i - len(rc_seq_1) + 1))
        s2 = '%s%s%s' % ('n' * (arrlen - i - len(seq_2) + 1), seq_2, 'n' * i)
        score = full_alignment_scoring(s1, s2)
        if score > max_score:
            max_score = score
            best_seq1 = seq_1
            best_seq2 = seq_2
            best_strands = strands_2

    return ((best_strands[0], Pattern(best_seq1)), (best_strands[1], Pattern(best_seq2)))


def full_alignment_scoring(seq_1, seq_2):
    assert len(seq_1) == len(seq_2)

    score = 0
    for i in range(len(seq_1)):
        if seq_1[i] == seq_2[i]:
            score += 1
