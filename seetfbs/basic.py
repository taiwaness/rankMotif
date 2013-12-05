import re
import logging
from .seqio import revcomp, parse_fasta_noheader


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
            self.pos_matches.get(seqid).append(list(xrange(hit_start, hit_end)))
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
                        self._match_position.add(
                            self.n_seqs, sequence, i, j.start(), j.end(), False)
                    elif j.group(2):
                        self._match_position.add(
                            self.n_seqs, rc_sequence, i, j.start(), j.end(), True)
                elif j.group(1):
                    self._match_position.add(
                        self.n_seqs, sequence, i, j.start(), j.end(), False)
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


class MergePattern(object):

    def __init__(self, reverse_complement):
        self._logger = logging.getLogger(self.__class__.__name__)
        self.reverse_complement = reverse_complement
        self.patterns = []
        self._strands = {}

    def add(self, pattern, strand):
        self.patterns.append(pattern)
        self._strands.update({pattern: strand})

    def extract_match_sequences(self, fasta_handle):
        pos_matches = {}
        match_sequences = []

        for pattern in self._strands:
            pattern.build_matchtable_pset(
                parse_fasta_noheader(fasta_handle), self.reverse_complement)

        for pattern in self._strands:
            for seqid, pos in pattern.matchtable_pset.pos_matches.iteritems():
                for i, j in enumerate(pos):
                    if seqid in pos_matches and j in pos_matches.get(seqid):
                        continue
                    match_strand, match_sequence = pattern.matchtable_pset.match_sequences.get(seqid)[i]
                    if match_strand == 2:
                        match_sequence = revcomp(match_sequence)
                    if seqid in pos_matches:
                        pos_matches.get(seqid).append(j)
                    else:
                        pos_matches.update({seqid: [j]})

                    match_sequences.append(match_sequence)

        return match_sequences


def merge_patterns(pattern_list, reverse_complement=False):
    """Merge patterns and return directional pattern objects
    and their relative strand directions"""
    merged_patterns = MergePattern(reverse_complement)
    if len(pattern_list) > 1:
        reference = pattern_list[0].sequence
        for pattern in pattern_list[1:]:
            reference = merge_sequences(reference, pattern.sequence, reverse_complement)[0]
        merged_patterns.add(Pattern(reference), 1)
        for pattern in pattern_list[1:]:
            seq_1, seq_2, strands = merge_sequences(reference, pattern.sequence, reverse_complement)
            merged_patterns.add(Pattern(seq_2), strands[1])
    else:
        merged_patterns.add(Pattern(pattern_list[0].sequence), 1)

    return merged_patterns


def merge_sequences(seq_1, seq_2, reverse_complement=False):
    arrlen = len(seq_1) + len(seq_2) - 1
    max_score = -1
    best_seq1 = None
    best_seq2 = None
    best_strands = None

    strands_1 = [1, 1]
    for i in xrange(arrlen - len(seq_1) + 1):
        for j in xrange(arrlen - len(seq_2) + 1):
            if j > i:
                left_1 = ''
                left_2 = 'n' * (j - i)
            elif i > j:
                left_1 = 'n' * (i - j)
                left_2 = ''
            else:
                left_1 = ''
                left_2 = ''

            if len(seq_2) + j - 1 > len(seq_1) + i - 1:
                right_1 = 'n' * ((len(seq_2) + j - 1) - (len(seq_1) + i - 1))
                right_2 = ''
            elif len(seq_1) + i - 1 > len(seq_2) + j - 1:
                right_1 = ''
                right_2 = 'n' * ((len(seq_1) + i - 1) - (len(seq_2) + j - 1))
            else:
                right_1 = ''
                right_2 = ''

            s1 = '%s%s%s' % (left_1, seq_1, right_1)
            s2 = '%s%s%s' % (left_2, seq_2, right_2)
            score = full_alignment_scoring(s1, s2)
            if score > max_score:
                max_score = score
                best_seq1 = s1
                best_seq2 = s2
                best_strands = strands_1

    if reverse_complement:
        strands_2 = [1, 2]
        seq_2 = revcomp(seq_2)

        for i in xrange(arrlen - len(seq_1) + 1):
            for j in xrange(arrlen - len(seq_2) + 1):
                if j > i:
                    left_1 = ''
                    left_2 = 'n' * (j - i)
                elif i > j:
                    left_1 = 'n' * (i - j)
                    left_2 = ''
                else:
                    left_1 = ''
                    left_2 = ''

                if len(seq_2) + j - 1 > len(seq_1) + i - 1:
                    right_1 = 'n' * ((len(seq_2) + j - 1) - (len(seq_1) + i - 1))
                    right_2 = ''
                elif len(seq_1) + i - 1 > len(seq_2) + j - 1:
                    right_1 = ''
                    right_2 = 'n' * ((len(seq_1) + i - 1) - (len(seq_2) + j - 1))
                else:
                    right_1 = ''
                    right_2 = ''

                s1 = '%s%s%s' % (left_1, seq_1, right_1)
                s2 = '%s%s%s' % (left_2, seq_2, right_2)
                score = full_alignment_scoring(s1, s2)
                if score > max_score:
                    max_score = score
                    best_seq1 = s1
                    best_seq2 = s2
                    best_strands = strands_2

    return (best_seq1, best_seq2, best_strands)


def full_alignment_scoring(seq_1, seq_2):
    assert len(seq_1) == len(seq_2), (seq_1, seq_2)

    score = 0
    for i in xrange(len(seq_1)):
        if seq_1[i] == seq_2[i] != 'n':
            score += 1

    return score
