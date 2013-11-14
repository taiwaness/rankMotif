import re
from seqio import revcomp


class Pattern(object):

    def __init__(self, sequence):
        self.sequence = sequence
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

    class MatchPosition(object):

        def __init__(self):
            self.seqid = set()
            self.pos_matches = {}
            self.pos_wildcards = {}
            self.pos_nonwilcards = {}
            self.match_sequences = {}

        def add(self, seqid, query, hit, hit_start, hit_end, reverse_match=False):
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
            if reverse_match:
                self.match_sequences.get(seqid).append(revcomp(hit[hit_start: hit_end]))
            else:
                self.match_sequences.get(seqid).append(hit[hit_start: hit_end])

    def __init__(self, reverse_complement=False):
        self.reverse_complement = reverse_complement
        self._match_position = self.MatchPosition()
        self.n_hitseqs = 0
        self.n_hitsites = 0
        self.n_seqs = 0

    def index(self, sequence, seqset, append=False):
        if not append:
            self._match_position = self.MatchPosition()
            self.n_hitseqs = 0
            self.n_hitsites = 0
            self.n_seqs = 0

        repl = [
            ('n', '[atcg]'),
            ('p', '[at]'),
            ('q', '[ac]'),
            ('r', '[ag]'),
            ('s', '[tc]'),
            ('u', '[tg]'),
            ('v', '[cg]'),
        ]
        seq_retype = sequence
        rcseq_retype = revcomp(sequence)
        for i in repl:
            seq_retype = seq_retype.replace(*i)
            rcseq_retype = rcseq_retype.replace(*i)

        p = re.compile('(%s)|(%s)' % (seq_retype, rcseq_retype), re.IGNORECASE)

        for i in seqset:
            self.n_seqs += 1
            has_match = False
            for j in p.finditer(i):
                has_match = True
                self.n_hitsites += 1
                if self.reverse_complement:
                    if j.group(1):
                        self.index.add(self.n_seqset, sequence, i, j.start(), j.end(), False)
                    elif j.group(2):
                        self.index.add(self.n_seqset, sequence, i, j.start(), j.end(), True)
                elif j.group(1):
                    self.index.add(self.n_seqset, sequence, i, j.start(), j.end(), False)
            if has_match:
                self.n_hitseqs += 1

    @property
    def seqid(self):
        return self._match_position.seqid

    @property
    def pos_wildcards(self):
        return self._match_position.pos_wildcards

    @property
    def pos_nonwilcards(self):
        return self._match_position.pos_nonwilcards

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


def merge_pattern(seq_1, seq_2):
    """Merge several pattarns and return the matching sequences"""

    p = {'a', 't'}
    q = {'a', 'c'}
    r = {'a', 'g'}
    s = {'t', 'c'}
    u = {'t', 'g'}
    v = {'c', 'g'}

    if len(seq_1) > len(seq_2):
        reference = seq_1
        scanner = seq_2
    elif len(seq_1) < len(seq_2):
        reference = seq_2
        scanner = seq_1
    else:
        merged_seq = []
        for i in range(len(seq_1)):
            if seq_1[i] == seq_2[i]:
                merged_seq.append(seq_1[i])
            elif seq_1[i] in p and seq_2[i] in p:
                merged_seq.append('p')
            elif seq_1[i] in q and seq_2[i] in q:
                merged_seq.append('q')
            elif seq_1[i] in r and seq_2[i] in r:
                merged_seq.append('r')
            elif seq_1[i] in s and seq_2[i] in s:
                merged_seq.append('s')
            elif seq_1[i] in u and seq_2[i] in u:
                merged_seq.append('u')
            elif seq_1[i] in v and seq_2[i] in v:
                merged_seq.append('v')
            else:
                merged_seq.append('n')
        return ''.join(merged_seq)

    best_score = -1
    best_match = None
    for i in range(len(reference) - len(scanner) + 1):
        score = alignment(scanner, reference[i: i + len(scanner)])
        if score > best_score:
            best_score = score
            best_match = ('%s%s%s' % ('n' * i, scanner, 'n' * (len(reference) - len(scanner) - i)), reference)

    return merge_pattern(*best_match)


def alignment(seq_1, seq_2):
    assert len(seq_1) == len(seq_2)

    score = 0
    for i in range(len(seq_1)):
        if seq_1[i] == seq_2[i]:
            score += 1
