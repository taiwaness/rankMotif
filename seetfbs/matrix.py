from basic import Pattern, MatchTable


def position_frequency_matrix(pattern):
    """Claculate position frequency matrix (PFM)
    of matching sequences of a pattern"""
    assert isinstance(pattern, Pattern)
    assert isinstance(pattern.matchtable, MatchTable)
    sequences = [i.lower() for i in pattern.matchtable.pindex.match_sequence.itervalues()]
    lrow = max([len(i) for i in sequences])

    matrix = {
        'a': [0] * lrow,
        't': [0] * lrow,
        'c': [0] * lrow,
        'g': [0] * lrow,
    }

    for s in sequences:
        for i, j in enumerate(s):
            matrix[j][i] += 1

    return matrix


def position_weight_matrix(pfm1, pfm2):
    pass
