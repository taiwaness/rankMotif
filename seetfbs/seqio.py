def revcomp(sequence):
    """Convert a DNA sequence into its reverse-complement
    counterpart"""
    table = {
        'a': 't',
        't': 'a',
        'c': 'g',
        'g': 'c',
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'n': 'n',
    }

    rcseq = []
    for i in sequence[::-1]:
        rcseq.append(table.get(i))

    return ''.join(rcseq)
