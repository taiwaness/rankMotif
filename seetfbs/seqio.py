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
    }

    rcseq = []
    for i in sequence[::-1]:
        if i in table:
            rcseq.append(table.get(i))
        else:
            rcseq.append(i)

    return ''.join(rcseq)
