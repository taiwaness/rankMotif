def revcomp(sequence):
    """Convert a DNA sequence into its reverse-complement
    counterpart"""
    table = {
        'a': 't',
        't': 'a',
        'c': 'g',
        'g': 'c',
        'n': 'n',
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'N': 'N',
    }

    rcseq = []
    for i in sequence[::-1]:
        if i in table:
            rcseq.append(table.get(i))
        else:
            raise Exception('[revcomp] Unsupported base: {0}'.format(i))

    return ''.join(rcseq)


def parse_fasta(handle):
    """handle: path to FASTA file or file object"""
    is_path = False
    if isinstance(handle, str):
        is_path = True
        fi = open(handle, 'r')
    else:
        fi = handle

    header = ''
    sequence = []
    while True:
        line = fi.readline()
        if line == '':
            if header and sequence:
                yield (header, ''.join(sequence))
            break

        line = line.strip()
        if line == '':
            continue
        elif line[0] == '>':
            if header and sequence:
                yield (header, ''.join(sequence))
                header = line[1:]
                sequence = []
            else:
                header = line[1:]
        else:
            sequence.append(line)

    if is_path:
        fi.close()


def parse_fasta_noheader(handle):
    """handle: path to FASTA file or file object"""
    is_path = False
    if isinstance(handle, str):
        is_path = True
        fi = open(handle, 'r')
    else:
        fi = handle

    header = ''
    sequence = []
    while True:
        line = fi.readline()
        if line == '':
            if header and sequence:
                yield ''.join(sequence)
            break

        line = line.strip()
        if line == '':
            continue
        elif line[0] == '>':
            if header and sequence:
                yield ''.join(sequence)
                header = line[1:]
                sequence = []
            else:
                header = line[1:]
        else:
            sequence.append(line)

    if is_path:
        fi.close()
