# !/usr/bin/env python3


def convert_cigar(cigar):
    """
    Converts the CIGAR string to a list of tuples
    Args:
        cigar: CIGAR string from pysam record

    Returns:
        list of tuples
    """
    cigar_char = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8, 'B': 9}
    cigartuples = []
    length = ''
    for c in cigar:
        if c.isdigit(): length += c
        else:
            cigartuples.append((cigar_char[c], int(length)))
            length = ''

    return cigartuples


def convert_cigartuples(cigartuples):
    """
    Converts the CIGAR tuples to a string
    Args:
        cigartuples: CIGAR tuples from pysam record

    Returns:
        CIGAR string
    """
    cigar_char = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 6: 'P', 7: '=', 8: 'X', 9: 'B'}
    cigar = ''
    for c in cigartuples:
        cigar += f'{c[1]}{cigar_char[c[0]]}'

    return cigar


def remove_softclips(cigar):
    """
    Removes the softclips from the CIGAR string

    Args:
        cigar(str): CIGAR string

    Returns:
        CIGAR string without softclips
    """

    new_cigar = ''
    length = ''
    for c in cigar:
        if c.isdigit(): length += c
        else:
            if c != 'S': new_cigar += f'{length}{c}'
            length = ''
    return new_cigar


def complement_cigar(cigar):
    """
    Reverses the CIGAR string
    Args:
        cigar: CIGAR string from pysam record

    Returns:
        list of tuples
    """
    new_cigar = ''
    cigar = remove_softclips(cigar)
    for c in cigar:
        if c == 'I': c = 'D'
        elif c == 'D': c = 'I'
        new_cigar += c
    return new_cigar
