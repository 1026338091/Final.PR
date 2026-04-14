"""
sequence_utils.py
Provides utility functions for biological sequence manipulation, 
specifically handling bi-directional strand orientation.
"""

def get_reverse_complement(seq):
    """
    Returns the reverse complement of a given DNA sequence.
    """
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    rev_comp = "".join([complement_map.get(base.upper(), 'N') for base in reversed(seq)])
    return rev_comp

def extract_kmers_both_strands(sequence, k):
    """
    Extracts all possible k-mers from a sequence, including their forward 
    and reverse-complement orientations.
    """
    kmers = []
    n = len(sequence)
    for i in range(n - k + 1):
        forward_kmer = sequence[i:i+k]
        reverse_kmer = get_reverse_complement(forward_kmer)
        kmers.append((forward_kmer, i, '+'))
        kmers.append((reverse_kmer, i, '-'))
    return kmers
