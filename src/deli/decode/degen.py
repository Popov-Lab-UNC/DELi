"""Functions for handling degenerate barcode counts based on UMI sequences"""

from numba import njit


@njit
def _get_dna_hamming_neighbors(seq) -> list[str]:
    """
    Given a DNA sequence, return all sequences with Hamming distance of 1

    Parameters
    ----------
    seq: str

    Returns
    -------
    list[str]
    """
    seqs = []
    for i, let in enumerate(seq):
        for new in ["A", "G", "C", "T"]:
            if new == let:
                continue
            seqs.append(seq[:i] + let + seq[i + 1 :])
    return seqs
