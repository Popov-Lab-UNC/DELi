"""DNA sequence alignment"""

from typing import List, Optional, Self, Tuple

from numba import njit


class SemiGlobalAlignment:
    """
    Contains the global alignment of two sequences.

    Parameters
    ----------
    seq1: str
        first sequence from alignment
    seq2: str
        second sequence from alignment
    alignment: List[Tuple[Optional[int], Optional[int]]]
        the alignment of the two squences

    Notes
    -----
    This class should never be manually instantiated, it should be
    generated using the class method `globally_align`
    """

    def __init__(
        self,
        seq1: str,
        seq2: str,
        alignment: List[Tuple[Optional[int], Optional[int]]],
    ):
        self.seq1 = seq1
        self.seq2 = seq2
        self.alignment = alignment

    def __iter__(self):
        """Iterate over the alignment indexes"""
        for idx_1, idx_2 in self.alignment:
            yield idx_1, idx_2

    @classmethod
    def globally_align(cls, seq1: str, seq2: str) -> Self:
        """
        Globally align two sequences to each other

        Parameters
        ----------
        seq1: str
            sequence to be aligned
        seq2: str
            sequence to align to

        Returns
        -------
        alignment: list[(int, int)]
            the resulting global alignment of the two sequences

        Notes
        -----
        Uses the Needleman–Wunsch algorithm[1]_ to calculate the global alignment

        The alignment will be returned as a list of tuples.
        The tuples contain two elements: the index of the first sequence and
        the index of the second sequence that it was aligned to.
        If there is gap, one of the two index values will be `None`.

        For example, alignment:

        ACTG--C
        -CTGCCC

        would look like this:
        [(0, None), (1, 0), (2, 1), (3, 2), (None, 3), (None, 4), (4, 5)]

        References
        ----------
        .. [1] Needleman, S. B.; Wunsch, C. D. A General Method Applicable to the
        Search for Similarities in the Amino Acid Sequence of Two Proteins.
        Journal of Molecular Biology 1970, 48 (3), 443–453.

        Examples
        --------
        >>> seq_1 = 'ACTG'
        >>> seq_2 = 'ATGC'
        >>> my_alignment = SemiGlobalAlignment.globally_align(seq_1, seq_2)
        [(0, 0), (1, None), (2, 1), (3, 2), (None, 3)]
        """
        return cls(seq1, seq2, _semi_global_align(seq1, seq2))

    def __str__(self) -> str:
        """Render the alignment as a a readable formatted string"""
        return (
            "".join("-" if i is None else self.seq1[i] for i, _ in self.alignment)
            + "\n"
            + "".join("-" if j is None else self.seq2[j] for _, j in self.alignment)
        )


@njit()
def _semi_global_align(seq1: str, seq2: str) -> List[Tuple[Optional[int], Optional[int]]]:
    """
    Numba accelerated semi-global alignment implementation

    Should use SemiGlobalAlignment.globally_align() to align two sequences
    this function should never be called outside of this file
    """
    n, m = len(seq1), len(seq2)

    table = {}
    direction_table = {}

    table[-1, -1] = 0
    for i in range(n):
        table[i, -1] = 0
    for j in range(m):
        table[-1, j] = 0

    for i in range(n):
        for j in range(m):
            choice = (
                table[i - 1, j - 1] + (seq1[i] == seq2[j]),
                table[i - 1, j] - 1,
                table[i, j - 1] - 1,
            )
            table[i, j], direction_table[i, j] = max(zip(choice, ((-1, -1), (-1, 0), (0, -1))))

    alignment: List[Tuple[Optional[int], Optional[int]]] = []
    element: Tuple[Optional[int], Optional[int]]

    # Find the maximum score at the last row or column (semi-global alignment)
    last_row_scores = [table[n - 1, j] for j in range(m)]
    last_col_scores = [table[i, m - 1] for i in range(n)]

    # Determine if we start traceback from the last row or last column
    if max(last_row_scores) >= max(last_col_scores):
        i, j = n - 1, last_row_scores.index(max(last_row_scores))
    else:
        i, j = last_col_scores.index(max(last_col_scores)), m - 1

    while i >= 0 and j >= 0:
        direction = direction_table[i, j]
        if direction == (-1, -1):
            element = i, j
        elif direction == (-1, 0):
            element = i, None
        else:
            element = None, j
        alignment.insert(0, element)
        di, dj = direction
        i, j = i + di, j + dj
    while i >= 0:
        alignment.insert(0, (i, None))
        i -= 1
    while j >= 0:
        alignment.insert(0, (None, j))
        j -= 1

    return alignment
