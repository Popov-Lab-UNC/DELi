"""DNA observed_seq alignment"""

import abc
from typing import no_type_check

from numba import njit


class Alignment:
    """
    Holds information about an alignment between two sequences
    """

    def __init__(self, seq1: str, seq2: str, alignment: list[tuple[int, int]]):
        """
        Initialize the alignment object

        Parameters
        ----------
        seq1: str
            the top observed_seq in the alignment (index 0)
        seq2: str
            the bottom observed_seq in the alignment (index 1)
        alignment: list[tuple[int, int]]
            the alignment as a list of tuples, each containing the index of from
            the two aligned sequences that match as that step.
            `-1` is used to represent a gap
            for example:
                seq1: ACTG--C
                seq2: -CCGCCC
            would look like this:
            [(0, -1), (1, 0), (2, 1), (3, 2), (-1, 3), (-1, 4), (4, 5)]
        """
        self.seq1: str = seq1
        self.seq2: str = seq2
        self.alignment: list[tuple[int, int]] = alignment

        self._alignment_lookup_seq1, self._alignment_lookup_seq2 = _get_alignment_lookups(
            self.alignment
        )

    def print_alignment(self) -> str:
        """
        Prints out the alignment in a human readable clustal-like format.

        The format is the following:
        "|": match
        "X": mismatch
        " "/"-": gap

        For example:
        ------CCAGT--CCTTTCCTGAGAGT
              |||||  |||X|
        GCTTGCCCAGTGGCCTCT---------

        Returns
        -------
        str
        """
        _row1 = "".join("-" if i == -1 else self.seq1[i] for i, _ in self.alignment)
        _row3 = "".join("-" if j == -1 else self.seq2[j] for _, j in self.alignment)

        _row2 = ""
        for i, j in self.alignment:
            if (i == -1) or (j == -1):
                _row2 += " "
            elif i == j:
                _row2 += "|"
            else:
                _row2 += "X"

        return f"{_row1}\n{_row2}\n{_row3}"

    def map_seq1_span_to_seq2_span(self, start_idx: int, stop_idx: int) -> tuple[int, int]:
        """
        Returns the span from seq2 that is aligned to the passed seq1 span

        Notes
        -----
        uses numerical indexes of the sequences

        Parameters
        ----------
        start_idx: int
            the start of the span in seq1
        stop_idx: int
            the end of the span in seq1

        Returns
        -------
        tuple[int, int]
            the span from seq2 aligned to the passed seq1 span
        """
        return self._alignment_lookup_seq1[start_idx], self._alignment_lookup_seq1[stop_idx]

    def map_seq2_span_to_seq1_span(self, start_idx: int, stop_idx: int) -> tuple[int, int]:
        """
        Returns the span from seq1 that is aligned to the passed seq2 span

        Notes
        -----
        uses numerical indexes of the sequences

        Parameters
        ----------
        start_idx: int
            the start of the span in seq2
        stop_idx: int
            the end of the span in seq2

        Returns
        -------
        tuple[int, int]
            the span from seq1 aligned to the passed seq2 span
        """
        return self._alignment_lookup_seq2[start_idx], self._alignment_lookup_seq2[stop_idx]

    def __iter__(self):
        """Iterate over the alignment indexes"""
        for idx_1, idx_2 in self.alignment:
            yield idx_1, idx_2

    def __repr__(self) -> str:
        """Represent the alignment in clustal-like format (see `print_alignment`)"""
        return self.print_alignment()


class SemiGlobalAlignment(Alignment):
    """Hold results of a SemiGlobal Alignment"""

    def __init__(self, seq1: str, seq2: str, alignment: list[tuple[int, int]], score: int):
        """
        Initialize the alignment object

        Parameters
        ----------
        seq1: str
            the top observed_seq in the alignment (index 0)
        seq2: str
            the bottom observed_seq in the alignment (index 1)
        alignment: list[tuple[int, int]]
            the alignment as a list of tuples, each containing the index of from
            the two aligned sequences that match as that step.
            `-1` is used to represent a gap
            for example:
                seq1: ACTG--C
                seq2: -CCGCCC
            would look like this:
            [(0, -1), (1, 0), (2, 1), (3, 2), (-1, 3), (-1, 4), (4, 5)]
        score: int
            the score of the alignment
        """
        super().__init__(seq1, seq2, alignment)

        self.score = score


class HybridSemiGlobalAlignment(SemiGlobalAlignment):
    """Hold results of a Hybrid SemiGlobal Alignment"""

    def __init__(
        self,
        seq1: str,
        seq2: str,
        alignment: list[tuple[int, int]],
        score: int,
        edit_cost: int,
    ):
        """
        Initialize the alignment object

        Parameters
        ----------
        seq1: str
            the top observed_seq in the alignment (index 0)
        seq2: str
            the bottom observed_seq in the alignment (index 1)
        alignment: list[tuple[int, int]]
            the alignment as a list of tuples, each containing the index of from
            the two aligned sequences that match as that step.
            `-1` is used to represent a gap
            for example:
                seq1: ACTG--C
                seq2: -CCGCCC
            would look like this:
            [(0, -1), (1, 0), (2, 1), (3, 2), (-1, 3), (-1, 4), (4, 5)]
        score: int
            the score of the alignment
        edit_cost: int
            the edit cost of the alignment
        """
        super().__init__(seq1, seq2, alignment, score)

        self.edit_cost = edit_cost


class Aligner(abc.ABC):
    """
    Base class for all objects that can generate alignments

    On initialization, will pre-compile the alignment code with
    numba to accelerate runtime.
    This might cause obj initialization time to take 2-3 seconds.
    """

    def __init__(self):
        """Initialize the aligner"""
        # dummy align to force compile
        self.align("CCC", "CCA")

    @abc.abstractmethod
    def align(self, seq1: str, seq2: str) -> Alignment:
        """
        Align two sequences

        Parameters
        ----------
        seq1: str
            the top observed_seq in the alignment (index 0)
        seq2: str
            the bottom observed_seq in the alignment (index 1)
        """
        raise NotImplementedError


class SemiGlobalAligner(Aligner):
    """Aligner to align two seqs using a semi-global alignment"""

    def align(self, seq1: str, seq2: str) -> SemiGlobalAlignment:
        """
        Semi-globally align two sequences

        In a semi-global alignment, gaps at the start of end of both sequences
        are not penalized.

        Notes
        -----
        Uses the Needleman–Wunsch algorithm[1]_ to calculate the semi global alignment.
        Gaps and mismatches are both penalized with -1.

        Parameters
        ----------
        seq1: str
            observed_seq to be aligned
        seq2: str
            observed_seq to align to

        Returns
        -------
        alignment: SemiGlobalAlignment
            the resulting alignment of the two sequences

        References
        ----------
        .. [1] Needleman, S. B.; Wunsch, C. D. A General Method Applicable to the
        Search for Similarities in the Amino Acid Sequence of Two Proteins.
        Journal of Molecular Biology 1970, 48 (3), 443–453.

        Examples
        --------
        >>> seq_1 = 'ACTG'
        >>> seq_2 = 'ATCC'
        >>> aligner = SemiGlobalAligner()
        >>> my_alignment = aligner.align(seq_1, seq_2)
        >>> print(my_alignment)
        ACTG-
        | |X
        A-TCC
        """
        _alignment, _top_score = _semi_global_align(seq1, seq2)
        return SemiGlobalAlignment(seq1, seq2, alignment=_alignment, score=_top_score)


class HybridSemiGlobalAligner(Aligner):
    """Contains the hybrid semi global alignment of two sequences."""

    def align(self, seq1: str, seq2: str) -> HybridSemiGlobalAlignment:
        """
        Align two sequences using a hybrid semi-global alignment

        Notes
        -----
        This is inspired by the cutadapt[1]_ approach for alignment.
        In this algorithm, the dynamic programming matrix is filling is
        guided by the "edit cost" which is the raw number of mismatches
        in the alignment. The directions to mover are also created from this
        table. At the same time, the a socre matrix is filled
        based on the status of the edit cost alignment.
        When picking an alignment, the score table is naviagted using the edit
        table directions.

        To avoid the aligner decided that a full gap between the two sequences
        is optimal, only the first observed_seq has not penalties for gaps.
        The second observed_seq will have all gaps penalized.

        This makes the Hyprid approach far better at aligning a small
        region to a much larger observed_seq. For example, for trying to find
        the adapter in region in a larger read (whicih is why cutadpat uses it)

        For the above reason, this should really only be used when demultiplexing
        in DELi, for example the LibraryCaller uses it. All other alignments
        should be done with the SemiGlobalAligner.

        Parameters
        ----------
        seq1: str
            observed_seq to be aligned
        seq2: str
            observed_seq to align to

        Returns
        -------
        alignment: HybridSemiGlobalAlignment
            the resulting global alignment of the two sequences

        References
        ----------
        .. [1] https://cutadapt.readthedocs.io/en/stable/algorithms.html
        """
        _alignment, _top_score, _edit_cost = _hybrid_semi_global_align(seq1, seq2)
        return HybridSemiGlobalAlignment(seq1, seq2, _alignment, _top_score, _edit_cost)


@no_type_check
# @njit()
def _hybrid_semi_global_align(
    ref_seq: str, adapt_seq: str
) -> tuple[list[tuple[int, int]], int, int]:
    """
    Numba accelerated hybrid semi-global alignment implementation

    This is heavily drawn from cutadapt.
    It is not a full semi-global alignment,
    it will only ignore penalties for gaps in the top observed_seq.
    This is because if it ignored both, the algorithm
    would always prefer a complete misalignment:
    -----GCGC
    GTCGA----

    This approach is better for alignments required for
    demultiplexing, as it prioritizes better matches
    earlier on in the observed_seq to the adapter tag

    Note: this should not be called directly
    """
    n, m = len(ref_seq), len(adapt_seq)

    edit_table = {}
    score_table = {}
    direction_table = {}

    edit_table[-1, -1] = 0
    score_table[-1, -1] = 0
    # turn off gap free penalty to ref_seq
    for i in range(n):
        edit_table[i, -1] = 0
        score_table[i, -1] = 0
    for j in range(m):
        edit_table[-1, j] = j
        score_table[-1, j] = -j

    for i in range(n):
        for j in range(m):
            base_pair_match = ref_seq[i] == adapt_seq[j]

            choice = (
                edit_table[i - 1, j - 1] + (not base_pair_match),
                edit_table[i - 1, j] + 1,
                edit_table[i, j - 1] + (1 if (j != (m - 1)) else 0),
            )

            if base_pair_match:
                scores = (
                    score_table[i - 1, j - 1] + 1,
                    score_table[i - 1, j] - 2,
                    score_table[i, j - 1] - 2,
                )
            else:
                scores = (
                    score_table[i - 1, j - 1] - 1,
                    score_table[i - 1, j] - 2,
                    score_table[i, j - 1] - 2,
                )

            edit_table[i, j], direction_table[i, j], score_table[i, j] = min(
                zip(choice, ((-1, -1), (-1, 0), (0, -1)), scores)
            )

    alignment = []

    # Find the maximum score at the last column (semi-global alignment just seq 2)
    last_col_scores = [score_table[i, m - 1] for i in range(n)]
    top_score = max(last_col_scores)

    i, j = last_col_scores.index(top_score), m - 1
    edit_cost = edit_table[i, j]

    # this bit of code will also add tailing Gap matches
    # not sure if that is really needed and would slow
    # down the alignment. for not we won't use it

    while i < n - 1:
        alignment.insert(0, (n - 1, -1))
        n -= 1
    while j < m - 1:
        alignment.insert(0, (-1, m - 1))
        m -= 1

    while i >= 0 and j >= 0:
        direction = direction_table[i, j]
        if direction == (-1, -1):
            element = i, j
        elif direction == (-1, 0):
            element = i, -1
        else:
            element = -1, j
        alignment.insert(0, element)
        di, dj = direction
        i, j = i + di, j + dj
    while i >= 0:
        alignment.insert(0, (i, -1))
        i -= 1
    while j >= 0:
        alignment.insert(0, (-1, j))
        j -= 1

    return alignment, top_score, edit_cost


@no_type_check
@njit()
def _semi_global_align(seq1: str, seq2: str) -> tuple[list[tuple[int, int]], int]:
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

    alignment: list[tuple[int, int]] = []
    element: tuple[int, int]

    # Find the maximum score at the last row or column (semi-global alignment)
    last_row_scores = [table[n - 1, j] for j in range(m)]
    last_col_scores = [table[i, m - 1] for i in range(n)]

    # Determine if we start traceback from the last row or last column
    max_last_row_score = max(last_row_scores)
    max_last_col_score = max(last_col_scores)

    top_score: int
    if max_last_row_score >= max_last_col_score:
        i, j = n - 1, last_row_scores.index(max(last_row_scores))
        top_score = max(last_row_scores)
    else:
        i, j = last_col_scores.index(max(last_col_scores)), m - 1
        top_score = max(last_col_scores)

    # this bit of code will also add tailing Gap matches
    # not sure if that is really needed and would slow
    # down the alignment. for not we won't use it

    while i < n - 1:
        alignment.insert(0, (n - 1, -1))
        n -= 1
    while j < m - 1:
        alignment.insert(0, (-1, m - 1))
        m -= 1

    while i >= 0 and j >= 0:
        direction = direction_table[i, j]
        if direction == (-1, -1):
            element = i, j
        elif direction == (-1, 0):
            element = i, -1
        else:
            element = -1, j
        alignment.insert(0, element)
        di, dj = direction
        i, j = i + di, j + dj
    while i >= 0:
        alignment.insert(0, (i, -1))
        i -= 1
    while j >= 0:
        alignment.insert(0, (-1, j))
        j -= 1

    return alignment, top_score


def _get_alignment_lookups(
    alignment: list[tuple[int, int]],
) -> tuple[dict[int, int], dict[int, int]]:
    """
    Generate a mapping between the index of seq1 to seq2 index and vice versa

    in the case of gaps, will use the most recent real index of the other observed_seq

    For example, the alignment:
    AG-CCA--
    AGGCCATG

    would generate the maps
    {0:0, 1:1, 2:3, 3:4, 4:5} for seq1->seq2
    {0:0, 1:1, 2:1, 3:2, 4:3, 5:4, 6:4, 7:4} for seq1->seq2

    This is useful for when you can to take the span of one
    seq and map it to the span it matched in the other
    """
    alignment_lookup_seq1 = {}
    _last_non_null_seq1 = 0
    alignment_lookup_seq2 = {}
    _last_non_null_seq2 = 0

    max_seq1_idx = 0
    max_seq2_idx = 0

    for seq1_idx, seq2_idx in alignment:
        _last_non_null_seq1 = max(seq1_idx, _last_non_null_seq1)
        _last_non_null_seq2 = max(seq2_idx, _last_non_null_seq2)
        if seq1_idx != -1:
            alignment_lookup_seq1[seq1_idx] = _last_non_null_seq2
            max_seq1_idx = seq1_idx
        if seq2_idx != -1:
            alignment_lookup_seq2[seq2_idx] = _last_non_null_seq1
            max_seq2_idx = seq2_idx

    # need to add a ghost lookup at the end cause indexing from 0
    alignment_lookup_seq1[max_seq1_idx + 1] = _last_non_null_seq2 + 1
    alignment_lookup_seq2[max_seq2_idx + 1] = _last_non_null_seq1 + 1

    return alignment_lookup_seq1, alignment_lookup_seq2
