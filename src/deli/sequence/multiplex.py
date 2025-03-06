"""handle demultiplexing reads"""

from typing import Union

from deli.dels import DELibraryPool
from deli.utils import get_min_levenshtein_distance

from .align import Aligner, HybridSemiGlobalAligner, HybridSemiGlobalAlignment, SemiGlobalAligner


class DemultiplexingError(Exception):
    """raised when an issue with demultiplexing occurs"""

    pass


class DemultiplexAlignment(HybridSemiGlobalAlignment):
    """
    Hold info on a multiplexed sequence alignment

    Attributes
    ----------
    adapt_overlap: str
        the region of the query sequence that is aligned adapter sequence
    """

    def __init__(
        self, seq1: str, seq2: str, alignment: list[tuple[int, int]], score: int, edit_cost: int
    ):
        """
        Initialize the alignment object

        Parameters
        ----------
        seq1: str
            the top sequence in the alignment (index 0)
        seq2: str
            the bottom sequence in the alignment (index 1)
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
        super().__init__(seq1, seq2, alignment, score, edit_cost)

        self.adapt_overlap: str = seq1[slice(self.map_seq2_span_to_seq1_span(0, len(seq2) - 1))]


class Demuliplexer:
    """
    Generic demultiplexer class

    Used to set up a set of possible adapters to demultiplex[1]_ on.

    Note: should not be used for DELibrary demultiplexing.
    Use the LibraryDemultiplexer instead.

    References
    ----------
    .. [1] https://edna.dnabarcoding101.org/bioinformatics.html
    """

    def __init__(self, tags: Union[list[str], dict[str, str]], align_mode: str = "semi"):
        self.tags: dict[str, str]
        if isinstance(tags, list):
            self.tags = {str(i): tag for i, tag in enumerate(tags)}
        else:
            self.tags = tags

        self.align_mode: Aligner
        if align_mode == "semi":
            self.align_mode = SemiGlobalAligner()
        elif align_mode == "hybrid":
            self.align_mode = HybridSemiGlobalAligner()
        else:
            raise ValueError(f"unrecognized demultiplex align_mode '{align_mode}'")

        ### VALIDATE ###
        if len(set(self.tags.values())) != len(self.tags):
            for tag_id1, tag1 in self.tags.items():
                for tag_id2, tag2 in enumerate(self.tags.items()):
                    if (tag_id1 != tag_id2) and (tag1 == tag2):
                        raise DemultiplexingError(
                            f"tags '{tag_id1}' and '{tag_id2}' have identical"
                            f"DNA tag: '{tag_id1}'"
                        )

        self.min_levenshtein_distance = get_min_levenshtein_distance(list(self.tags.values()))


class LibraryDemultiplexer(Demuliplexer):
    """Used to demultiplex a DEL pool into specific single DEL sets"""

    def __init__(self, library_pool: DELibraryPool):
        """
        Initialize the library demultiplexer

        Parameters
        ----------
        library_pool: DELibraryPool
            the pool of possible libraries
        """
        super().__init__(
            tags={lib.library_id: lib.library_tag for lib in library_pool}, align_mode="hybrid"
        )
