"""code for demultiplexing in DELi"""

from typing import Optional

from cutadapt.adapters import MultipleAdapters, SingleMatch
from dnaio import SequenceRecord


class Demultiplexer:
    """
    Base class for all demultiplexers

    This is a very thin wrapper to cutadapt to enable easier maintainability.
    To understand the logarithm used to do the demultiplexing,
    read the cutadapt documentation
    """

    def __init__(self, adapters: MultipleAdapters):
        """
        Initialize the Demultiplexer

        Parameters
        ----------
        adapters: MultipleAdapters
            the adapters to use fr demultiplexing
        """
        self.adapters = adapters

    def _demultiplex(self, sequence: SequenceRecord) -> SingleMatch:
        """Given a observed_seq, match it to the adapters"""
        return self.adapters.match_to(sequence.sequence)

    @staticmethod
    def _compare_matches(
        match1: tuple[SequenceRecord, Optional[SingleMatch]],
        match2: tuple[SequenceRecord, Optional[SingleMatch]],
    ) -> tuple[SequenceRecord, Optional[SingleMatch]]:
        """
        Compare two matches (seq and its match) to see which one is better

        If either match is None, will return the other.
        If both are None, will return match1 (still with the None match).
        If score is tied, will return match1.
        Otherwise will return match with best (highest) score.

        Parameters
        ----------
        match1: tuple[SequenceRecord, Optional[SingleMatch]]
            the first observed_seq and its match
            match can be `None` if no match is found
        match2: tuple[SequenceRecord, Optional[SingleMatch]]
            the second observed_seq and its match
            match can be `None` if no match is found

        Returns
        -------
        tuple[SequenceRecord, Optional[SingleMatch]]
            the best scoring match and its accompanying observed_seq
        """
        if match1[1] is None:
            return match2
        elif match2[1] is None:
            return match1
        elif match1[1].score >= match2[1].score:
            return match1
        else:
            return match2
