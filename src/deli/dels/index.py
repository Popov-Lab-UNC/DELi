"""Code for handling indexes"""

import json
from typing import Iterator, List, Optional, Self, Union

from Levenshtein import distance

from deli.constants import MAX_INDEX_RISK_DIST_THRESHOLD


class Index:
    """
    Object that contains the index id and corresponding DNA tag
    """

    def __init__(self, index_id: str, dna_tag: str):
        """
        Initializes the object with the index id and corresponding DNA tag

        Parameters
        ----------
        index_id: str
            id of the index
        dna_tag: str
            the DNA tag of the index
        """
        self.index_id = index_id
        self.dna_tag = dna_tag


class IndexSet:
    """
    Holds a set of indexes

    Notes
    -----
    Useful for multiplexed run where index decoding is needed during calling
    """

    def __init__(self, index_set: List[Index]):
        """
        initializes the object with the index set

        Parameters
        ----------
        index_set: List[Index]
            the index objects that make up the index set
        """
        self.index_set = index_set

    def __len__(self) -> int:
        """Return the number of indexes in the IndexSet"""
        return len(self.index_set)

    def __iter__(self) -> Iterator[Index]:
        """Iterate all indexes in the IndexSet"""
        return iter(self.index_set)

    @classmethod
    def from_json(cls, path: str) -> Self:
        """
        Load a index set from a JSON file

        Notes
        -----
        See the "De-multiplexing with DELi" docs for more info

        Parameters
        ----------
        path: str
            path to json string

        Returns
        -------
        IndexSet

        """
        data = json.load(open(path))
        return cls(index_set=[Index(**d) for d in data])


def get_min_index_distance(included_index: Optional[Union[list[Index], IndexSet]]) -> int:
    """
    Determine the minimum Levenshtein distance between all Indexes

    Notes
    -----
    Can use a IndexSet or a list of Index objects

    Parameters
    ----------
    included_index: list[Index] or IndexSet
        the index_ids to use

    Returns
    -------
    min_distance: int
        the minimum Levenshtein distance between the passed indexes
    """
    # if only one index just give it a constant threshold
    if len(included_index) == 1:
        return MAX_INDEX_RISK_DIST_THRESHOLD
    if included_index is None:
        return 0
    _index_sequences = [i.dna_tag for i in included_index]
    return min(
        [distance(s1, s2) for s1 in _index_sequences for s2 in _index_sequences if s1 != s2]
    )
