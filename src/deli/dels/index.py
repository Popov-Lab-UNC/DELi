"""Code for handling indexes"""

import json
from typing import Iterator, List, Optional, Self, Union

from Levenshtein import distance

from deli.configure import accept_deli_data

from .base import DeliDataLoadableMixin


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


class IndexSet(DeliDataLoadableMixin):
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

    def __getitem__(self, index: int) -> Index:
        """Get the Index at the passed index from the set"""
        return self.index_set[index]

    @classmethod
    @accept_deli_data("indexes", "json")
    def load(cls, path: str) -> Self:
        """
        Load a index set from the DELi data directory

        Notes
        -----
        This is decorated by `accept_deli_data`
        which makes this function actually take
          path_or_name: str
          deli_config: DeliConfig

        `path_or_name` can be the full path to the file
        or it can be the name of the object to load

        See `Storing DEL info` in docs for more details


        Parameters
        ----------
        path: str
            path of the index set to load

        Returns
        -------
        IndexSet
        """
        return cls.from_json(path)

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
        return cls(index_set=[Index(index_id=key, dna_tag=val) for key, val in data.items()])


def get_min_index_distance(included_indexes: Optional[Union[list[Index], IndexSet]]) -> int:
    """
    Determine the minimum Levenshtein distance between all Indexes

    Notes
    -----
    Can use a IndexSet or a list of Index objects

    Parameters
    ----------
    included_indexes: list[Index] or IndexSet
        the index_ids to use

    Returns
    -------
    min_distance: int
        the minimum Levenshtein distance between the passed indexes
    """
    # if only one index just give it a constant threshold
    if included_indexes is None or len(included_indexes) <= 1:
        return 0
    _index_sequences = [i.dna_tag for i in included_indexes]
    return min(
        [distance(s1, s2) for s1 in _index_sequences for s2 in _index_sequences if s1 != s2]
    )
