"""Code for handling indexes"""

import os.path
from os import PathLike
from typing import Iterator, List, Optional, Self, Union, overload

from Levenshtein import distance

from deli.configure import DeliDataLoadable, accept_deli_data_name


class DELIndexError(Exception):
    """exception for issues when load indexes"""

    pass


class Index(DeliDataLoadable):
    """
    Object that contains the index id and corresponding DNA bases
    """

    def __init__(self, index_id: str, dna_tag: str, sample_name: Optional[str] = None):
        """
        Initializes the object with the index id and corresponding DNA bases

        Parameters
        ----------
        index_id: str
            id of the index
        dna_tag: str
            the DNA bases of the index
        sample_name: Optional[str], default = None
            name of a the sample linked to the experiment
        """
        self.index_id = index_id
        self.dna_tag = dna_tag
        self.sample_name = sample_name

    def __len__(self):
        """Gets the length of the index DNA bases"""
        return len(self.dna_tag)

    def __eq__(self, other):
        """Two indexes are equal if they have the same id and DNA bases"""
        if isinstance(other, Index):
            return (self.index_id == other.index_id) and (self.dna_tag == other.dna_tag)

    @classmethod
    @accept_deli_data_name("indexes", "txt")
    def load(cls, path: Union[str, PathLike]) -> Self:
        """Load an Index from a txt file"""
        index_id = os.path.basename(path).split(".")[0]
        index_tag = open(path, "r").readlines()[0].strip()
        _tmp = cls(index_id=index_id, dna_tag=index_tag)
        _tmp.loaded_from = path
        return _tmp


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
        self._check_validity()

    def __len__(self) -> int:
        """Return the number of indexes in the IndexSet"""
        return len(self.index_set)

    def __iter__(self) -> Iterator[Index]:
        """Iterate all indexes in the IndexSet"""
        return iter(self.index_set)

    @overload
    def __getitem__(self, index: int) -> Index: ...

    @overload
    def __getitem__(self, index: slice) -> List[Index]: ...

    def __getitem__(self, index: Union[int, slice]) -> Union[Index, List[Index]]:
        """Get the Index at the passed index from the set"""
        return self.index_set[index]

    def __add__(self, other) -> Self:
        """Merge two index sets together"""
        if not isinstance(other, self.__class__):
            raise TypeError(f"cannot add type {type(other)} to {self.__class__}")

        return self.__class__(self.index_set + other.index_set)

    def __bool__(self) -> bool:
        """IndexSet is False when it is empty"""
        return len(self) > 0

    # @classmethod
    # @accept_deli_data_name("indexes", "json")
    # def load(cls, path: str) -> Self:
    #     """
    #     Load a index set from the DELi data directory
    #
    #     Notes
    #     -----
    #     This is decorated by `accept_deli_data`
    #     which makes this function actually take
    #       path_or_name: str
    #       deli_config: DeliConfig
    #
    #     `path_or_name` can be the full path to the file
    #     or it can be the name of the object to load
    #
    #     See `Storing DEL info` in docs for more details
    #
    #
    #     Parameters
    #     ----------
    #     path: str
    #         path of the index set to load
    #
    #     Returns
    #     -------
    #     IndexSet
    #     """
    #     _cls = cls.from_json(path)
    #     _cls.loaded_from = path
    #     return _cls
    #
    # @classmethod
    # def from_json(cls, path: str) -> Self:
    #     """
    #     Load a index set from a JSON file
    #
    #     Notes
    #     -----
    #     See the "De-multiplexing with DELi" docs for more info
    #
    #     Parameters
    #     ----------
    #     path: str
    #         path to json string
    #
    #     Returns
    #     -------
    #     IndexSet
    #
    #     """
    #     data = json.load(open(path))
    #     return cls(index_set=[Index(index_id=key, dna_tag=val) for key, val in data.items()])

    def _check_validity(self):
        """Checks that there are no duplicate or conflicts in index set"""
        _ids = []
        _tags = []
        for index in self.index_set:
            # check id uniqueness
            if index.index_id in _ids:
                if index.dna_tag == self.index_set[_ids.index(index.index_id)]:
                    raise DELIndexError("identical indexes found in index set")
                else:
                    raise DELIndexError(
                        "multiple indexes have the same id; index_ids must be unique"
                    )
            else:
                _ids.append(index.index_id)

            # check the bases uniqueness
            if index.dna_tag in _tags:
                _idx = _tags.index(index.dna_tag)
                raise DELIndexError(
                    f"index {index.index_id} and index {self.index_set[_idx]} "
                    f"have the same dna bases: {index.dna_tag}"
                )
            else:
                _tags.append(index.index_id)

    def has_index_with_name(self, index_name: str) -> bool:
        """
        Returns True if the Set has a index with the passed ID/name

        Parameters
        ----------
        index_name: str
            name/id of the index

        Returns
        -------
        bool
        """
        return any([idx.index_id == index_name for idx in self.index_set])

    def get_index_with_name(self, index_name: str) -> Index:
        """
        Returns Index with given name from the set

        Parameters
        ----------
        index_name: str
            name/id of the index

        Raises
        ------
        KeyError
            if the passed index name is not in the IndexSet

        Returns
        -------
        Index
        """
        for index in self.index_set:
            if index.index_id == index_name:
                return index
        raise KeyError(f"cannot find index with name '{index_name}' in IndexSet")


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
