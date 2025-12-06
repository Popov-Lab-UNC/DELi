"""code for degenerating barcode reads based on UMI tags"""

import abc
import dataclasses
from collections import defaultdict
from typing import Generic, TypeVar

from typing_extensions import Self

from .decoder import DecodedDELCompound, DecodedCompound, DecodedToolCompound
from .umi import UMI


# class DELCounter(abc.ABC):
#     """Base class for all decoded DEL counters"""
#
#     @abc.abstractmethod
#     def get_degen_count(self) -> int:
#         """Return the degen count"""
#         raise NotImplementedError()
#
#     @abc.abstractmethod
#     def get_raw_count(self) -> int:
#         """Return the raw count for the degenerator"""
#         raise NotImplementedError()
#
#     @abc.abstractmethod
#     def to_dict(self) -> dict:
#         """Convert the counter to a dictionary representation"""
#         raise NotImplementedError()
#
#
# class DELIdCounter(DELCounter):
#     """
#     Basic DEL ID counter; mainly for compatibility
#
#     Attributes
#     ----------
#     count: int
#         the current count for this DEL ID
#         will be both the degen and raw count
#     """
#
#     def __init__(self):
#         """Initialize the DEL ID counter"""
#         self.count: int = 0
#
#     def __add__(self, other):
#         """
#         Add two DELIdCounters together by summing the counts
#
#         Parameters
#         ----------
#         other: DELIdCounter
#             the other counter to add
#
#         Returns
#         -------
#         DELIdCounter
#         """
#         if isinstance(other, DELIdCounter):
#             new_counter = DELIdCounter()
#             new_counter.count = self.count + other.count
#             return new_counter
#         raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")
#
#     def add_id(self):
#         """Increment the ID counter by 1"""
#         self.count += 1
#
#     def get_degen_count(self) -> int:
#         """
#         Return the degen count
#
#         Notes
#         -----
#         For DELIdCounters, the degen count is the same as the
#         raw count.
#
#         Returns
#         -------
#         int
#             the degen count
#         """
#         return self.count
#
#     def get_raw_count(self) -> int:
#         """
#         Return the raw count for the degen counter
#
#         Returns
#         -------
#         int
#             the raw count
#         """
#         return self.count
#
#     def to_dict(self) -> dict:
#         """Covert the counter to a JSON compatible dictionary representation"""
#         return {"raw_count": self.count}
#
#
# class UMIDegenerate(abc.ABC):
#     """Base class for UMI Degenerators"""
#
#     @abc.abstractmethod
#     def add_umi(self, umi: UMI) -> bool:
#         """
#         Given a UMI, add it to the set and track its count
#
#         Will return a bool indicating if the UMI was added
#
#         Parameters
#         ----------
#         umi: UMI
#             the UMI to add
#
#         Returns
#         -------
#         bool
#             True if UMI is added (was novel)
#             False if it is not added (not novel)
#         """
#         raise NotImplementedError()
#
#
# class UMICounter(UMIDegenerate, DELCounter):
#     """
#     Counts the number of unique UMI occurrences
#
#     Attributes
#     ----------
#     umis: set[str]
#         the set of UMI already observed
#     """
#
#     def __init__(self):
#         """Initialize a UMI Counter"""
#         self.umis: set[str] = set()
#         self._raw_count = 0
#
#     def __add__(self, other):
#         """
#         Add two UMICounters together by merging the UMI set and summing raw counts
#
#         Parameters
#         ----------
#         other: UMICounter
#             the other counter to add
#
#         Returns
#         -------
#         UMICounter
#         """
#         if isinstance(other, UMICounter):
#             new_counter = UMICounter()
#             new_counter.umis = self.umis.union(other.umis)
#             new_counter._raw_count = self._raw_count + other._raw_count
#             return new_counter
#         raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")
#
#     def add_umi(self, umi: UMI) -> bool:
#         """
#         Given a UMI, add it to the set and track its count
#
#         Will return a bool indicating if the UMI was added
#
#         Parameters
#         ----------
#         umi: UMI
#             the UMI to add
#
#         Returns
#         -------
#         bool
#             True if UMI is added (was novel)
#             False if it is not added (not novel)
#         """
#         self._raw_count += 1
#         _found_match: bool = False
#         umi_ascii_code = "Null"  # for the IDE checker
#         for umi_ascii_code in umi.to_ascii_code():
#             if umi_ascii_code in self.umis:
#                 _found_match = True
#                 break
#         if not _found_match:
#             # add the last umi_ascii_code generated
#             self.umis.add(umi_ascii_code)
#             return True
#         return False
#
#     def get_degen_count(self) -> int:
#         """
#         Return the degen count
#
#         Returns
#         -------
#         int
#             the degen count
#         """
#         return len(self.umis)
#
#     def get_raw_count(self) -> int:
#         """
#         Return the raw count for the degen counter
#
#         Returns
#         -------
#         int
#             the raw count
#         """
#         return self._raw_count
#
#     def to_dict(self) -> dict:
#         """
#         Convert the counter to a dictionary representation
#
#         Returns
#         -------
#         dict
#             the dictionary representation of the counter
#         """
#         return {"raw_count": self._raw_count, "umis": list(self.umis)}


# TODO: build a better UMI clustering algorithm
# class UMICluster(UMIDegenerate, DELCounter):
#     """
#     Count UMI occurrences with clustering
#
#     Essentially, this will treat all UMIs
#     within some minimum distance to each other
#     as the same (in a greedy fashion).
#     Useful if you have a high error rate in sequencing
#     """
#
#     def __init__(self, min_dist: int = 2):
#         """
#         Initialize a UMICluster
#
#         Parameters
#         ----------
#         min_dist: int, default = 2
#             the minimum levenshtein distance between two UMIs
#             to be considered unique
#         """
#         self.min_dist = min_dist
#         self.umis: list[UMI] = list()
#         self._raw_count = 0
#
#     def __add__(self, other):
#         """
#         Add two UMIClusters together by merging the UMI clusters and summing raw counts
#
#         Parameters
#         ----------
#         other: UMICluster
#             the other counter to add
#
#         Returns
#         -------
#         UMICluster
#         """
#         if isinstance(other, UMICluster):
#             if self.min_dist != other.min_dist:
#                 raise ValueError(
#                     f"cannot add UMIClusters with different min_dist values: "
#                     f"{self.min_dist} and {other.min_dist}"
#                 )
#             new_counter = UMICluster(min_dist=self.min_dist)
#             _all_umis = self.umis + other.umis
#             for _umi in _all_umis:
#                 new_counter.add_umi(_umi)
#             new_counter._raw_count = self._raw_count + other._raw_count
#             return new_counter
#         raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")
#
#     def add_umi(self, umi: UMI) -> bool:
#         """
#         Given a UMI, add it to the set and track its count
#
#         Will return a bool indicating if the UMI was added
#
#         Parameters
#         ----------
#         umi: UMI
#             the UMI to add
#
#         Returns
#         -------
#         bool
#             True if UMI is added (was novel)
#             False if it is not added (not novel)
#         """
#         self._raw_count += 1
#         for existing_umi in self.umis:
#             # if too close, fail to add new UMI
#             if levenshtein_distance(existing_umi.umi_tag, umi.umi_tag) < self.min_dist:
#                 return False
#         self.umis.append(umi)
#         return True
#
#     def get_degen_count(self) -> int:
#         """
#         Return the degen count
#
#         Returns
#         -------
#         int
#             the degen count
#         """
#         return len(self.umis)
#
#     def get_raw_count(self) -> int:
#         """
#         Return the raw count for the degen counter
#
#         Returns
#         -------
#         int
#             the raw count
#         """
#         return self._raw_count
#
#     def to_dict(self) -> dict:
#         """
#         Convert the counter to a dictionary representation
#
#         Returns
#         -------
#         dict
#             the dictionary representation of the counter
#         """
#         return {"raw_count": self._raw_count, "umis": list(self.umis)}


# class DELIdUmiCounter(DELCounter):
#     """
#     Count the number of unique UMI reads for a single DEL ID
#
#     Supports both a UMI clustering mode and basic mode.
#     Clustering mode will consider any two UMIs within a given
#     minimum distance of each other as the same.
#     If you have high error rates in your sequencing
#     use the clustering mode. Otherwise, use the basic mode.
#
#     Notes
#     -----
#     UMI stands for "unique molecular identifier"
#     UMI degen is crucial for effective noise reduction
#     of the PCR amplification step in DEL as it removes
#     uneven amplification of a single read, since now
#     we know that each unique read has a unique UMI
#     """
#
#     def __init__(self):
#         """
#         Initialize DELIdUmiCounter
#         """
#         self.umis = UMICounter()
#
#     def __add__(self, other):
#         """
#         Add two DELIdUmiCounters together by merging the UMI clusters and summing raw counts
#
#         Parameters
#         ----------
#         other: DELIdUmiCounter
#             the other counter to add
#
#         Returns
#         -------
#         DELIdUmiCounter
#         """
#         if isinstance(other, DELIdUmiCounter):
#             new_counter = DELIdUmiCounter()
#             new_counter.umis = self.umis + other.umis
#             return new_counter
#         raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")
#
#     def add_umi(self, umi: UMI) -> bool:
#         """
#         Add a UMI to the umi counter
#
#         Will automatically handle the UMI already existing
#
#         Parameters
#         ----------
#         umi: UMI
#             the UMI to add
#
#         Returns
#         -------
#         bool
#             `True` if UMI is added has not been seen yet, else `False`
#         """
#         return self.umis.add_umi(umi)
#
#     def get_degen_count(self) -> int:
#         """
#         Return the degen count
#
#         Returns
#         -------
#         int
#             the degen count
#         """
#         return self.umis.get_degen_count()
#
#     def get_raw_count(self) -> int:
#         """
#         Return the raw count for the degen counter
#
#         Returns
#         -------
#         int
#             the raw count
#         """
#         return self.umis.get_raw_count()
#
#     def to_dict(self) -> dict:
#         """
#         Convert the counter to a dictionary representation
#
#         Returns
#         -------
#         dict
#             the dictionary representation of the counter
#         """
#         return self.umis.to_dict()
#
#
# # to help with type hinting
# DEL_COUNTER_TYPE = TypeVar("DEL_COUNTER_TYPE", bound=DELCounter)
#
#
# class DELCollectionCounter(abc.ABC, Generic[DEL_COUNTER_TYPE]):
#     """Base class for all DELibraryCollection degen counters"""
#
#     del_counter: dict[str, defaultdict[DecodedDELCompound, DEL_COUNTER_TYPE]]
#
#     def __len__(self):
#         """Get total number of unique DELs in the counter"""
#         return sum([len(barcodes) for barcodes in self.del_counter.values()])
#
#     @abc.abstractmethod
#     def __add__(self, other):
#         """
#         Add two DELCollectionCounter together by merging their counters
#
#         Parameters
#         ----------
#         other: DELCollectionCounter
#             the other counter to add
#
#         Returns
#         -------
#         DELCollectionCounter
#         """
#         raise NotImplementedError()
#
#     @abc.abstractmethod
#     def count_barcode(self, barcode: DecodedDELCompound) -> bool:
#         """
#         Given a barcode, add it to the current degen count
#
#         Will handle separating based on DEL ID and UMI
#         correction (if UMI is used)
#
#         Parameters
#         ----------
#         barcode: DecodedDELCompound
#             the barcode to add to the counter
#
#         Returns
#         -------
#         bool
#             `True` is added barcode is new (not degenerate), else `False`
#         """
#         raise NotImplementedError()
#
#     @abc.abstractmethod
#     def to_json(self, path: str | Path, compress: bool = False):
#         """
#         Convert the counter to a JSON serializable dict
#
#         Returns
#         -------
#         dict
#             the JSON serializable dict representation of the counter
#         """
#         raise NotImplementedError()
#
#
# class DELCollectionIdUmiCounter(DELCollectionCounter):
#     """
#     Degen counter for library collections that have UMI tags
#
#     Handles Degen by DEL ID and UMI.
#     Each observed DEL compound will have a count based
#     on its degen reads and a count based on how many
#     raw reads of the given compound there were
#     """
#
#     def __init__(self):
#         """
#         Initialize DELIdUmiCounter
#         """
#         # god forgive me for this one
#         self.del_counter: dict[str, defaultdict[DecodedDELCompound, DELIdUmiCounter]] = defaultdict(
#             lambda: defaultdict(partial(DELIdUmiCounter))
#         )
#
#     def __getstate__(self):
#         """
#         Covert state to generic dict with no partial functions when requesting state
#         """
#         state = self.__dict__.copy()
#         state["del_counter"] = {k: dict(v) for k, v in self.del_counter.items()}
#         return state
#
#     def __setstate__(self, state):
#         """
#         Set state for del_counter back to nested defaultdict when reloading state from getstate
#         """
#         self.__dict__.update(state)
#         self.del_counter = defaultdict(
#             lambda: defaultdict(partial(DELIdUmiCounter)),
#             {k: defaultdict(DELIdUmiCounter, v) for k, v in state["del_counter"].items()},
#         )
#
#     def __add__(self, other):
#         """
#         Add two DELCollectionIdUmiCounter objects together
#
#         Notes
#         -----
#         Adding two counters loaded from a pickle could cause
#         memory to grow unexpectedly
#
#         Parameters
#         ----------
#         other: DELCollectionIdUmiCounter
#             the other counter to add
#
#         Returns
#         -------
#         DELCollectionIdUmiCounter
#             the new counter with the sum of the two
#         """
#         if isinstance(other, DELCollectionIdUmiCounter):
#             new_counter = deepcopy(self)
#             for library_id in other.del_counter.keys():
#                 for compound in other.del_counter[library_id].keys():
#                     if compound in new_counter.del_counter[library_id]:
#                         new_counter.del_counter[library_id][compound] += other.del_counter[library_id][compound]
#                     else:
#                         new_counter.del_counter[library_id][compound] = deepcopy(
#                             other.del_counter[library_id][compound]
#                         )
#             return new_counter
#         raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")
#
#     def count_barcode(self, barcode: DecodedDELCompound) -> bool:
#         """
#         Given a barcode, add it to the current degen count
#
#         Will handle separating based on DEL ID and UMI
#         correction
#
#         Parameters
#         ----------
#         barcode: DecodedDELCompound
#             the barcode to add to the counter
#
#         Returns
#         -------
#         bool
#             `True` is added barcode is new (not degenerate), else `False`
#
#         Raises
#         ------
#         RuntimeError
#             if the decoded barcode is missing a UMI
#         """
#         if barcode.umi is None:
#             raise RuntimeError("cannot UMI degen on read missing UMI")
#         else:
#             return self.del_counter[barcode.library.library_id][barcode].add_umi(barcode.umi)
#
#     def to_json(self, path: str | Path, compress: bool = False):
#         """
#         Convert the counter to a JSON serializable dict
#
#         Parameters
#         ----------
#         path: str | Path
#             path to save the JSON file to
#         compress: bool, default = False
#             if True, will compress the JSON file using gzip
#             will be encoded with utf-8
#
#         Returns
#         -------
#         dict
#             the JSON serializable dict representation of the counter
#         """
#         import json
#
#         # collect data
#         _data: dict[str, dict[str, object]] = {}
#         for lib_ids, dels in self.del_counter.items():
#             _data[lib_ids] = {}
#             for compound, counter in dels.items():
#                 _id = compound.compound_id
#                 _info = {
#                     "lib_id": compound.library.library_id,
#                     "bb_ids": [bb.bb_id for bb in compound.building_blocks],
#                     "raw_count": counter.umis.get_raw_count(),
#                     "umis": list(counter.umis.umis),
#                 }
#                 _data[lib_ids][_id] = _info
#
#         if compress:
#             with gzip.open(path, "wt", encoding="utf-8") as zipfile:
#                 json.dump(_data, zipfile)
#         else:
#             json.dump(_data, open(path, "w"))
#
#
# class DELCollectionIdCounter(DELCollectionCounter):
#     """
#     Degen counter for library collectionss without UMI tags
#
#     Handles Degen by DEL ID
#     Each observed DEL compound will have a count based
#     on its degen reads and a count based on how many
#     raw reads of the given compound there were
#     """
#
#     def __init__(self):
#         self.del_counter: defaultdict[str, defaultdict[DecodedDELCompound, DELIdCounter]] = defaultdict(
#             lambda: defaultdict(DELIdCounter)
#         )
#
#     def __getstate__(self):
#         """
#         Covert state to generic dict with no partial functions when requesting state
#         """
#         state = self.__dict__.copy()
#         state["del_counter"] = {k: dict(v) for k, v in self.del_counter.items()}
#         return state
#
#     def __setstate__(self, state):
#         """
#         Set state for del_counter back to nested defaultdict when reloading state from getstate
#         """
#         self.__dict__.update(state)
#         self.del_counter = defaultdict(
#             lambda: defaultdict(DELIdCounter),
#             {k: defaultdict(DELIdUmiCounter, v) for k, v in state["del_counter"].items()},
#         )
#
#     def __add__(self, other):
#         """
#         Add two DELCollectionIdCounter together by merging the counters
#
#         Parameters
#         ----------
#         other: DELCollectionIdCounter
#             the other counter to add
#
#         Returns
#         -------
#         DELCollectionIdCounter
#         """
#         if isinstance(other, DELCollectionIdCounter):
#             new_counter = DELCollectionIdCounter()
#             for lib_id in self.del_counter.keys() | other.del_counter.keys():
#                 for barcode in self.del_counter[lib_id].keys() | other.del_counter[lib_id].keys():
#                     new_counter.del_counter[lib_id][barcode] = (
#                         self.del_counter[lib_id][barcode] + other.del_counter[lib_id][barcode]
#                     )
#             return new_counter
#         raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")
#
#     def count_barcode(self, barcode: DecodedDELCompound) -> bool:
#         """
#         Given a barcode, add it to the current degen count
#
#         Will handle separating based on DEL ID
#
#         Parameters
#         ----------
#         barcode: DecodedDELCompound
#             the barcode to add to the counter
#
#         Returns
#         -------
#         bool
#             `True` is added barcode is new (not degenerate), else `False`
#         """
#         self.del_counter[barcode.library.library_id][barcode].add_id()
#         return True  # always not degenerate with this counter
#
#     def to_json(self, path: str | Path, compress: bool = False):
#         """
#         Convert the counter to a JSON serializable dict
#
#         Parameters
#         ----------
#         path: str | Path
#             path to save the JSON file to
#         compress: bool, default = False
#             if True, will compress the JSON file using gzip
#
#         Returns
#         -------
#         dict
#             the JSON serializable dict representation of the counter
#         """
#         import json
#
#         # collect data
#         _data: dict[str, dict[str, object]] = {}
#         for lib_ids, dels in self.del_counter.items():
#             _data[lib_ids] = {}
#             for compound, counter in dels.items():
#                 _id = compound.compound_id
#                 _info = {
#                     "lib_id": compound.library.library_id,
#                     "bb_ids": [bb.bb_id for bb in compound.building_blocks],
#                     "raw_count": counter.get_raw_count(),
#                 }
#                 _data[lib_ids][_id] = _info
#
#         if compress:
#             with gzip.open(path, "wt", encoding="utf-8") as zipfile:
#                 json.dump(_data, zipfile)
#         else:
#             json.dump(_data, open(path, "w"))


class CompoundUMICounter:
    """
    Keeps track of degenerate compound count for a single compound

    Notes
    -----
    if the UMI is "null" or `None`, it will be counted
    as a *unique* compound, since there is no way to
    determine if it is a duplicate or not

    Attributes
    ----------
    counter: defaultdict[str, int]
        the counter for the compound
    """

    def __init__(self):
        self.counter: dict[str, int] = {"null": 0, "None": 0}

    def __add__(self, other):
        """add two CompoundUMICounters together by merging their counters"""
        if isinstance(other, CompoundUMICounter):
            new_counter = CompoundUMICounter()
            # merge counters
            for umi, count in self.counter.items():
                new_counter.counter[umi] = new_counter.counter.get(umi, 0) + count
            for umi, count in other.counter.items():
                new_counter.counter[umi] = new_counter.counter.get(umi, 0) + count
            return new_counter
        raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")

    @property
    def raw_count(self) -> int:
        """Get the raw count for the compound (number of times this compound was seen)"""
        return sum(lib_counter for lib_counter in self.counter.values())

    @property
    def degen_count(self) -> int:
        """Get the degenerate count for the compound (number of unique compounds seen)"""
        # exclude null and None from the unique count, but add their counts
        # since we are treating any compound missing a UMI as unique
        return len(self.counter) - 2 + self.counter.get("null", 0) + self.counter.get("None", 0)

    def add_umi(self, umi: str | UMI | None):
        """
        Add a UMI to the compound counter

        Parameters
        ----------
        umi: Any
            The UMI to add to the counter.
            Will cast to string

        Returns
        -------
        bool
            True if UMI was novel, else False
        """
        self.counter[str(umi)] += 1

C = TypeVar("C", bound=DecodedCompound)
class DegenCounter(abc.ABC, Generic[C]):
    """Base class for all degeneration counters"""

    @abc.abstractmethod
    def add_compound(self, compound: C):
        """
        Add a decoded compound to the degenerator

        Parameters
        ----------
        compound: DecodedCompound
            the decoded DEL compound to add
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def get_total_raw_count(self) -> int:
        """
        Get the total raw count

        Returns
        -------
        int
            the total raw count
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def get_total_degen_count(self) -> int:
        """
        Get the total degenerate count

        Returns
        -------
        int
            the total degenerate count
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def __add__(self, other: Self) -> Self:
        """
        Add two DegenCounters together by merging their counters

        Parameters
        ----------
        other: DegenCounter
            the other counter to add

        Returns
        -------
        DegenCounter
        """
        raise NotImplementedError()


class CompoundDegenCounter(DegenCounter[DecodedCompound]):
    """
    Counts degenerate compounds for any decoded compound

    Notes
    -----
    unlike `DELDegenCounter` and `ToolDegenCounter`, this will use the
    decoded compound directly as the key, rather than converting to a string
    representation. This is useful when you are degenerating in the same
    process as decoding, and you are not writing to a file after degeneration.
    This way, the information about the compound itself (building block SMILES, any
    enumerated SMILES, etc.) will be preserved. Otherwise, it needs to be reloaded
    later to recover that info.

    This is primarily used for the 'decode' CLI and python API when degeneration
    will occur in the same process as decoding.

    Attributes
    ----------
    counter: defaultdict[DecodedCompound, CompoundUMICounter]
        the counter for the compounds
    """
    def __init__(self):
        self.counter: defaultdict[DecodedCompound, CompoundUMICounter] = defaultdict(CompoundUMICounter)

    def __add__(self, other):
        """Add two CompoundDegenCounters together by merging their counters"""
        if isinstance(other, CompoundDegenCounter):
            new_counter = CompoundDegenCounter()
            # merge counters
            for compound, compound_counter in self.counter.items():
                new_counter.counter[compound] = new_counter.counter.get(compound, CompoundUMICounter()) + compound_counter
            for compound, compound_counter in other.counter.items():
                new_counter.counter[compound] = new_counter.counter.get(compound, CompoundUMICounter()) + compound_counter
            return new_counter
        raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")

    def add_compound(self, compound: DecodedCompound):
        self.counter[compound].add_umi(compound.umi)

    def get_total_raw_count(self) -> int:
        """Get the total raw count for the library"""
        return sum(compound_counter.raw_count for compound_counter in self.counter.values())

    def get_total_degen_count(self) -> int:
        """Get the total degenerate count for the library"""
        return sum(compound_counter.degen_count for compound_counter in self.counter.values())


class DELDegenCounter(DegenCounter[DecodedDELCompound]):
    """
    Counts degenerate DEL compounds for a single DEL

    Notes
    -----
    This assumes that all compounds added to the counter
    belong to the same DEL library. Therefore, the counter
    key is just the building block IDs.

    Attributes
    ----------
    counter: defaultdict[tuple[str, ...], CompoundUMICounter]
        the counter for the compounds
        the key is a tuple of building block IDs
    """
    def __init__(self):
        self.counter: defaultdict[tuple[str, ...], CompoundUMICounter] = defaultdict(CompoundUMICounter)

    def __add__(self, other):
        """Add two DELDegenCounter together by merging their counters"""
        if isinstance(other, DELDegenCounter):
            new_counter = DELDegenCounter()
            # merge counters
            for bb_id_tuple, compound_counter in self.counter.items():
                new_counter.counter[bb_id_tuple] = new_counter.counter.get(bb_id_tuple, CompoundUMICounter()) + compound_counter
            for bb_id_tuple, compound_counter in other.counter.items():
                new_counter.counter[bb_id_tuple] = new_counter.counter.get(bb_id_tuple, CompoundUMICounter()) + compound_counter
            return new_counter
        raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")

    def add_compound(self, compound: DecodedDELCompound | tuple[tuple[str, ...], str]):
        """
        Add a decoded DEL compound to the degenerator

        Will also accept the tuple of building block IDs directly
        attached to the UMI string

        Parameters
        ----------
        compound: DecodedDELCompound | tuple[tuple[str, ...], str]
            the decoded DEL compound to add
            can be the building block ID tuple directly with UMI
        """
        if isinstance(compound, tuple):
            self.counter[tuple(compound[0])].add_umi(compound[1])
        else:
            self.counter[tuple(compound.building_block_ids)].add_umi(compound.umi)

    def get_total_raw_count(self) -> int:
        """Get the total raw count for the library"""
        return sum(compound_counter.raw_count for compound_counter in self.counter.values())

    def get_total_degen_count(self) -> int:
        """Get the total degenerate count for the library"""
        return sum(compound_counter.degen_count for compound_counter in self.counter.values())


class ToolDegenCounter(DegenCounter[DecodedToolCompound]):
    """
    Counts degenerate tool compounds

    Notes
    -----
    Tool compounds are basically libraries with
    1 compound in them (that tool compounds ID)
    so degeneration is just counting unique
    occurrences of each tool compound UMI, rather
    than also tracking building block IDs

    Attributes
    ----------
    counter: CompoundUMICounter
        the counter for the tool compound
    """
    def __init__(self):
        self.counter = CompoundUMICounter()

    def __add__(self, other):
        """Add two ToolDegenCounters together by merging their counters"""
        if isinstance(other, ToolDegenCounter):
            new_counter = ToolDegenCounter()
            new_counter.counter = self.counter + other.counter
            return new_counter
        raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")

    def add_compound(self, compound: DecodedToolCompound | str):
        """
        Add a decoded tool compound to the degenerator

        Can also accept the UMI string directly
        (since there is only one possible compound per in the degen so id not needed)

        Parameters
        ----------
        compound: DecodedToolCompound | str
            the decoded tool compound to add
            or the UMI string directly
        """
        if isinstance(compound, str):
            self.counter.add_umi(compound)
        else:
            self.counter.add_umi(compound.umi)

    def get_total_raw_count(self) -> int:
        """Get the total raw count for the tool compound"""
        return self.counter.raw_count

    def get_total_degen_count(self) -> int:
        """Get the total degenerate count for the tool compound"""
        return self.counter.degen_count


@dataclasses.dataclass(frozen=True)
class DegenSettings:
    """
    Settings for running degeneration

    Degeneration is the process of collapsing reads based on UMI tags
    to get a more accurate count of unique molecules, avoiding noise
    created by PCR amplification and sequencing bias.

    Parameters
    ----------
    umi_clustering: bool, default = False
        if True, will use UMI clustering for UMI degeneration
    ignore_missing_umi: bool, default = False
        if True, will ignore reads missing UMI tags during degeneration;
        they will not be counted at all
    """

    umi_clustering: bool = False
    def to_file(self, path: str):
        """
        Save settings to a YAML file

        Parameters
        ----------
        path : str
            path to save settings to
        """
        import yaml
        yaml.dump(dataclasses.asdict(self), open(path, "w"))

    @classmethod
    def from_file(cls, path: str) -> "DegenSettings":
        """
        Load settings from a YAML file

        Will first check if there is a "decode_settings" key
        and load settings from that sub dict.
        Otherwise, will load from the YAML file keys

        Parameters
        ----------
        path : str
            Path to YAML file

        Returns
        -------
        DecodingSettings

        Raises
        ------
        RuntimeError
            if valid decode settings cannot be loaded from the passed YAML file
        """
        import yaml
        _data = yaml.safe_load(open(path, "r"))
        if "degen_settings" not in _data:
            try:
                return cls(**yaml.safe_load(open(path, "r")))
            except Exception as e:
                raise RuntimeError(f"Failed to load degen settings from {path}") from e
        else:
            try:
                return cls(**_data["degen_settings"])
            except Exception as e:
                raise RuntimeError(f"Failed to load degen settings from {path}") from e


class Degenerator(abc.ABC):
    """Base class for all degeneration runners"""

    @abc.abstractmethod
    def degen_decoded_compound(self, compound: DecodedCompound):
        """
        Degenerate a decoded compound

        Parameters
        ----------
        compound: DecodedCompound
            the decoded DEL compound to add
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def combine_degenerator(self, other: Self) -> Self:
        """
        Merge (add) another degenerator into this one

        Parameters
        ----------
        other: Degenerator
            the other degenerator to merge
        """
        raise NotImplementedError()


    def __add__(self, other: Self) -> Self:
        """Add two degenerators together by combining their counters"""
        return self.combine_degenerator(other)



class CompoundSelectionDegenerator(Degenerator):
    """
    Degenerate decoded barcodes based on UMI tags

    Unlike the SelectionDecoder, this is not embarrassingly parallel since degeneration
    requires tracking all observed UMIs for each DEL compound. This means it will not be able
    to operate on the fly and will need to keep everything in memory until all decoded
    compounds have been processed. If it stops before this, the degeneration results will
    have incorrect counts. Counts will never be lower than the actual degenerate counts, but
    they may be higher.

    This also means that if two Degenerator instances are run for the same
    selection (e.g. in parallel seeing different sets of decoded compounds), the results
    cannot just be simply be added together since the UMI sets will not be the same.
    `DELSelectionDegenerator` implement custom adding logic to do this merge correctly.
    They can also be saved and loaded from JSON files to enable merge separate processes.

    Notes
    -----
    This implementation of the degenerator will use the DecodedCompound objects
    directly as keys rather than converting to string representations. This is to
    preserve all information about the compounds during degeneration, rather than
    needing to reload them later. If you do not need this information (for example,
    if you are just writing to a file after degeneration), consider using

    Parameters
    ----------
    degen_settings: DegenSettings
        the settings for degeneration
        see `DegenSettings` for details

    Attributes
    ----------
    counter: defaultdict[str, CompoundDegenCounter]
        the counter for each library ID
    """
    def __init__(self, degen_settings: DegenSettings):
        self.settings = degen_settings
        self.counter: defaultdict[str, CompoundDegenCounter] = defaultdict(CompoundDegenCounter)

    def degen_decoded_compound(self, compound: DecodedCompound):
        """
        Degenerate a decoded compound

        Parameters
        ----------
        compound: DecodedCompound
            the decoded DEL compound to add
        """
        library_id = compound.get_library_id()
        self.counter[library_id].add_compound(compound)

    def combine_degenerator(self, other: Self) -> Self:
        """
        Merge (add) another degenerator into this one

        Parameters
        ----------
        other: Degenerator
            the other degenerator to merge

        Returns
        -------
        CompoundSelectionDegenerator
            the merged degenerator
        """
        new_degenerator = CompoundSelectionDegenerator(self.settings)
        new_degenerator.merge_degenerator(self)
        new_degenerator.merge_degenerator(other)
        return new_degenerator

    def merge_degenerator(self, other: Self):
        """
        Merge (add) another degenerator into this one in place

        Parameters
        ----------
        other: Degenerator
            the other degenerator to merge
        """
        for lib_id, counter in other.counter.items():
            self.counter[lib_id] = self.counter[lib_id] + counter


class SelectionDegenerator(Degenerator):
    """
    Degenerate decoded barcodes based on UMI tags

    Unlike the SelectionDecoder, this is not embarrassingly parallel since degeneration
    requires tracking all observed UMIs for each DEL compound. This means it will not be able
    to operate on the fly and will need to keep everything in memory until all decoded
    compounds have been processed. If it stops before this, the degeneration results will
    have incorrect counts. Counts will never be lower than the actual degenerate counts, but
    they may be higher.

    This also means that if two SelectionDegenerator instances are run for the same
    selection (e.g. in parallel seeing different sets of decoded compounds), the results
    cannot just be simply be added together since the UMI sets will not be the same.
    `DELSelectionDegenerator` implement custom adding logic to do this merge correctly.
    They can also be saved and loaded from JSON files to enable merge separate processes.

    Degenerators will group compounds by their library ID automatically.
    This will avoid excessive memory usage encoding the same library multiple times.

    The SelectionDegenerator is also keyed using the tuple of building block id
    strings rather than the full DEL compound ID. This is to avoid issues where
    extracting the building block ids and compound ideas from string representations
    can be error-prone. Also, the compound ID hold the library ID, which is already
    being used to group the compounds in the degenerator making this key redundant.

    Notes
    -----
    This implementation of the degenerator uses string representations of the
    decoded compounds as keys rather than the DecodedCompound objects directly.
    This means it can degenerate from decoded compound files. It will also
    convert any decoded compound it degenerates into a string representation.
    This means that some information about the compounds (building block SMILES,
    any enumerated SMILES, etc.) will be lost during degeneration unless the
    compounds are reloaded later. If you want to preserve this information during
    degeneration, consider using `CompoundSelectionDegenerator`.

    Parameters
    ----------
    degen_settings: DegenSettings
        the settings for degeneration
        see `DegenSettings` for details

    Attributes
    ----------
    counter: dict[str, DELDegenCounter | ToolDegenCounter]
        the counter for each library ID
        will be either a DELDegenCounter or ToolDegenCounter
        based on the compound type
    """
    def __init__(self, degen_settings: DegenSettings):
        self.settings = degen_settings
        self.counter: dict[str, DELDegenCounter | ToolDegenCounter] = {}

    def degen_decoded_compound(self, compound: DecodedDELCompound | DecodedToolCompound | tuple[tuple[str, ...], str]):
        # handle decoded compound object input
        if not isinstance(compound, tuple):
            lib_id = compound.get_library_id()
            if lib_id not in self.counter:
                if isinstance(compound, DecodedDELCompound):
                    self.counter[lib_id] = DELDegenCounter()
                elif isinstance(compound, DecodedToolCompound):
                    self.counter[lib_id] = ToolDegenCounter()
                else:
                    raise RuntimeError(f"This error should be unreachable, please report a bug!")
            self.counter[lib_id].add_compound(compound)
        # handle string based input
        else:
            compound_id, umi = compound
            if len(compound_id) == 1:  # assume tool compound
                if compound_id[0] not in self.counter:
                    self.counter[compound_id[0]] = ToolDegenCounter()
                self.counter[compound_id[0]].add_compound(umi)
            else:
                lib_id = compound_id[0]
                if lib_id not in self.counter:
                    self.counter[lib_id] = DELDegenCounter()
                self.counter[lib_id].add_compound((compound[1:], umi))

    def combine_degenerator(self, other: Self) -> Self:
        """
        Merge (add) another degenerator into this one

        Parameters
        ----------
        other: Degenerator
            the other degenerator to merge

        Returns
        -------
        SelectionDegenerator
            the merged degenerator
        """
        new_degenerator = SelectionDegenerator(self.settings)
        new_degenerator.merge_degenerator(self)
        new_degenerator.merge_degenerator(other)
        return new_degenerator

    def merge_degenerator(self, other: Self):
        """
        Merge (add) another degenerator into this one in place

        Parameters
        ----------
        other: Degenerator
            the other degenerator to merge
        """
        for lib_id, counter in other.counter.items():
            if lib_id in self.counter:
                self.counter[lib_id] = self.counter[lib_id] + counter
            else:
                self.counter[lib_id] = counter
