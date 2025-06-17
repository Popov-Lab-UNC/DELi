"""code for degenerating barcode reads"""

import abc
import gzip
from collections import defaultdict
from copy import deepcopy
from functools import partial
from pathlib import Path
from typing import Generic, TypeVar

from Levenshtein import distance as levenshtein_distance

from .decoder import DecodedDELCompound
from .umi import UMI


class DELCounter(abc.ABC):
    """Base class for all decoded DEL counters"""

    @abc.abstractmethod
    def get_degen_count(self) -> int:
        """Return the degen count"""
        raise NotImplementedError()

    @abc.abstractmethod
    def get_raw_count(self) -> int:
        """Return the raw count for the degenerator"""
        raise NotImplementedError()

    @abc.abstractmethod
    def to_dict(self) -> dict:
        """Convert the counter to a dictionary representation"""
        raise NotImplementedError()


class DELIdCounter(DELCounter):
    """
    Basic DEL ID counter; mainly for compatibility

    Attributes
    ----------
    count: int
        the current count for this DEL ID
        will be both the degen and raw count
    """

    def __init__(self):
        """Initialize the DEL ID counter"""
        self.count: int = 0

    def __add__(self, other):
        """
        Add two DELIdCounters together by summing the counts

        Parameters
        ----------
        other: DELIdCounter
            the other counter to add

        Returns
        -------
        DELIdCounter
        """
        if isinstance(other, DELIdCounter):
            new_counter = DELIdCounter()
            new_counter.count = self.count + other.count
            return new_counter
        raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")

    def add_id(self):
        """Increment the ID counter by 1"""
        self.count += 1

    def get_degen_count(self) -> int:
        """
        Return the degen count

        Notes
        -----
        For DELIdCounters, the degen count is the same as the
        raw count.

        Returns
        -------
        int
            the degen count
        """
        return self.count

    def get_raw_count(self) -> int:
        """
        Return the raw count for the degen counter

        Returns
        -------
        int
            the raw count
        """
        return self.count

    def to_dict(self) -> dict:
        """Covert the counter to a JSON compatible dictionary representation"""
        return {"raw_count": self.count}


class UMIDegenerate(abc.ABC):
    """Base class for UMI Degenerators"""

    @abc.abstractmethod
    def add_umi(self, umi: UMI) -> bool:
        """
        Given a UMI, add it to the set and track its count

        Will return a bool indicating if the UMI was added

        Parameters
        ----------
        umi: UMI
            the UMI to add

        Returns
        -------
        bool
            True if UMI is added (was novel)
            False if it is not added (not novel)
        """
        raise NotImplementedError()


class UMICounter(UMIDegenerate, DELCounter):
    """
    Counts the number of unique UMI occurrences

    Attributes
    ----------
    umis: set[str]
        the set of UMI already observed
    """

    def __init__(self):
        """Initialize a UMI Counter"""
        self.umis: set[str] = set()
        self._raw_count = 0

    def __add__(self, other):
        """
        Add two UMICounters together by merging the UMI set and summing raw counts

        Parameters
        ----------
        other: UMICounter
            the other counter to add

        Returns
        -------
        UMICounter
        """
        if isinstance(other, UMICounter):
            new_counter = UMICounter()
            new_counter.umis = self.umis.union(other.umis)
            new_counter._raw_count = self._raw_count + other._raw_count
            return new_counter
        raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")

    def add_umi(self, umi: UMI) -> bool:
        """
        Given a UMI, add it to the set and track its count

        Will return a bool indicating if the UMI was added

        Parameters
        ----------
        umi: UMI
            the UMI to add

        Returns
        -------
        bool
            True if UMI is added (was novel)
            False if it is not added (not novel)
        """
        self._raw_count += 1
        umi_code = umi.to_ascii_code()
        if umi_code in self.umis:
            return False
        self.umis.add(umi_code)
        return True

    def get_degen_count(self) -> int:
        """
        Return the degen count

        Returns
        -------
        int
            the degen count
        """
        return len(self.umis)

    def get_raw_count(self) -> int:
        """
        Return the raw count for the degen counter

        Returns
        -------
        int
            the raw count
        """
        return self._raw_count

    def to_dict(self) -> dict:
        """
        Convert the counter to a dictionary representation

        Returns
        -------
        dict
            the dictionary representation of the counter
        """
        return {"raw_count": self._raw_count, "umis": list(self.umis)}


class UMICluster(UMIDegenerate, DELCounter):
    """
    Count UMI occurrences with clustering

    Essentially, this will treat all UMIs
    within some minimum distance to each other
    as the same (in a greedy fashion).
    Useful if you have a high error rate in sequencing
    """

    def __init__(self, min_dist: int = 2):
        """
        Initialize a UMICluster

        Parameters
        ----------
        min_dist: int, default = 2
            the minimum levenshtein distance between two UMIs
            to be considered unique
        """
        self.min_dist = min_dist
        self.umis: list[UMI] = list()
        self._raw_count = 0

    def __add__(self, other):
        """
        Add two UMIClusters together by merging the UMI clusters and summing raw counts

        Parameters
        ----------
        other: UMICluster
            the other counter to add

        Returns
        -------
        UMICluster
        """
        if isinstance(other, UMICluster):
            if self.min_dist != other.min_dist:
                raise ValueError(
                    f"cannot add UMIClusters with different min_dist values: "
                    f"{self.min_dist} and {other.min_dist}"
                )
            new_counter = UMICluster(min_dist=self.min_dist)
            _all_umis = self.umis + other.umis
            for _umi in _all_umis:
                new_counter.add_umi(_umi)
            new_counter._raw_count = self._raw_count + other._raw_count
            return new_counter
        raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")

    def add_umi(self, umi: UMI) -> bool:
        """
        Given a UMI, add it to the set and track its count

        Will return a bool indicating if the UMI was added

        Parameters
        ----------
        umi: UMI
            the UMI to add

        Returns
        -------
        bool
            True if UMI is added (was novel)
            False if it is not added (not novel)
        """
        self._raw_count += 1
        for existing_umi in self.umis:
            # if too close, fail to add new UMI
            if levenshtein_distance(existing_umi.umi_tag, umi.umi_tag) < self.min_dist:
                return False
        self.umis.append(umi)
        return True

    def get_degen_count(self) -> int:
        """
        Return the degen count

        Returns
        -------
        int
            the degen count
        """
        return len(self.umis)

    def get_raw_count(self) -> int:
        """
        Return the raw count for the degen counter

        Returns
        -------
        int
            the raw count
        """
        return self._raw_count

    def to_dict(self) -> dict:
        """
        Convert the counter to a dictionary representation

        Returns
        -------
        dict
            the dictionary representation of the counter
        """
        return {"raw_count": self._raw_count, "umis": list(self.umis)}


class DELIdUmiCounter(DELCounter):
    """
    Count the number of unique UMI reads for a single DEL ID

    Supports both a UMI clustering mode and basic mode.
    Clustering mode will consider any two UMIs within a given
    minimum distance of each other as the same.
    If you have high error rates in your sequencing
    use the clustering mode. Otherwise, use the basic mode.

    Notes
    -----
    UMI stands for "unique molecular identifier"
    UMI degen is crucial for effective noise reduction
    of the PCR amplification step in DEL as it removes
    uneven amplification of a single read, since now
    we know that each unique read has a unique UMI
    """

    def __init__(self, umi_clustering: bool = False, min_umi_cluster_dist: int = 2):
        """
        Initialize DELIdUmiCounter

        Parameters
        ----------
        umi_clustering: bool, default = False
            turn on UMI clustering
        min_umi_cluster_dist: int, default = 2
            minimum levenshtein distance between two UMIs to be considered unique
            ignored if `umi_clustering` is False
        """
        self.umis: UMICluster | UMICounter
        if umi_clustering:
            self.umis = UMICluster(min_dist=min_umi_cluster_dist)
        else:
            self.umis = UMICounter()

    def __add__(self, other):
        """
        Add two DELIdUmiCounters together by merging the UMI clusters and summing raw counts

        Parameters
        ----------
        other: DELIdUmiCounter
            the other counter to add

        Returns
        -------
        DELIdUmiCounter
        """
        if isinstance(other, DELIdUmiCounter):
            new_counter = DELIdUmiCounter()
            new_counter.umis = self.umis + other.umis
            return new_counter
        raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")

    def add_umi(self, umi: UMI) -> bool:
        """
        Add a UMI to the umi counter

        Will automatically handle the UMI already existing

        Parameters
        ----------
        umi: UMI
            the UMI to add

        Returns
        -------
        bool
            `True` if UMI is added has not been seen yet, else `False`
        """
        return self.umis.add_umi(umi)

    def get_degen_count(self) -> int:
        """
        Return the degen count

        Returns
        -------
        int
            the degen count
        """
        return self.umis.get_degen_count()

    def get_raw_count(self) -> int:
        """
        Return the raw count for the degen counter

        Returns
        -------
        int
            the raw count
        """
        return self.umis.get_raw_count()

    def to_dict(self) -> dict:
        """
        Convert the counter to a dictionary representation

        Returns
        -------
        dict
            the dictionary representation of the counter
        """
        return self.umis.to_dict()


# to help with type hinting
DEL_COUNTER_TYPE = TypeVar("DEL_COUNTER_TYPE", bound=DELCounter)


class DELCollectionCounter(abc.ABC, Generic[DEL_COUNTER_TYPE]):
    """Base class for all DELibraryCollection degen counters"""

    del_counter: dict[str, defaultdict[DecodedDELCompound, DEL_COUNTER_TYPE]]

    def __len__(self):
        """Get total number of unique DELs in the counter"""
        return sum([len(barcodes) for barcodes in self.del_counter.values()])

    @abc.abstractmethod
    def __add__(self, other):
        """
        Add two DELCollectionCounter together by merging their counters

        Parameters
        ----------
        other: DELCollectionCounter
            the other counter to add

        Returns
        -------
        DELCollectionCounter
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def count_barcode(self, barcode: DecodedDELCompound) -> bool:
        """
        Given a barcode, add it to the current degen count

        Will handle separating based on DEL ID and UMI
        correction (if UMI is used)

        Parameters
        ----------
        barcode: DecodedDELCompound
            the barcode to add to the counter

        Returns
        -------
        bool
            `True` is added barcode is new (not degenerate), else `False`
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def to_json(self, path: str | Path, compress: bool = False):
        """
        Convert the counter to a JSON serializable dict

        Returns
        -------
        dict
            the JSON serializable dict representation of the counter
        """
        raise NotImplementedError()


class DELCollectionIdUmiCounter(DELCollectionCounter):
    """
    Degen counter for library collections that have UMI tags

    Handles Degen by DEL ID and UMI.
    Each observed DEL compound will have a count based
    on its degen reads and a count based on how many
    raw reads of the given compound there were
    """

    def __init__(self, umi_clustering: bool = False, min_umi_cluster_dist: int = 2):
        """
        Initialize DELIdUmiCounter

        Parameters
        ----------
        umi_clustering: bool, default = False
            turn on UMI clustering
        min_umi_cluster_dist: int, default = 2
            the minimum levenshtein distance between two UMI
            ignored if `umi_clustering` is False
        """
        # god forgive me for this one
        self.del_counter: dict[str, defaultdict[DecodedDELCompound, DELIdUmiCounter]] = (
            defaultdict(
                lambda: defaultdict(
                    partial(
                        DELIdUmiCounter,
                        umi_clustering=umi_clustering,
                        min_umi_cluster_dist=min_umi_cluster_dist,
                    )
                )
            )
        )
        self.umi_clustering = umi_clustering
        self.min_umi_cluster_dist = min_umi_cluster_dist

    def __getstate__(self):
        """
        Covert state to generic dict with no partial functions when requesting state
        """
        state = self.__dict__.copy()
        state["del_counter"] = {k: dict(v) for k, v in self.del_counter.items()}
        return state

    def __setstate__(self, state):
        """
        Set state for del_counter back to nested defaultdict when reloading state from getstate
        """
        self.__dict__.update(state)
        self.del_counter = defaultdict(
            lambda: defaultdict(
                partial(
                    DELIdUmiCounter,
                    umi_clustering=self.umi_clustering,
                    min_umi_cluster_dist=self.min_umi_cluster_dist,
                )
            ),
            {k: defaultdict(DELIdUmiCounter, v) for k, v in state["del_counter"].items()},
        )

    def __add__(self, other):
        """
        Add two DELCollectionIdUmiCounter objects together

        Notes
        -----
        Adding two counters loaded from a pickle could cause
        memory to grow unexpectedly

        Parameters
        ----------
        other: DELCollectionIdUmiCounter
            the other counter to add

        Returns
        -------
        DELCollectionIdUmiCounter
            the new counter with the sum of the two
        """
        if isinstance(other, DELCollectionIdUmiCounter):
            if self.umi_clustering != other.umi_clustering:
                raise ValueError(
                    f"cannot add DELCollectionIdUmiCounter with "
                    f"different `umi_clustering` values: "
                    f"{self.umi_clustering} and {other.umi_clustering}"
                )
            if self.min_umi_cluster_dist != other.min_umi_cluster_dist:
                raise ValueError(
                    f"cannot add DELCollectionIdUmiCounter with different "
                    f"`min_umi_cluster_dist` values: "
                    f"{self.min_umi_cluster_dist} and {other.min_umi_cluster_dist}"
                )

            new_counter = deepcopy(self)
            for library_id in other.del_counter.keys():
                for compound in other.del_counter[library_id].keys():
                    if compound in new_counter.del_counter[library_id]:
                        new_counter.del_counter[library_id][compound] += other.del_counter[
                            library_id
                        ][compound]
                    else:
                        new_counter.del_counter[library_id][compound] = deepcopy(
                            other.del_counter[library_id][compound]
                        )
            return new_counter
        raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")

    def count_barcode(self, barcode: DecodedDELCompound) -> bool:
        """
        Given a barcode, add it to the current degen count

        Will handle separating based on DEL ID and UMI
        correction

        Parameters
        ----------
        barcode: DecodedDELCompound
            the barcode to add to the counter

        Returns
        -------
        bool
            `True` is added barcode is new (not degenerate), else `False`

        Raises
        ------
        RuntimeError
            if the decoded barcode is missing a UMI
        """
        if barcode.umi is None:
            raise RuntimeError("cannot UMI degen on read missing UMI")
        else:
            return self.del_counter[barcode.library.library_id][barcode].add_umi(barcode.umi)

    def to_json(self, path: str | Path, compress: bool = False):
        """
        Convert the counter to a JSON serializable dict

        Parameters
        ----------
        path: str | Path
            path to save the JSON file to
        compress: bool, default = False
            if True, will compress the JSON file using gzip
            will be encoded with utf-8

        Returns
        -------
        dict
            the JSON serializable dict representation of the counter
        """
        import json

        # collect data
        _data: dict[str, dict[str, object]] = {}
        for lib_ids, dels in self.del_counter.items():
            _data[lib_ids] = {}
            for compound, counter in dels.items():
                _id = compound.compound_id
                _info = {
                    "lib_id": compound.library.library_id,
                    "bb_ids": [bb.bb_id for bb in compound.building_blocks],
                    "raw_count": counter.umis.get_raw_count(),
                    "umis": list(counter.umis.umis),
                }
                _data[lib_ids][_id] = _info

        if compress:
            with gzip.open(path, "wt", encoding="utf-8") as zipfile:
                json.dump(_data, zipfile)
        else:
            json.dump(_data, open(path, "w"))


class DELCollectionIdCounter(DELCollectionCounter):
    """
    Degen counter for library collectionss without UMI tags

    Handles Degen by DEL ID
    Each observed DEL compound will have a count based
    on its degen reads and a count based on how many
    raw reads of the given compound there were
    """

    def __init__(self):
        self.del_counter: defaultdict[str, defaultdict[DecodedDELCompound, DELIdCounter]] = (
            defaultdict(lambda: defaultdict(DELIdCounter))
        )

    def __getstate__(self):
        """
        Covert state to generic dict with no partial functions when requesting state
        """
        state = self.__dict__.copy()
        state["del_counter"] = {k: dict(v) for k, v in self.del_counter.items()}
        return state

    def __setstate__(self, state):
        """
        Set state for del_counter back to nested defaultdict when reloading state from getstate
        """
        self.__dict__.update(state)
        self.del_counter = defaultdict(
            lambda: defaultdict(DELIdCounter),
            {k: defaultdict(DELIdUmiCounter, v) for k, v in state["del_counter"].items()},
        )

    def __add__(self, other):
        """
        Add two DELCollectionIdCounter together by merging the counters

        Parameters
        ----------
        other: DELCollectionIdCounter
            the other counter to add

        Returns
        -------
        DELCollectionIdCounter
        """
        if isinstance(other, DELCollectionIdCounter):
            new_counter = DELCollectionIdCounter()
            for lib_id in self.del_counter.keys() | other.del_counter.keys():
                for barcode in self.del_counter[lib_id].keys() | other.del_counter[lib_id].keys():
                    new_counter.del_counter[lib_id][barcode] = (
                        self.del_counter[lib_id][barcode] + other.del_counter[lib_id][barcode]
                    )
            return new_counter
        raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")

    def count_barcode(self, barcode: DecodedDELCompound) -> bool:
        """
        Given a barcode, add it to the current degen count

        Will handle separating based on DEL ID

        Parameters
        ----------
        barcode: DecodedDELCompound
            the barcode to add to the counter

        Returns
        -------
        bool
            `True` is added barcode is new (not degenerate), else `False`
        """
        self.del_counter[barcode.library.library_id][barcode].add_id()
        return True  # always not degenerate with this counter

    def to_json(self, path: str | Path, compress: bool = False):
        """
        Convert the counter to a JSON serializable dict

        Parameters
        ----------
        path: str | Path
            path to save the JSON file to
        compress: bool, default = False
            if True, will compress the JSON file using gzip

        Returns
        -------
        dict
            the JSON serializable dict representation of the counter
        """
        import json

        # collect data
        _data: dict[str, dict[str, object]] = {}
        for lib_ids, dels in self.del_counter.items():
            _data[lib_ids] = {}
            for compound, counter in dels.items():
                _id = compound.compound_id
                _info = {
                    "lib_id": compound.library.library_id,
                    "bb_ids": [bb.bb_id for bb in compound.building_blocks],
                    "raw_count": counter.get_raw_count(),
                }
                _data[lib_ids][_id] = _info

        if compress:
            with gzip.open(path, "wt", encoding="utf-8") as zipfile:
                json.dump(_data, zipfile)
        else:
            json.dump(_data, open(path, "w"))
