"""code for degenerating barcode reads"""

import abc
import os
from collections import defaultdict
from collections.abc import Iterator
from functools import partial
from typing import Literal

from Levenshtein import distance as levenshtein_distance

from .decoder import DecodedBarcode
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
        if umi.umi_tag in self.umis:
            return False
        self.umis.add(umi.umi_tag)
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


class DELibraryPoolCounter(abc.ABC):
    """Base class for all DELibraryPool degen counters"""

    del_counter: dict

    @abc.abstractmethod
    def count_barcode(self, barcode: DecodedBarcode) -> bool:
        """
        Given a barcode, add it to the current degen count

        Will handle separating based on DEL ID and UMI
        correction (if UMI is used)

        Parameters
        ----------
        barcode: DecodedBarcode
            the barcode to add to the counter

        Returns
        -------
        bool
            `True` is added barcode is new (not degenerate), else `False`
        """
        raise NotImplementedError()

    def _get_barcode_rows(self) -> Iterator[tuple[str, str, str, str]]:
        """
        Covert the barcode-count key value pair to a file row

        Should return a tuple with the following order:
        - DEL ID
        - SMILES (empty string if not available)
        - Raw count
        - UMI corrected count (same as raw count if no UMI)

        Yields
        ------
        tuple[str, str, str, str]
            the row tuple for a given DEL ID
        """
        for decoded_barcode, counter in self.del_counter.items():
            yield (
                decoded_barcode.id,
                decoded_barcode.get_smiles(""),
                str(counter.get_raw_count()),
                str(counter.get_degen_count()),
            )

    def to_file(self, path: str | os.PathLike, file_format: Literal["csv", "tsv"] = "tsv") -> None:
        """
        Write the Pool Counter to a human-readable file

        Each row will be a unique DEL ID with 4 columns:
        - DEL ID
        - SMILES (empty string if not available)
        - Raw count
        - UMI corrected count (same as raw count if no UMI)

        Parameters
        ----------
        path: str | os.PathLike
            path to save file to
        file_format: Literal["csv", "tsv"] = "tsv"
            which file format to write to
        """
        delimiter: str
        if file_format == "tsv":
            delimiter = "\t"
        elif file_format == "csv":
            delimiter = ","
        else:
            raise ValueError(f"file format '{file_format}' not recognized")

        with open(path, "w") as f:
            # write header
            f.write(delimiter.join(["DEL_ID", "SMILES", "RAW_COUNT", "UMI_CORRECTED_COUNT"]))
            for _row in self._get_barcode_rows():
                f.write(delimiter.join(_row) + "\n")


class DELibraryPoolIdUmiCounter(DELibraryPoolCounter):
    """
    Degen counter for library pools that have UMI tags

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
        self.del_counter: defaultdict[DecodedBarcode, DELIdUmiCounter] = defaultdict(
            partial(
                DELIdUmiCounter,
                umi_clustering=umi_clustering,
                min_umi_cluster_dist=min_umi_cluster_dist,
            )
        )

    def count_barcode(self, barcode: DecodedBarcode) -> bool:
        """
        Given a barcode, add it to the current degen count

        Will handle separating based on DEL ID and UMI
        correction

        Parameters
        ----------
        barcode: DecodedBarcode
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
            return self.del_counter[barcode].add_umi(barcode.umi)


class DELibraryPoolIdCounter(DELibraryPoolCounter):
    """
    Degen counter for library pools without UMI tags

    Handles Degen by DEL ID
    Each observed DEL compound will have a count based
    on its degen reads and a count based on how many
    raw reads of the given compound there were
    """

    def __init__(self):
        self.del_counter: defaultdict[DecodedBarcode, DELIdCounter] = defaultdict(DELIdCounter)

    def count_barcode(self, barcode: DecodedBarcode) -> bool:
        """
        Given a barcode, add it to the current degen count

        Will handle separating based on DEL ID

        Parameters
        ----------
        barcode: DecodedBarcode
            the barcode to add to the counter

        Returns
        -------
        bool
            `True` is added barcode is new (not degenerate), else `False`
        """
        self.del_counter[barcode].add_id()
        return True  # always not degenerate with this counter
