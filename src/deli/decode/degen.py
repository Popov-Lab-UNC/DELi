"""code for degenerating barcode reads"""

import abc
from collections import defaultdict
from functools import partial

from Levenshtein import distance as levenshtein_distance

from .decoder import DecodedBarcode
from .umi import UMI


class UMIDegenerate(abc.ABC):
    """Base class for UMI Degenerators"""

    raw_count: int = 0

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

    @abc.abstractmethod
    def get_count(self) -> int:
        """Return the degen count"""
        raise NotImplementedError()

    def get_raw_count(self) -> int:
        """Return the raw count for the degenerator"""
        return self.raw_count


class UMICounter(UMIDegenerate):
    """
    Counts the number of unique UMI occurrences

    Attributes
    ----------
    raw_count: int
        the number of times a UMI add was attempted
    """

    def __init__(self):
        """Initialize a UMI Counter"""
        self.umis: set[str] = set()

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
        self.raw_count += 1
        if umi.umi_tag in self.umis:
            return False
        self.umis.add(umi.umi_tag)
        return True

    def get_count(self) -> int:
        """
        Return the degen count

        Returns
        -------
        int
            the degen count
        """
        return len(self.umis)


class UMICluster(UMIDegenerate):
    """
    Count UMI occurrences with clustering

    Essentially this will treat all UMIs
    within some minimum distance to each other
    as the same (in a greedy fashion).
    Useful if you have a high error rate in sequencing

    Attributes
    ----------
    raw_count: int
        the number of times a UMI add was attempted
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
        self.raw_count += 1
        for existing_umi in self.umis:
            # if too close, fail to add new UMI
            if levenshtein_distance(existing_umi.umi_tag, umi.umi_tag) < self.min_dist:
                return False
        self.umis.append(umi)
        return True

    def get_count(self) -> int:
        """
        Return the degen count

        Returns
        -------
        int
            the degen count
        """
        return len(self.umis)


class DELIdUmiCounter:
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
        self.umis: UMIDegenerate
        if umi_clustering:
            self.umis = UMICluster(min_dist=min_umi_cluster_dist)
        else:
            self.umis = UMICounter()

    def add_umi(self, umi: UMI):
        """
        Add a UMI to the umi counter

        Will automatically handle the UMI already existing

        Parameters
        ----------
        umi: UMI
            the UMI to add
        """
        self.umis.add_umi(umi)


class DELibraryPoolCounter(abc.ABC):
    """Base class for all DELibraryPool degen counters"""

    @abc.abstractmethod
    def count_barcode(self, barcode: DecodedBarcode):
        """
        Given a barcode, add it to the current degen count

        Will handle separating based on DEL ID and UMI
        correction (if UMI is used)

        Parameters
        ----------
        barcode: DecodedBarcode
            the barcode to add to the counter
        """
        raise NotImplementedError()


class DELibraryPoolIdUmiCounter:
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
        self.del_counter: defaultdict[str, DELIdUmiCounter] = defaultdict(
            partial(
                DELIdUmiCounter,
                umi_clustering=umi_clustering,
                min_umi_cluster_dist=min_umi_cluster_dist,
            )
        )

    def count_barcode(self, barcode: DecodedBarcode):
        """
        Given a barcode, add it to the current degen count

        Will handle separating based on DEL ID and UMI
        correction

        Parameters
        ----------
        barcode: DecodedBarcode
            the barcode to add to the counter
        """
        if barcode.umi is None:
            raise RuntimeError("cannot UMI degen on read missing UMI")
        else:
            self.del_counter[barcode.get_id()].add_umi(barcode.umi)


class DELibraryPoolIdCounter:
    """
    Degen counter for library pools without UMI tags

    Handles Degen by DEL ID
    Each observed DEL compound will have a count based
    on its degen reads and a count based on how many
    raw reads of the given compound there were
    """

    def __init__(self):
        self.del_counter = defaultdict(int)

    def count_barcode(self, barcode: DecodedBarcode):
        """
        Given a barcode, add it to the current degen count

        Will handle separating based on DEL ID

        Parameters
        ----------
        barcode: DecodedBarcode
            the barcode to add to the counter
        """
        self.del_counter[barcode.get_id()] += 1
