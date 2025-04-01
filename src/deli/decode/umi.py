"""umi functions and classes"""

import abc

from Levenshtein import distance as levenshtein_distance


class UMI:
    """Call class for UMIs"""

    def __init__(self, umi_tag: str):
        """
        Initialize a UMICall

        Parameters
        ----------
        umi_tag: str
            the DNA sequence of the umi tag
        """
        self.umi_tag = umi_tag

    def __str__(self):
        """Return the DNA sequence of the UMI as a string"""
        return self.umi_tag


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


class UMICounter:
    """Counts the number of unique UMI occurrences"""

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
        if umi.umi_tag in self.umis:
            return False
        self.umis.add(umi.umi_tag)
        return True


class UMICluster:
    """
    Count UMI occurrences with clustering

    Essentially this will treat all UMIs
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
        for existing_umi in self.umis:
            # if too close, fail to add new UMI
            if levenshtein_distance(existing_umi.umi_tag, umi.umi_tag) < self.min_dist:
                return False
        self.umis.append(umi)
        return True
