"""class/functions for UMI barcode sections"""

from typing_extensions import Self

from Levenshtein import distance


class Umi:
    """Class for holding info on the UMI barcode section"""

    def __init__(self, umi: str):
        """
        Initialize the UMI barcode section

        Parameters
        ----------
        umi: str
            the DNA bases of the umi
        """
        self.umi = umi

    def __eq__(self, other) -> bool:
        """Check if two UMIs have the same barcode"""
        if isinstance(other, Umi):
            return self.umi == other.umi
        return False

    def __len__(self) -> int:
        """Get the length of the UMI barcode"""
        return len(self.umi)

    def is_nearby(self, other: Self, cutoff: int = 1) -> bool:
        """
        If the other UMI barcode is "nearby" based on a distance cutoff

        Notes
        -----
        Nearby means the levenshtein distance between the two UMI barcodes
        is equal to or below the cutoff

        This is used to cluster UMIs to help minimize possible noise
        cause by misread in the UMI region of the barcode

        You can calculate the rough odds that two UMIs will have a
        levenshtein distance of X between them with:

        ```
        (len(my_umi) * 3) / (len(my_umi) ** 4)
        ```

        This assumes only a single SNP occurs and not a indel or more
        than 1 mutation. It is likely lower than the true
        chance, but not by much as double misreads are rare

        Most Umi's are around a length of 8 to 10. The odds are roughly:
        - 8:  24/4096  = 0.6%
        - 9:  36/6561  = 0.5%
        - 10: 40/10000 = 0.4%

        Parameters
        ----------
        other: Umi
            the Umi you want to check is nearby
        cutoff: int, default = 1
            the levenshtein distance cutoff

        Returns
        -------
        bool
        """
        if distance(self.umi, other.umi) <= cutoff:
            return True
        else:
            return False
