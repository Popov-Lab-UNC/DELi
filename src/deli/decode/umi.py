"""umi functions and classes"""


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
