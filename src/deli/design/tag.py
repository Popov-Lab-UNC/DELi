"""handles hamming encoded bases generation"""


class Tag:
    """Hold a list of bases that make a bases"""

    def __init__(self, bases: str):
        """
        Initialize a Tag object

        Parameters
        ----------
        bases: List[Union[Base, str]]
            the bases in an ordered list
            if a str, will be converted to a Base
        """
        self.bases = bases

    def __repr__(self) -> str:
        """Represent the bases as a string of the single letter bases"""
        return self.__str__()

    def __str__(self) -> str:
        """Return the bases as a string of the single letter bases"""
        return self.bases

    def __len__(self) -> int:
        """Return length of the bases (number of bases)"""
        return len(self.bases)

    def __eq__(self, other) -> bool:
        """Return true if the tags are the same"""
        return str(self) == str(other) if isinstance(other, Tag) else False
