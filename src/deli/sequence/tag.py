"""handles hamming encoded tag generation"""

from enum import Enum
from typing import Any, List, Self, Union


class Base(Enum):
    """Enum class for single nucleotides"""

    A = 0
    T = 1
    C = 2
    G = 3

    @classmethod
    def _missing_(cls, value) -> Self:
        """Create a enum from a int or str"""
        if isinstance(value, int):
            return cls(value % 4)
        if value == "A":
            return cls(0)
        elif value == "T":
            return cls(1)
        elif value == "C":
            return cls(2)
        elif value == "G":
            return cls(3)
        else:
            raise ValueError(f"unrecognized enum value {value}")

    @property
    def bin(self) -> str:
        """Convert to the two digit binary rep (00, 01, 10, 11)"""
        bin_repr = bin(self.value)[2:]
        return "0" * (2 - len(bin_repr)) + bin_repr

    def __str__(self) -> str:
        """Return the single letter base as a str"""
        return self.name

    def __add__(self, other: Any) -> int:
        """Add two bases or a base and an int"""
        if isinstance(other, int):
            return self.value + other
        elif isinstance(other, self.__class__):
            return self.value + other.value
        raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")

    def __radd__(self, other: Any) -> int:
        """Right add"""
        return self.__add__(other)

    def __mod__(self, other: int) -> int:
        """Return the module of the base by some other number"""
        return self.value % other

    def __sub__(self, other):
        """Subtract the other base/int and return the base that results"""
        if isinstance(other, int):
            other_num = other
        elif isinstance(other, self.__class__):
            other_num = other.value
        else:
            raise TypeError(
                f"unsupported operand type(s) for -: '{type(self)}' and '{type(other)}'"
            )
        return self.__class__((self.value - other_num) % 4)

    def __rsub__(self, other):
        """Subtract the other base/int and return the base that results"""
        self.__sub__(other)

    def __eq__(self, other) -> bool:
        """Return true if both base are the same else false"""
        return self.value == other.value if isinstance(other, self.__class__) else False


class Tag:
    """Hold a list of bases that make a tag"""

    def __init__(self, bases: Union[List[Base], List[str]]):
        """
        Initialize a Tag object

        Parameters
        ----------
        bases: List[Union[Base, str]]
            the bases in an ordered list
            if a str, will be converted to a Base
        """
        self.bases: List[Base] = [Base(_) if not isinstance(_, Base) else _ for _ in bases]

    def __repr__(self) -> str:
        """Represent the tag as a string of the single letter bases"""
        return self.__str__()

    def __str__(self) -> str:
        """Return the tag as a string of the single letter bases"""
        return "".join(base.name for base in self.bases)

    def __len__(self) -> int:
        """Return length of the tag (number of bases)"""
        return len(self.bases)

    @property
    def bin(self) -> str:
        """Binary representation of a tag."""
        return "".join(base.bin for base in self.bases)

    def __eq__(self, other) -> bool:
        """Return true if the tags are the same"""
        return str(self) == str(other) if isinstance(other, Tag) else False
