"""Code for building dna tag sets"""

import math
from typing import Optional

from deli._hamming import BaseQuaternaryHamming


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


class QuaternaryHammingTagFactory(BaseQuaternaryHamming):
    """Used to generate a set of DNA tags that is Hamming encoded"""

    def __init__(
        self,
        desired_number_of_tags: int,
        desired_tag_size: Optional[int] = None,
        use_extra_parity: bool = True,
    ):
        """
        Initialize a QuaternaryHammingTagFactory for DNA tag set generation

        Parameters
        ----------
        desired_number_of_tags: int
            the number of tags you will need
        desired_tag_size: Optional[int], default=None
            the size you would like the tag
            must be at least big enough to cover desired number of tags
        use_extra_parity: bool, default=True
            use the extra parity bit encoding to allow for detection of
            double errors

        Raises
        ------
        ValueError
            if the desired number of tags requires a tag size greater than
            the desired tag size

        Returns
        -------
        Self
        """
        super().__init__()
        self.use_extra_parity = use_extra_parity
        self.message_size = math.ceil(math.log(desired_number_of_tags, 4))

        self.parity = 0
        while (2**self.parity) < (self.message_size + self.parity + 1):
            self.parity += 1

        self.total_size = self.message_size + self.parity + self.use_extra_parity

        if desired_tag_size is None:
            _desired_tag_size = self.total_size
        else:
            _desired_tag_size = desired_tag_size

        if _desired_tag_size < self.total_size:
            raise ValueError(
                f"cannot build hamming tag set with {desired_number_of_tags} "
                f"tags with desired length of {_desired_tag_size}"
            )

        self.max_tags = 4**self.message_size

    def _encode(self, bases: list[int]) -> list[int]:
        """Given a list of ints, encode it into a hamming vector"""
        _bases = bases.copy()
        full_message = [
            _bases.pop(0) if not (i & (i - 1) == 0) else 0
            for i in range(self.total_size - self.use_extra_parity + 1)
        ]

        parity = [self._get_sub_parity(full_message, _) for _ in range(self.parity, 0, -1)]

        for i in range(self.parity):
            full_message[2**i] = parity[i]

        if self.use_extra_parity:
            full_message[0] = self._get_global_parity(full_message)
        else:
            full_message = full_message[1:]

        return full_message
