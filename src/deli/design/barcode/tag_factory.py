"""Code for building dna tag sets"""

import math
import random
from itertools import product
from typing import Optional

from deli._hamming import BaseQuaternaryHamming
from deli.configure import get_deli_config


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

    def gc_content(self) -> float:
        """
        Calculate the GC content of the tag as a fraction

        Returns
        -------
        float
            the GC content of the tag as a fraction between 0 and 1
        """
        gc_count = sum(1 for base in self.bases if base in "GC")
        return gc_count / len(self) if len(self) > 0 else 0.0

    def num_identical_neighbors(self) -> int:
        """
        Get the number of identical neighboring bases in the tag

        Returns
        -------
        int
            the number of identical neighboring bases in the tag
        """
        count = 0
        for i in range(1, len(self.bases)):
            if self.bases[i] == self.bases[i - 1]:
                count += 1
        return count

    def shannon_entropy(self) -> float:
        """
        Calculate the Shannon entropy of the tag

        Returns
        -------
        float
            the Shannon entropy of the tag, calculated as -sum(p_i * log2(p_i)) for each base
        """
        base_counts = {base: self.bases.count(base) for base in set(self.bases)}
        total_bases = len(self.bases)
        entropy = -sum((count / total_bases) * math.log2(count / total_bases) for count in base_counts.values())
        return entropy

    def basepair_counts(self) -> dict[str, int]:
        """
        Count how many times each base appears in the tag

        Returns
        -------
        dict[str, int]
            a dictionary mapping each base to its count in the tag
        """
        return {base: self.bases.count(base) for base in set(self.bases)}


class QuaternaryHammingTagFactory(BaseQuaternaryHamming):
    """Used to generate a set of DNA tags that is Hamming encoded"""

    def __init__(
        self,
        desired_tag_size: Optional[int] = None,
        desired_number_of_tags: Optional[int] = None,
        use_extra_parity: bool = False,
    ):
        """
        Initialize a QuaternaryHammingTagFactory for DNA tag set generation

        Parameters
        ----------
        desired_tag_size: Optional[int], default=None
            the size you would like the tag (for example 8 for an 8 base tag)
        desired_number_of_tags: Optional[int], default=None
            the number of tags you will need
        use_extra_parity: bool, default=False
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

        if desired_number_of_tags is not None:
            if desired_tag_size is None:
                _desired_tag_size = float("inf")
            else:
                _desired_tag_size = float(desired_tag_size)

            # determine the message size and parity size needed to support the desired number of tags
            self.message_size = math.ceil(math.log(desired_number_of_tags, 4))
            self.parity = 0
            while (2**self.parity) < (self.message_size + self.parity + 1):
                self.parity += 1
            self.total_size = self.message_size + self.parity + self.use_extra_parity

            if self.total_size > _desired_tag_size:
                raise ValueError(
                    f"cannot build hamming tag set with {desired_number_of_tags} "
                    f"tags with desired length of {desired_tag_size}" + "and extra parity bit" * use_extra_parity
                )
        elif desired_tag_size is not None:
            self.parity = math.ceil(math.log(desired_tag_size, 2))
            self.message_size = desired_tag_size - self.parity - use_extra_parity
            self.total_size = desired_tag_size
        else:
            raise ValueError("must specify at least one of desired_tag_size or desired_number_of_tags")

        self.max_tags = 4**self.message_size

    def _encode(self, bases: list[int]) -> list[int]:
        """Given a list of ints, encode it into a hamming vector"""
        _bases = bases.copy()
        full_message = [
            _bases.pop(0) if not (i & (i - 1) == 0) else 0 for i in range(self.total_size - self.use_extra_parity + 1)
        ]

        parity = [self._get_sub_parity(full_message, _) for _ in range(self.parity, 0, -1)]

        for i in range(self.parity):
            full_message[2**i] = parity[i]

        if self.use_extra_parity:
            full_message[0] = self._get_global_parity(full_message)
        else:
            full_message = full_message[1:]

        return full_message

    def _generate_tag(self, combo: list[int]) -> Tag:
        """Given a list of int encoded bases, generate the hamming encoded tag"""
        bases = [get_deli_config().int_2_nuc[x] for x in self._encode(list(combo))]
        return Tag("".join(bases))

    def _get_tag_combinations(self) -> list[list[int]]:
        """Return all possible combinations of message bits"""
        return [
            list(combo) for combo in product([0, 1, 2, 3], repeat=self.message_size)
        ]  # not very memory efficient...

    def generate_all_tags(self) -> list[Tag]:
        """
        Generate the full set of tags in the factory

        Returns
        -------
        List[Tag]
            the full set of tags in the factory
        """
        tag_combinations = self._get_tag_combinations()
        return [self._generate_tag(combo) for combo in tag_combinations]

    def generate_tags(self, num_tags: int, shuffle: bool = True) -> list[Tag]:
        """
        Generate a random set of tags from the factory

        Parameters
        ----------
        num_tags: int
            the number of tags to generate
        shuffle: bool, default=True
            whether to shuffle the generated tags

        Returns
        -------
        List[Tag]
            a random set of tags from the factory
        """
        if num_tags > self.max_tags:
            raise ValueError(f"cannot generate {num_tags} unique tags with this factory; max tags is {self.max_tags}")

        tag_combos = self._get_tag_combinations()

        if shuffle:
            random.shuffle(tag_combos)

        return [self._generate_tag(combo) for combo in tag_combos[:num_tags]]
