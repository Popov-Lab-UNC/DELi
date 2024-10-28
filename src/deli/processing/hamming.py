"""Handles hamming correction of DEL barcode tags"""

from math import ceil
from typing import Iterator, List, Union

from deli.sequence.tag import Base, Tag


class DecodeError(Exception):
    """raised when the decoder fails to decode a barcode"""

    pass


class QuaternaryHammingEncoder:
    """Generates and decodes quaternary hamming codes"""

    def __init__(self, total_length: int):
        """
        Initialize the Quaternary Hamming Encoder.

        Parameters
        ----------
        total_length: int
            the length (num nucleotides) of the encoded tags
        """
        self.total_length = total_length
        for i in range(2, 10):
            if total_length <= (2**i):
                self.parity = i
                self.message_length = self.total_length - self.parity - 1
                break
        self.max_tags = 4**self.message_length

    def _get_sub_parity(self, tag_list: list[Base], order: int = 0) -> int:
        """Get the sub parity of the tag"""
        if len(tag_list) != self.total_length:
            raise ValueError("tag not correct length for encoder")
        if order == 0:
            return (4 - sum(tag_list)) % 4
        new_list = []
        s = ceil(self.total_length / (2**order))
        for i in range((2**order)):
            if i % 2 == 1:
                new_list.extend(tag_list[s * i : s * (i + 1)])
        return (4 - sum(new_list)) % 4

    @staticmethod
    def _compute_parity(base_list: List[Base]) -> Base:
        """Get the full parity of the tag"""
        return Base(4 - sum(base_list) % 4)

    def encode(self, tag_number: int) -> Tag:
        """
        Generate a given hamming code

        Parameters
        ----------
        tag_number: int
            the tag number in the code series to generate

        Returns
        -------
        Tag
        """
        if tag_number >= self.max_tags:
            raise ValueError("tag number exceeds range of encoder")
        message: Iterator[Base] = iter(
            [Base((tag_number % (4 ** (i + 1))) // (4**i)) for i in range(0, self.message_length)][
                ::-1
            ]
        )
        full_message: List[Base] = [
            message.__next__() if not (i & (i - 1) == 0) else Base(0)
            for i in range(self.total_length)
        ]
        parity: List[Base] = [
            Base(self._get_sub_parity(list(full_message), _)) for _ in range(self.parity, 0, -1)
        ]
        for i in range(self.parity):
            full_message[2**i] = parity[i]
        full_message[0] = Base(self._get_sub_parity(list(full_message), 0))
        return Tag(full_message)

    def decode(self, tag: Union[List[Base], Tag]) -> Tag:
        """
        Decode a tag and correct any correctable errors

        Parameters
        ----------
        tag: Union[List[Base], Tag]
            the tag to decode/correct

        Raises
        ------
        DecodeError
            when more than 2 errors are detected in the barcode
            and it is not correctable

        Returns
        -------
        Tag
        """
        if isinstance(tag, Tag):
            _tag = tag.bases
        else:
            _tag = tag

        sub_parity = [self._get_sub_parity(_tag, order=p) for p in range(self.parity, 0, -1)][::-1]
        overall_parity = self._get_sub_parity(_tag, order=0)

        if len(set(sub_parity)) > 2:
            raise DecodeError("two or more errors detected")
        if (len(set(sub_parity)) == 2) and (0 not in sub_parity):
            raise DecodeError("two or more errors detected")
        if (overall_parity == 0) and (len(set(sub_parity)) == 2):
            raise DecodeError("two or more errors detected")
        corrected_tag = _tag.copy()
        if max(sub_parity) != 0:
            pos_corr = int("".join([str(1 if c else 0) for c in sub_parity]), 2)
            corrected_tag[pos_corr] -= 4 - max(sub_parity)
            overall_parity = self._get_sub_parity(corrected_tag, order=0)
        elif (
            overall_parity != 0
        ):  # need to correct the parity bit in the case that this is the one that is misread
            corrected_tag[0] = Base((4 - sum(corrected_tag[1:])) % 4)
            overall_parity = self._get_sub_parity(corrected_tag, order=0)
        if overall_parity != 0:
            raise DecodeError("two or more errors detected")
        return Tag(corrected_tag)
