"""generate hamming encoded tag sets for DEL barcode design"""

from typing import Generator, Optional

from deli.processing.hamming import QuaternaryHammingEncoder

from .tag import Tag


class TagFactory:
    """Class used to generate sets of hamming encoded tags"""

    def __init__(self, encoder: QuaternaryHammingEncoder):
        """
        Initialize a TagFactory

        Parameters
        ----------
        encoder: QuaternaryHammingEncoder
            the encoder used to generate the hamming encoded tag set
        """
        self.next_tag_number = 0
        self.encoder = encoder
        self.sequence_length = self.encoder.total_length

    @property
    def next_tag_number(self):
        """Get the next tag number from the set"""
        ntn = self._next_tag_number
        self._next_tag_number += 1
        return ntn

    @next_tag_number.setter
    def next_tag_number(self, value):
        """Set the tag number"""
        self._next_tag_number = value

    def reset_tag_number(self):
        """Reset tag number to zero"""
        self.next_tag_number = 0

    def create_tag(self, tag_number: Optional[int] = None) -> Tag:
        """Create an encoded tag from the tag number or the next tag number"""
        if tag_number is None:
            tag_number = self.next_tag_number
        return self.encoder.encode(tag_number)

    def create_tags(self, num: Optional[int] = None) -> Generator[Tag, None, None]:
        """Create N tags"""
        if num is None:
            num = self.encoder.max_tags
        for _ in range(num):
            yield self.create_tag()
