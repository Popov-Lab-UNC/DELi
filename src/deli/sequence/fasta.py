"""reading and parsing fasta/q files"""

import os
from typing import Iterator, List, Self, Union

from .utils import reverse_compliment


class FastaSequence(object):
    """The sequence and metadata from a FASTA sequence"""

    def __init__(self, read_tag: str, read_id: int, sequence: str):
        """
        Initialize a FastaSequence object

        Parameters
        ----------
        read_tag: str
            the read bases given in the FASTA file
        read_id: int
            a numerical idx assigned to the read based on location in file
            starts from 0
        sequence: str
            the sequence
        """
        self.read_tag = read_tag
        self.read_id = read_id
        self.sequence = sequence

    def subset(self, start: int, end: int) -> Self:
        """
        get copy of the object sequence subset between start and end

        Parameters
        ----------
        start: int
            the start of the subset (inclusive)
        end: int
            the end of the subset (exclusive)

        Returns
        -------
        FastqSequence
        """
        return self.__class__(self.read_tag, self.read_id, self.sequence[start:end])

    def revcomp(self) -> Self:
        """Get the reverse complement of the sequence"""
        return self.__class__(self.read_tag, self.read_id, reverse_compliment(self.sequence))

    def __str__(self):
        """Return attribute as strings seperated by whitespace"""
        return self.sequence

    def __len__(self):
        """Returns length of the sequence"""
        return len(self.sequence)


class FastqSequence(FastaSequence):
    """
    The sequence and metadata from a FASTQ sequence

    Attributes
    ----------
    quality: list[int]
        the quality of the sequence in a numerical list format
        0 is worst
    """

    def __init__(self, read_tag: str, read_id: int, sequence: str, quality: str):
        """
        Initialize a FastqSequence object

        Parameters
        ----------
        read_tag: str
            the read bases given in the FASTA file
        read_id: int
            a numerical idx assigned to the read based on location in file
            starts from 0
        sequence: str
            the sequence
        quality: str
            the quality of the sequence in FASTQ format
        """
        super().__init__(read_tag, read_id, sequence)
        self._quality = quality

    @property
    def quality(self) -> List[int]:
        """Numerical quality of fastq sequence"""
        return [ord(q) - 33 for q in self._quality]

    def get_read_quality(self, idx: int) -> int:
        """
        Given a idx in the sequence get the quality of that read

        Parameters
        ----------
        idx: int
            the index of the sequence

        Returns
        -------
        quality: int
            the quality of the read at that index
        """
        return ord(self._quality[idx]) - 33

    def get_average_quality(self) -> float:
        """
        Get the average quality of the whole sequence

        Returns
        -------
        quality: float
            the average quality of the whole sequence
        """
        return sum(self.quality) / len(self)

    def subset(self, start: int, end: int) -> Self:
        """
        get copy of the object sequence subset between start and end

        Parameters
        ----------
        start: int
            the start of the subset (inclusive)
        end: int
            the end of the subset (exclusive)

        Returns
        -------
        FastqSequence
        """
        return self.__class__(
            self.read_tag, self.read_id, self.sequence[start:end], self._quality[start:end]
        )

    def revcomp(self) -> Self:
        """Get the reverse complement of the sequence"""
        return self.__class__(
            self.read_tag,
            self.read_id,
            reverse_compliment(self.sequence),
            quality=self._quality[::-1],
        )


def yield_fastq(fastq_file: Union[str, os.PathLike]) -> Iterator[FastqSequence]:
    """
    Read a FASTQ file and yield the sequences

    Parameters
    ----------
    fastq_file: str
        path to the FASTQ file

    Returns
    -------
    Iterator[FastqSequence]
    """
    _id_counter = 0
    with open(fastq_file) as file:
        for line in file:
            if line.startswith("@"):
                seq = file.readline().strip()
                _ = file.readline()  # junk line
                quality = file.readline().strip()
                yield FastqSequence(
                    read_tag=line.strip(), read_id=_id_counter, sequence=seq, quality=quality
                )
                _id_counter += 1


def read_fastq(fastq_file: Union[str, os.PathLike]) -> List[FastqSequence]:
    """
    Read a FASTQ file and return a list of the sequences

    Parameters
    ----------
    fastq_file: str
        path to the FASTQ file

    Returns
    -------
    List[FastqSequence]
    """
    return list(yield_fastq(fastq_file))
