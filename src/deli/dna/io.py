"""for reading in sequence files"""

import abc
import os
from collections.abc import Iterator
from pathlib import Path
from typing import Union

import dnaio
from dnaio import SequenceRecord


class Reader(abc.ABC):
    """base class for all seq reader objects"""

    @abc.abstractmethod
    def __iter__(self) -> Iterator[SequenceRecord]:
        """Yield all sequences in the file"""
        pass

    def read_sequences(self) -> Iterator[SequenceRecord]:
        """Yield all sequences in the fill, wrapper to __iter__"""
        return self.__iter__()


class SequenceReader(Reader):
    """Read sequence data from a file on the fly"""

    def __init__(self, sequence_file: Union[str, Path]):
        """
        Initialize a SequenceReader object.

        Parameters
        ----------
        sequence_file: Union[str, Path]
            path to the sequence file to read
        """
        self.sequence_file = Path(sequence_file)

    def __iter__(self) -> Iterator[SequenceRecord]:
        """Yield all sequences in the file"""
        with dnaio.open(self.sequence_file) as f:
            for record in f:
                yield record


class SequenceDirectoryReader:
    """Read all sequences from all files in a directory"""

    def __init__(self, sequence_dir: Union[str, Path]):
        """
        Initialize a SequenceDirectoryReader object.

        Parameters
        ----------
        sequence_dir: Union[str, Path]
            path to the directory with sequence files to read
        """
        self.sequence_dir = Path(sequence_dir)

    def get_sequence_files(self) -> list[Path]:
        """
        Get all the fastq files from the directory.

        Returns
        -------
        list[Path]
        """
        _files: list[Path] = list()
        for file in os.listdir(self.sequence_dir):
            if file.endswith(".fastq"):
                _files.append(Path(os.path.join(self.sequence_dir, file)))
        return _files

    def __iter__(self) -> Iterator[SequenceRecord]:
        """Yield all sequences from all files in the directory"""
        for sequence_file in self.get_sequence_files():
            with dnaio.open(sequence_file) as f:
                for record in f:
                    yield record
