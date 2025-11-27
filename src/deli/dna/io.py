"""for reading in sequence files"""

import abc
import glob
import os
from collections.abc import Iterator
from pathlib import Path
from typing import Sequence

import dnaio
from dnaio import SequenceRecord
from dnaio.exceptions import UnknownFileFormat


def _check_if_seq_file(file_name: str) -> bool:
    """
    Check if a file is a sequence file based on its extension.

    Covers common extensions for fasta and fastq files.

    This is lifted from dnaio since we are just wrapping that library.

    Parameters
    ----------
    file_name: str
        the name of the file to check

    Returns
    -------
    bool
        True if the file is a sequence file, False otherwise
    """
    name = file_name.lower()
    for ext in (".gz", ".xz", ".bz2", ".zst"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    name, ext = os.path.splitext(name)
    if ext in [".fasta", ".fa", ".fna", ".csfasta", ".csfa", ".fastq", ".fq"]:
        return True
    else:
        return False


class SequenceIOError(Exception):
    """for errors related to sequence IO"""

    pass


class SequenceReader(abc.ABC):
    """base class for all seq reader objects"""

    def __iter__(self) -> Iterator[SequenceRecord]:
        """Yield all sequences in the file"""
        return self.iter_seqs()

    @abc.abstractmethod
    def iter_seqs(self) -> Iterator[SequenceRecord]:
        """
        Yield all sequences from all files covered by the reader

        Yields
        ------
        SequenceRecord
        """
        return self.__iter__()

    @abc.abstractmethod
    def iter_seqs_with_filenames(self) -> Iterator[tuple[Path, SequenceRecord]]:
        """
        Yield all sequences and origin file path from all files covered by the reader

        Yields
        ------
        tuple[Path, SequenceRecord]
            the source file path and the sequence record originating from it
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def get_sequence_files(self) -> tuple[Path, ...]:
        """
        Get all sequence files covered by this reader

        Returns
        -------
        tuple[Path, ...]
        """
        raise NotImplementedError()


class SingleFileSequenceReader(SequenceReader):
    """
    Read sequence data from a file on the fly

    Parameters
    ----------
    sequence_file: os.PathLike
        path to the sequence file to read

    Attributes
    ----------
    sequence_file: Path
        sequence file path
    """

    def __init__(self, sequence_file: os.PathLike):
        self.sequence_file: Path = Path(sequence_file)

        # validate file quickly
        try:
            with dnaio.open(self.sequence_file):
                pass
        except UnknownFileFormat as e:
            raise SequenceIOError(f"File {sequence_file} is not a recognized sequence file") from e

    def iter_seqs(self) -> Iterator[SequenceRecord]:
        """
        Yield all sequences from the file

        Yields
        ------
        SequenceRecord
        """
        with dnaio.open(self.sequence_file) as f:
            for record in f:
                yield record

    def iter_seqs_with_filenames(self) -> Iterator[tuple[Path, SequenceRecord]]:
        """
        Yield all sequences from the file, paired with the filename

        Yields
        ------
        tuple[Path, SequenceRecord]
        """
        with dnaio.open(self.sequence_file) as f:
            for record in f:
                yield self.sequence_file, record

    def get_sequence_files(self) -> tuple[Path, ...]:
        """
        Get all sequence files covered by this reader

        Returns
        -------
        tuple[Path, ...]
        """
        return (self.sequence_file,)


class MultiFileSequenceReader(SequenceReader):
    """
    Read all sequences from a bulk list of files

    Notes
    -----
    Files must be single fastq files (in any compression).
    They cannot be directories or glob patterns.

    Files must also end in a standard fasta/fastq file extension.

    Parameters
    ----------
    sequence_files: Sequence[os.PathLike]
        paths to the sequence files to read

    Attributes
    ----------
    sequence_files: tuple[Path, ...]
    """

    def __init__(self, sequence_files: Sequence[os.PathLike]):
        self.sequence_files: tuple[Path, ...] = tuple(Path(sequence_file) for sequence_file in sequence_files)
        self._sequence_readers: tuple[SingleFileSequenceReader, ...] = tuple(
            SingleFileSequenceReader(_file) for _file in self.sequence_files
        )

    def iter_seqs(self) -> Iterator[SequenceRecord]:
        """
        Yield all sequences from all files covered by the reader

        Yields
        ------
        SequenceRecord
        """
        for reader in self._sequence_readers:
            for record in reader:
                yield record

    def iter_seqs_with_filenames(self) -> Iterator[tuple[Path, SequenceRecord]]:
        """
        Yield all sequences and source file from all files covered by the reader

        Yields
        ------
        Path, SequenceRecord
            the source file path and the sequence record originating from it
        """
        for reader in self._sequence_readers:
            for source, record in reader.iter_seqs_with_filenames():
                yield source, record

    def get_sequence_files(self) -> tuple[Path, ...]:
        """
        Get all sequence files covered by this reader

        Returns
        -------
        tuple[Path, ...]
        """
        return self.sequence_files


class SequenceDirectoryReader(MultiFileSequenceReader):
    """
    Read all sequences from all files in a given directory

    Notes
    -----
    Only files with commonly used fasta/fastq file extensions will be detected
    and read. If you are using non-standard extensions, please use
    `MultiFileSequenceReader` instead.

    Parameters
    ----------
    sequence_dir: os.PathLike
        path to the directory with sequence files to read

    Attributes
    ----------
    sequence_dir: Path
        Path object to the directory with sequence files
    sequence_files: tuple[Path, ...]
        all sequence files detected in the directory
    """

    def __init__(self, sequence_dir: os.PathLike):
        self.sequence_dir: Path = Path(sequence_dir)

        possible_files = glob.glob(f"{self.sequence_dir}/*")

        if len(possible_files) == 0:
            raise SequenceIOError(f"No sequence files detected in directory {self.sequence_dir}")

        super().__init__(sequence_files=[Path(file) for file in possible_files if _check_if_seq_file(file)])


class SequenceGlobReader(MultiFileSequenceReader):
    """
    Read all sequences from all files from a glob pattern

    matches any valid glob pattern that python's
    `glob.glob` function can handle

    Notes
    -----
    Only files with commonly used fasta/fastq file extensions will be detected
    and read.

    Parameters
    ----------
    sequence_glob: os.PathLike
        path to the directory with sequence files to read

    Attributes
    ----------
    sequence_glob: Path
        Path object to the directory with sequence files
    sequence_files: tuple[Path, ...]
        all sequence files detected in the directory
    """

    def __init__(self, sequence_glob: os.PathLike):
        self.sequence_glob: Path = Path(sequence_glob)

        super().__init__(
            sequence_files=[Path(file) for file in glob.glob(f"{self.sequence_glob}") if _check_if_seq_file(file)]
        )


def get_reader(sequence_inputs: os.PathLike | Sequence[os.PathLike]) -> MultiFileSequenceReader:
    """
    Get a sequence reader from a list of sequence inputs

    Inputs can be any combination of files, directories or globs.

    Parameters
    ----------
    sequence_inputs: os.PathLike | Sequence[os.PathLike]
        path(s) to sequence file(s), directory(ies) or glob pattern(s)

    Returns
    -------
    MultiFileSequenceReader
        reader covering all sequences from all inputs
    """
    if not isinstance(sequence_inputs, Sequence):
        sequence_inputs = (sequence_inputs,)

    _files: list[Path] = []
    for sequence_input in sequence_inputs:
        path = Path(sequence_input)
        if path.is_file():
            # single file
            _files.append(path)
        elif path.is_dir():
            # directory
            _files.extend([Path(file) for file in glob.glob(f"{path}/*") if _check_if_seq_file(file)])
        else:
            # glob
            _files.extend([Path(file) for file in glob.glob(f"{path}") if _check_if_seq_file(file)])

    return MultiFileSequenceReader(sequence_files=_files)
