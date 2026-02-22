"""Module for tracking statistics during decoding and summarizing the results after decoding"""

import json
import os
from collections import defaultdict
from dataclasses import dataclass, field

from deli.dels.library import Library


class DecodeStatistics:
    """
    Track the statistics of the decoding run

    This includes information about number of sequences attempted, number of sequences decoded,
    failures during decoding, and number of sequences/compounds decoded per library.

    Note: only `num_seqs_decoded_per_lib` will be tracked during decoding.
    Both `num_compounds_decoded_per_lib` and `num_molecules_decoded_per_lib` can only
    be properly calculated after collection and counting, which happen separately from raw
    sequence decoding.

    Attributes
    ----------
    num_seqs_read: int
        the number of sequences read during decoding
    num_seqs_decoded_per_lib: dict[str, int]
        the number of sequences decoded per library
    num_failed_too_short: int
        the number of decoded failed because barcode read was too short
    num_failed_too_long: int
        the number of decoded failed because barcode read was too long
    num_failed_alignment: int
        the number of decoded failed because initial alignment failed
    num_failed_library_call: int
        the number of decoded failed because the library was not called
    num_failed_building_block_call: int
        the number of decoded failed because a building block was not called
    num_failed_ambiguous_building_block_call: int
        the number of decoded failed because a building block call was ambiguous
    num_failed_umi: int
        the number of decoded failed because UMI calling failed
    """

    def __init__(self):
        """Initialize a DecodeStatistics object"""
        self.num_seqs_read: int = 0
        self.num_seqs_decoded_per_lib: defaultdict[str, int] = defaultdict(int)

        # track the unique failures
        self.num_failed_too_short: int = 0
        self.num_failed_too_long: int = 0
        self.num_failed_alignment: int = 0
        self.num_failed_library_call: int = 0
        self.num_failed_building_block_call: int = 0
        self.num_failed_ambiguous_building_block_call: int = 0
        self.num_failed_umi: int = 0

    def __add__(self, other) -> "DecodeStatistics":
        """Add two DecodeStatistics objects together"""
        if not isinstance(other, DecodeStatistics):
            raise TypeError(f"unsupported operand type(s) for +: {type(self)} and {type(other)}")
        result = DecodeStatistics()
        result.num_seqs_read = self.num_seqs_read + other.num_seqs_read
        result.num_failed_too_short = self.num_failed_too_short + other.num_failed_too_short
        result.num_failed_too_long = self.num_failed_too_long + other.num_failed_too_long
        result.num_failed_library_call = self.num_failed_library_call + other.num_failed_library_call
        result.num_failed_alignment = self.num_failed_alignment + other.num_failed_alignment
        result.num_failed_building_block_call = (
            self.num_failed_building_block_call + other.num_failed_building_block_call
        )
        result.num_failed_ambiguous_building_block_call = (
            self.num_failed_ambiguous_building_block_call + other.num_failed_ambiguous_building_block_call
        )
        result.num_failed_umi = self.num_failed_umi + other.num_failed_umi

        # Merge defaultdicts
        result.num_seqs_decoded_per_lib = defaultdict(
            int,
            {
                k: self.num_seqs_decoded_per_lib[k] + other.num_seqs_decoded_per_lib[k]
                for k in set(self.num_seqs_decoded_per_lib) | set(other.num_seqs_decoded_per_lib)
            },
        )
        return result

    def __str__(self) -> str:
        """Convert the statistic object to a string (new line separated)"""
        return "\n".join([f"{key}={val}\n" for key, val in self.__dict__.items()])

    def __repr__(self) -> str:
        """Represent the statistic object as string ('; ' separated)"""
        return "; ".join([f"{key}={val}\n" for key, val in self.__dict__.items()])

    def with_libraries(self, libraries: list[Library]) -> "DecodeStatistics":
        """
        Add libraries to the statistics object for better reporting

        Will set the number of sequences decoded and compounds decoded for each library to 0 if not already present

        Parameters
        ----------
        libraries: list[Library]

        Returns
        -------
        DecodeStatistics
        """
        return self.with_library_ids([lib.library_id for lib in libraries])

    def with_library_ids(self, library_ids: list[str]) -> "DecodeStatistics":
        """
        Add libraries to the statistics object for better reporting

        Will set the number of sequences decoded and compounds decoded for each library to 0 if not already present

        Parameters
        ----------
        library_ids: list[str]

        Returns
        -------
        DecodeStatistics
        """
        result = self
        for lib_id in library_ids:
            if lib_id not in result.num_seqs_decoded_per_lib:
                result.num_seqs_decoded_per_lib[lib_id] = 0
        return result

    @property
    def num_seqs_decoded(self) -> int:
        """Number of sequences decoded successfully in total"""
        return sum(self.num_seqs_decoded_per_lib.values())

    def to_file(self, out_path: str | os.PathLike):
        """
        Write the statistics to a file

        Notes
        -----
        Will be in JSON format

        Parameters
        ----------
        out_path: str or os.PathLike
            path to write the statistics to
        """
        with open(out_path, "w") as out_file:
            json.dump(self.__dict__, out_file, indent=4)

    @classmethod
    def from_file(cls, path: str | os.PathLike) -> "DecodeStatistics":
        """
        Read in a Statistics Object from a file

        Must be in JSON format

        Parameters
        ----------
        path: str or os.PathLike
            path to read the statistics from

        Returns
        -------
        DecodeStatistics
        """
        result = cls()
        data = json.load(open(path, "r"))

        result.num_seqs_read = data.get("num_seqs_read", 0)

        result.num_seqs_decoded_per_lib = defaultdict(
            int,
            {str(key): int(val) for key, val in data.get("num_seqs_decoded_per_lib", {}).items()},
        )

        # collect error info
        result.num_failed_too_short = data.get("num_failed_too_short", 0)
        result.num_failed_too_long = data.get("num_failed_too_long", 0)
        result.num_failed_alignment = data.get("num_failed_alignment", 0)
        result.num_failed_library_call = data.get("num_failed_library_call", 0)
        result.num_failed_building_block_call = data.get("num_failed_building_block_call", 0)
        result.num_failed_ambiguous_building_block_call = data.get("num_failed_ambiguous_building_block_call", 0)
        result.num_failed_umi = data.get("num_failed_umi", 0)
        return result


@dataclass(frozen=True)
class DecodeSummary:
    """
    A summary of the number of sequences, compounds, and molecules decoded in a selection

    This data is often needed for computing various enrichment metrics.
    This data is also frozen, as it is not safe to modify the values contained within
    do to the loss of information about specific compounds/molecules when aggregating
    the summary info.

    Notes
    -----
    Summary objects should only be created by reading a file
    (ideally created by the `deli decode summarize` command).

    Attributes
    ----------
    total_seqs_read: int
        the total number of sequences read during decoding
    total_seqs_decoded: int
        the total number of sequences decoded in the selection
    seqs_decoded_per_lib: dict[str, int]
        the number of sequences decoded per library in the selection
    total_compounds_decoded: int
        the total number of compounds decoded in the selection
    compounds_decoded_per_lib: dict[str, int]
        the number of compounds decoded per library in the selection
    total_molecules_decoded: int
        the total number of molecules decoded in the selection
    molecules_decoded_per_lib: dict[str, int]
        the number of molecules decoded per library in the selection
    """

    num_seqs_read: int
    num_seqs_decoded: int = field(init=False)
    seqs_decoded_per_lib: dict[str, int]
    num_compounds_decoded: int = field(init=False)
    compounds_decoded_per_lib: dict[str, int]
    num_molecules_decoded: int = field(init=False)
    molecules_decoded_per_lib: dict[str, int]

    def __post_init__(self):
        """Precompute the calculation of total sequences, compounds, and molecules decoded for easier access later"""
        object.__setattr__(self, "num_seqs_decoded", sum(self.seqs_decoded_per_lib.values()))
        object.__setattr__(self, "num_compounds_decoded", sum(self.compounds_decoded_per_lib.values()))
        object.__setattr__(self, "num_molecules_decoded", sum(self.molecules_decoded_per_lib.values()))

    def with_libraries(self, libraries: list[Library]) -> "DecodeSummary":
        """
        Add libraries to the summary object for better reporting

        Will set the number of sequences decoded and compounds decoded for each library to 0 if not already present

        Parameters
        ----------
        libraries: list[Library]

        Returns
        -------
        DecodeSummary
        """
        return self.with_library_ids([lib.library_id for lib in libraries])

    def with_library_ids(self, library_ids: list[str]) -> "DecodeSummary":
        """
        Add libraries to the summary object for better reporting

        Will set the number of sequences decoded and compounds decoded for each library to 0 if not already present

        Parameters
        ----------
        library_ids: list[str]

        Returns
        -------
        DecodeSummary
        """
        result = self
        for lib_id in library_ids:
            if lib_id not in result.seqs_decoded_per_lib:
                result.seqs_decoded_per_lib[lib_id] = 0
            if lib_id not in result.compounds_decoded_per_lib:
                result.compounds_decoded_per_lib[lib_id] = 0
            if lib_id not in result.molecules_decoded_per_lib:
                result.molecules_decoded_per_lib[lib_id] = 0
        return result

    def to_file(self, out_path: str | os.PathLike):
        """
        Write the statistics to a file

        Notes
        -----
        Will be in JSON format

        Parameters
        ----------
        out_path: str or os.PathLike
            path to write the statistics to
        """
        with open(out_path, "w") as out_file:
            json.dump(self.__dict__, out_file, indent=4)

    @classmethod
    def from_file(cls, path: str | os.PathLike) -> "DecodeSummary":
        """
        Read in a Summary Object from a file

        Must be in JSON format

        Parameters
        ----------
        path: str or os.PathLike
            path to read the summary from

        Returns
        -------
        DecodeSummary
        """
        data = json.load(open(path, "r"))

        num_seqs_read = data["num_seqs_read"]

        num_seqs_decoded_per_lib = defaultdict(
            int,
            {str(key): int(val) for key, val in data["num_seqs_decoded_per_lib"].items()},
        )

        num_compounds_decoded_per_lib = defaultdict(
            int,
            {str(key): int(val) for key, val in data["num_compounds_decoded_per_lib"].items()},
        )

        num_molecules_decoded_per_lib = defaultdict(
            int,
            {str(key): int(val) for key, val in data["num_molecules_decoded_per_lib"].items()},
        )

        return cls(
            num_seqs_read=num_seqs_read,
            seqs_decoded_per_lib=num_seqs_decoded_per_lib,
            compounds_decoded_per_lib=num_compounds_decoded_per_lib,
            molecules_decoded_per_lib=num_molecules_decoded_per_lib,
        )
