"""decoding statistics tracking"""

import json
import os
from collections import defaultdict
from typing import Self


class Statistics:
    """Base class for all statistic classes"""

    def __str__(self) -> str:
        """Convert the statistic object to a string (new line seperated)"""
        return "\n".join([f"{key}={val}\n" for key, val in self.__dict__.items()])

    def __repr__(self) -> str:
        """Represent the statistic object as string ('; ' seperated)"""
        return "; ".join([f"{key}={val}\n" for key, val in self.__dict__.items()])

    def to_file(self, out_path: str | os.PathLike):
        """
        Write the statistics to a file

        Will be in JSON format

        Parameters
        ----------
        out_path: str or os.PathLike
            path to write the statistics to
        """
        json.dump(self.__dict__, open(out_path, "w"))

    @classmethod
    def from_file(cls, path: str | os.PathLike) -> Self:
        """
        Read in a Statistics Object from a file

        Must be in JSON format

        Parameters
        ----------
        path: str or os.PathLike
            path to read the statistics from

        Returns
        -------
        Self
        """
        return cls(**json.load(open(path, "r")))


class DecodeStatistics(Statistics):
    """
    Track the statistics of the decoding run

    Will count how many sequences are read in,
    how many are decoded, and how many remain after degen.
    Will split decode/degen by library.
    It Will also track the number of times decoding failed,
    and for what reasons

    Attributes
    ----------
    num_seqs_read: int
        the number of sequences read during decoding
    num_seqs_decoded_per_lib: defaultdict[str, int]
        the number of sequences decoded per library
    num_seqs_degen_per_lib: defaultdict[str, int]
        the number of sequences degened per library
    num_failed_too_short: int
        the number of decoded failed because sequence read was too short
    num_failed_too_long: int
        the number of decoded failed because sequence read was too long
    num_failed_library_call: int
        the number of decoded failed because the library was not called
    num_failed_library_match_too_short: int
        the number of decoded failed because the library match was too short
    num_failed_building_block_call: int
        the number of decoded failed because a building block was not called
    num_failed_alignment: int
        the number of decoded failed because alignment failed
    """

    def __init__(self):
        """Initialize a DecodeStatistics object"""
        self.num_seqs_read: int = 0
        self.num_seqs_decoded_per_lib: defaultdict[str, int] = defaultdict(int)
        self.num_seqs_degen_per_lib: defaultdict[str, int] = defaultdict(int)

        # track the unique failures
        self.num_failed_too_short: int = 0
        self.num_failed_too_long: int = 0
        self.num_failed_library_call: int = 0
        self.num_failed_library_match_too_short: int = 0
        self.num_failed_building_block_call: int = 0
        self.num_failed_alignment: int = 0
