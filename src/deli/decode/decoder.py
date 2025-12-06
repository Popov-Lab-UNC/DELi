"""code for calling DNA barcodes"""

import abc
import dataclasses
import json
import os
from collections import defaultdict
from dataclasses import asdict
from typing import Any, Generic, Iterator, Optional, TypeVar, Literal

from dnaio import SequenceRecord
from tqdm import tqdm

from deli.dels.barcode import BarcodeSection
from deli.dels.base import Library
from deli.dels.building_block import TaggedBuildingBlock
from deli.dels.combinatorial import DELibrary
from deli.dels.compound import DELCompound
from deli.dels.tool_compounds import DopedToolCompound, TaggedToolCompound, TaggedToolCompoundLibrary, ToolCompound
from deli.enumeration.enumerator import EnumerationRunError
from deli.selection import DELSelection

from .base import FailedDecodeAttempt
from .barcode_calling import AmbiguousBarcodeCall, BarcodeCaller, FailedBarcodeLookup, ValidCall, get_barcode_caller
from .library_demultiplex import AlignedSeq, FailedStaticAlignment, LibraryDemultiplexer, get_library_demultiplexer_type
from .umi import UMI
from ..dna import SequenceReader

MAX_RETRIES = 50  # number of alignments to attempt before giving up


class DecodeStatistics:
    """
    Track the statistics of the decoding run

    Will count how many sequences are read in,
    how many are decoded.
    It Will also track the number of times decoding failed,
    and for what reasons.

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

    @property
    def num_seqs_decoded(self) -> int:
        """Number of sequences decoded successfully in total"""
        return sum(self.num_seqs_decoded_per_lib.values())

    def to_file(self, out_path: str | os.PathLike, include_read_lengths: bool = False):
        """
        Write the statistics to a file

        Notes
        -----
        Will be in JSON format

        Parameters
        ----------
        out_path: str or os.PathLike
            path to write the statistics to
        include_read_lengths: bool, default = False
            if True, will include the read lengths in the file
        """
        if include_read_lengths:
            json.dump(self.__dict__, open(out_path, "w"))
        else:
            _dict = self.__dict__.copy()
            del _dict["seq_lengths"]
            json.dump(_dict, open(out_path, "w"))

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


class DecodedCompound(abc.ABC):
    """Base class for compounds that were decoded from a DNA read"""

    umi: UMI | None

    def has_umi(self) -> bool:
        """Check if the decoded compound has a UMI"""
        return self.umi is not None

    @abc.abstractmethod
    def to_cube_row_dict(self) -> dict[str, str]:
        """
        Convert the decoded compound to a dictionary of strings for writing to a cube

        Dictionary keys can be one of the following:
        - "library_id":
            the library ID
        - "compound_id":
            the compound ID
        - "smiles":
            the enumerated SMILES of the compound
        - "BB##_id":
            the id for the building block from cycle ##
        - "BB##_smiles":
            the smiles for the building block from cycle ##

        When writing to a cube, if any key is missing it will be filed with "null"
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def get_smiles(self) -> str:
        """
        Get the SMILES of the decoded compound

        Note
        ----
        Should enumerate the compound if able and needed

        Returns
        -------
        str
            the enumerated SMILES of the compound

        Raises
        ------
        EnumerationRunError
            if enumeration fails
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def get_library_id(self) -> str:
        """
        Get the library ID of the decoded compound

        Returns
        -------
        str
            the library ID
        """
        raise NotImplementedError()


class DecodedDELCompound(DELCompound, DecodedCompound):
    """
    Holds information about a decoded barcode

    Can handle both successful and failed decodes

    Notes
    -----
    IMPORTANT: two Decoded Compounds are considered identical if their
    DEL_IDs are the same, even if they have different UMIs.
    This might not seem intuitive, but is needed to help with
    degeneration (based on UMI) that happens later

    Parameters
    ----------
    library_call : ValidCall[DELibrary]
        the called library
    building_block_calls : list[ValidCall[TaggedBuildingBlock]]
        the building block calls *in order of the library's building block sections*
    umi: UMI or `None`
        the UMI for the read
        if not using umi or no umi in the barcode, use a `None`

    Attributes
    ----------
    compound_id: str
        the unique ID of the decoded compound
        composed of the library ID and building block IDs
    library: DELibrary
        the library associated with the decoded compound
    building_blocks: list[TaggedBuildingBlock]
        the building blocks associated with the decoded compound
    """

    def __init__(
        self,
        library_call: ValidCall[DELibrary],
        building_block_calls: list[ValidCall[TaggedBuildingBlock]],
        umi: UMI | None = None,
    ):
        self.library_call = library_call
        self.building_block_calls = building_block_calls
        self.umi = umi
        super().__init__(library=library_call.obj, building_blocks=[bb_call.obj for bb_call in building_block_calls])

    def to_cube_row_dict(self, enumerate_compounds: bool = False) -> dict[str, str]:
        """
        Convert the decoded DEL compound to a dictionary of strings for writing to a cube

        Dictionary keys will be:
        - "library_id":
            the library ID
        - "compound_id":
            the unique DEL compound ID
        - "smiles":
            the enumerated SMILES of the compound
            only if `enumerate_compounds` is True
        - "BB##_id":
            the id for the building block from cycle ##
        - "BB##_smiles":
            the smiles for the building block from cycle ##

        There will be one "BB##_id" and "BB##_smiles" entry for each building block in the compound

        Parameters
        ----------
        enumerate_compounds: bool, default = False
            if True, will attempt to enumerate SMILES of the compound
            and add to output dictionary with key "smiles"

        Returns
        -------
        dict[str, str]
        """
        row_dict: dict[str, str] = {
            "library_id": self.library.library_id,
            "compound_id": self.compound_id,
        }

        if enumerate_compounds:
            try:
                row_dict["smiles"] = self.library.enumerator.enumerate_by_bbs(self.building_blocks)
            except EnumerationRunError:
                row_dict["smiles"] = "ENUMERATION_FAILED"

        for idx, bb in enumerate(self.building_blocks, start=1):
            row_dict[f"BB{idx:02d}_id"] = bb.bb_id
            row_dict[f"BB{idx:02d}_smiles"] = bb.smi if bb.has_smiles() else "null"
        return row_dict

    def get_smiles(self) -> str:
        """
        Get the SMILES of the decoded compound

        Returns
        -------
        str
            the enumerated SMILES of the compound

        Raises
        ------
        EnumerationRunError
            if enumeration fails
        """
        return self.enumerate().smi

    def get_library_id(self) -> str:
        """
        Get the library ID of the decoded compound

        Returns
        -------
        str
            the library ID
        """
        return self.library.library_id


class DecodedToolCompound(DecodedCompound):
    """
    Holds information about a decoded barcode for a tool compound

    This covers both tool compounds added to selections *and*
    compounds doped into DEL libraries.

    Parameters
    ----------
    tool_library_call: ValidCall[TaggedToolCompoundLibrary] | ValidCall[DELibrary]
        the decoded tool compound library
    tool_compound_call: ValidCall[TaggedToolCompound] | ValidCall[DopedToolCompound]
        the decoded tool compound
    """

    def __init__(
        self,
        tool_library_call: ValidCall[TaggedToolCompoundLibrary] | ValidCall[DELibrary],
        tool_compound_call: ValidCall[TaggedToolCompound] | ValidCall[DopedToolCompound],
    ):
        self.tool_library_call = tool_library_call
        self.tool_compound_call = tool_compound_call
        self.tool_compound: ToolCompound = tool_compound_call.obj

    def to_cube_row_dict(self) -> dict[str, str]:
        """
        Convert the decoded tool compound to a dictionary of strings for writing to a cube

        Dictionary keys will be:
        - "library_id":
            the library ID
        - "compound_id":
            the unique DEL compound ID
        - "smiles":
            the SMILES of the compound (if present)

        Returns
        -------
        dict[str, str]
        """
        if self.tool_compound.has_smiles():
            return {
                "library_id": "ToolCompound",
                "compound_id": self.tool_compound.compound_id,
                "smiles": self.tool_compound.smi,
            }
        else:
            return {
                "library_id": "ToolCompound",
                "compound_id": self.tool_compound.compound_id,
            }

    def get_smiles(self) -> str:
        """
        Get the SMILES of the decoded tool compound

        Returns
        -------
        str
            the SMILES of the compound

        Raises
        ------
        ValueError
            if the tool compound has no SMILES
        """
        if self.tool_compound.has_smiles():
            return self.tool_compound.smi
        else:
            raise ValueError(f"Tool compound {self.tool_compound.compound_id} has no SMILES")

    def get_library_id(self) -> str:
        """
        Get the library ID of the decoded compound

        Returns
        -------
        str
            the library ID
        """
        return self.tool_compound.compound_id


class ReadTooShort(FailedDecodeAttempt):
    """returned if read to decode was too short (based on settings)"""

    def __init__(self, sequence):
        super().__init__(sequence=sequence, reason="Read too short")


class ReadTooLong(FailedDecodeAttempt):
    """returned if read to decode was too long (based on settings)"""

    def __init__(self, sequence):
        super().__init__(sequence=sequence, reason="Read too long")


class FailedBuildingBlockCall(FailedDecodeAttempt):
    """returned if building block call failed"""

    def __init__(self, sequence, barcode: str, bb_section_name: str):
        super().__init__(
            sequence=sequence,
            reason=f"failed to call building block barcode '{barcode}' for section '{bb_section_name}'",
        )


class AmbiguousBuildingBlockBarcode(FailedDecodeAttempt):
    """returned if building block barcode is ambiguous"""

    def __init__(self, sequence, barcode: str, bb_section_name: str):
        super().__init__(
            sequence=sequence,
            reason=f"ambiguous barcode for building block barcode '{barcode}' in section '{bb_section_name}'",
        )


class ImpossibleBuildingBlockBarcode(FailedDecodeAttempt):
    """returned if building block barcode is not covered by the library"""

    def __init__(self, sequence, bb_ids: list[str], library_id: str):
        super().__init__(
            sequence=sequence,
            reason=f"set of called bb_ids '{bb_ids}' is not covered by library '{library_id}'",
        )


class AlignmentFailed(FailedDecodeAttempt):
    """returned if alignment appears to be poor"""

    def __init__(self, sequence):
        super().__init__(sequence=sequence, reason="alignment of read to barcode is poor")


class DELCollectionDecoder:
    """
    Decodes reads into DEL compounds

    Parameters
    ----------
    library_demultiplexer: LibraryDemultiplexer
        the library demultiplexer to use for calling libraries
    wiggle: bool, default = False
        If true, allow for wiggling the tag to find the best match
    max_read_length: int or None, default = None
        maximum length of a read to be considered for decoding
        if above the max, decoding will fail
        if `None` will default to 5x the min_read_length
    min_read_length: int or None, default = None
        minimum length of a read to be considered for decoding
        if below the min, decoding will fail
        if `None` will default to the smallest min match length of
        any library in the collection considered for decoding
    default_error_correction_mode_str: str, default = "levenshtein_dist:1,asymmetrical"
        If a decodable section does not have an error correction mode defined,
        use this one by default.
    decode_statistics: DecodeStatistics or None, default = None
        the statistic tracker for the decoding run
        if None, will initialize a new, empty statistic object

    Attributes
    ----------
    library_decoders: dict[DELibrary, DELibraryDecoder]
        the library decoders for each library in the demultiplexer
    """

    def __init__(
        self,
        library_demultiplexer: LibraryDemultiplexer,
        wiggle: bool = False,
        max_read_length: Optional[int] = None,
        min_read_length: Optional[int] = None,
        default_error_correction_mode_str: str = "levenshtein_dist:1,asymmetrical",
        decode_statistics: Optional[DecodeStatistics] = None,
    ):
        self.library_demultiplexer = library_demultiplexer
        self.decode_statistics: DecodeStatistics = (
            decode_statistics if decode_statistics is not None else DecodeStatistics()
        )

        # build library decoders
        self.library_decoders: dict[Library, LibraryDecoder[Any]] = {
            library: DELibraryDecoder(
                library=library, wiggle=wiggle, default_error_correction_mode_str=default_error_correction_mode_str
            )
            if isinstance(library, DELibrary)
            else ToolCompoundDecoder(
                library=library, wiggle=wiggle, default_error_correction_mode_str=default_error_correction_mode_str
            )
            for library in self.library_demultiplexer.all_libraries
        }

        # set the min/max lengths
        self._min_read_length: int
        if isinstance(min_read_length, int):
            self._min_read_length = min_read_length
        else:
            self._min_read_length = (
                min([lib.barcode_schema.min_length for lib in library_demultiplexer.libraries]) - 10
            )  # allow some slack for seq errors

        self._max_read_length: int
        if isinstance(max_read_length, int):
            self._max_read_length = max_read_length
        else:
            self._max_read_length = 5 * self._min_read_length

    def decode_read(self, sequence: SequenceRecord) -> DecodedCompound | FailedDecodeAttempt:
        """
        Given a raw read, decode it to determine itself compound

        Notes
        -----
        If decoding fails, will return a FailedDecodeAttempt object
        explaining why it failed

        Parameters
        ----------
        sequence: SequenceRecord
            the raw read to decode

        Returns
        -------
        DecodedCompound or FailedDecodeAttempt
        """
        # check lengths
        if len(sequence) > self._max_read_length:
            self.decode_statistics.num_failed_too_long += 1
            return ReadTooLong(sequence)
        if len(sequence) < self._min_read_length:
            self.decode_statistics.num_failed_too_short += 1
            return ReadTooShort(sequence)

        _lib_call = self.library_demultiplexer.demultiplex(sequence)
        if isinstance(_lib_call, FailedDecodeAttempt):
            if isinstance(_lib_call, FailedStaticAlignment):
                self.decode_statistics.num_failed_alignment += 1
            else:
                self.decode_statistics.num_failed_library_call += 1
            return _lib_call
        else:
            library_call, alignment_iter = _lib_call
            library_decoder = self.library_decoders[library_call.obj]

            try:
                decoded_compound = self._attempt_to_decode_compound(library_call, library_decoder, alignment_iter)
            except BarcodeDecodingError:
                self.decode_statistics.num_failed_alignment += 1
                return FailedDecodeAttempt(sequence, "Decoding was not attempted")

            if isinstance(decoded_compound, FailedDecodeAttempt):
                if isinstance(decoded_compound, FailedBuildingBlockCall):
                    self.decode_statistics.num_failed_building_block_call += 1
                elif isinstance(decoded_compound, AmbiguousBuildingBlockBarcode):
                    self.decode_statistics.num_failed_ambiguous_building_block_call += 1
                elif isinstance(decoded_compound, UMIContainsAmbiguity):
                    self.decode_statistics.num_failed_umi += 1
                elif isinstance(decoded_compound, ImpossibleBuildingBlockBarcode):
                    self.decode_statistics.num_failed_building_block_call += 1
            return decoded_compound

    @staticmethod
    def _attempt_to_decode_compound(
        library_call: ValidCall[Any],
        library_decoder: "LibraryDecoder[Any]",
        alignment_iter: Iterator[tuple[AlignedSeq, float]],
    ) -> FailedDecodeAttempt | DecodedCompound:
        """
        Attempt to decode a compound using the provided library decoder and alignment iterator

        Parameters
        ----------
        library_call: ValidCall[Library]
            the called library
        library_decoder: LibraryDecoder
            the library decoder to use for decoding
        alignment_iter: Iterator[tuple[AlignedSeq, float]]
            the iterator of aligned sequences and their scores

        Returns
        -------
        DecodedCompound or FailedDecodeAttempt
        """
        # this is a placeholder to satisfy type checkers
        decoded_compound: None | DecodedCompound | FailedDecodeAttempt = None
        for attempt, (seq_alignment, _) in enumerate(alignment_iter):
            if attempt >= MAX_RETRIES:
                break
            decoded_compound = library_decoder.decode_sequence(library_call, seq_alignment)
            if isinstance(decoded_compound, DecodedCompound):
                return decoded_compound
        if decoded_compound is None:
            raise BarcodeDecodingError("No decoding attempts were made")
        return decoded_compound


class BarcodeDecodingError(Exception):
    """raised when there is an error during barcode decoding"""

    pass


class UMIContainsAmbiguity(FailedDecodeAttempt):
    """Returned if UMI contains ambiguous bases ('N')"""

    def __init__(self, sequence: SequenceRecord, umi_sequence: str):
        super().__init__(sequence=sequence, reason=f"UMI contains ambiguous bases: {umi_sequence}")


L = TypeVar("L", DELibrary, TaggedToolCompoundLibrary)
OBJ = TypeVar("OBJ")


class LibraryDecoder(abc.ABC, Generic[L]):
    """Base class for a library decoder"""

    def __init__(self, library: L, wiggle: bool = False, default_error_correction_mode_str="disable"):
        self.library: L = library
        self.wiggle = wiggle
        self.default_error_correction_mode_str = default_error_correction_mode_str

        self._callers = self._get_callers()

        # this distance is always defined as the distance between the last index of the
        # closest bb section and the start of the umi section
        self._dist_from_closest_sec_to_umi: int = 0
        self._direction_from_closest_sec_to_umi: int = 0
        self._umi_length: int = -1
        self._closest_sec_idx: int = -1
        self._closest_sec_stop_idx: int = -1
        self._closest_sec_start_idx: int = -1
        if self.library.barcode_schema.has_umi():
            umi_section = self.library.barcode_schema.get_section("umi")

            smallest_dist: int = 100000  # large number
            for i, bb_section in enumerate(self._callers.keys()):
                _dist = self.library.barcode_schema.get_length_between_sections(
                    bb_section.section_name,
                    "umi",
                    include_direction=True,
                )
                if abs(_dist) < smallest_dist:
                    smallest_dist = abs(_dist)
                    direction = self.library.barcode_schema.get_direction_of_sections(
                        bb_section.section_name,
                        "umi",
                    )

                    if direction == -1:
                        self._dist_from_closest_sec_to_umi = _dist - umi_section.get_overhang_length()
                    elif direction == 1:
                        self._dist_from_closest_sec_to_umi = _dist - bb_section.get_overhang_length()
                    else:
                        raise RuntimeError(
                            "direction between building block and UMI sections is always defined; "
                            "this should never happen please report a bug"
                        )

                    self._direction_from_closest_sec_to_umi = direction
                    self._closest_sec_idx = i
            self._umi_length = len(umi_section)

    @abc.abstractmethod
    def _get_callers(self) -> dict[BarcodeSection, BarcodeCaller]:
        """
        Get the barcode callers for each decodable section
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def decode_sequence(
        self, library_call: ValidCall[L], aligned_sequence: AlignedSeq
    ) -> DecodedCompound | FailedDecodeAttempt:
        """
        Given a library call, decode the barcode sections

        Parameters
        ----------
        library_call: ValidCall[L]
            the called library
        aligned_sequence: AlignedSeq
            the aligned sequence to decode the building blocks from

        Returns
        -------
        DecodedCompound or FailedDecodeAttempt
            the decoded compound if successful,
            else a FailedDecodeAttempt explaining why it failed
        """
        raise NotImplementedError()

    def _call_decodable_sections(
        self, aligned_sequence: AlignedSeq, section: BarcodeSection, caller: BarcodeCaller[OBJ]
    ) -> tuple[ValidCall[OBJ], int, int] | FailedDecodeAttempt:
        """
        Call the decodable sections from the aligned sequence

        Parameters
        ----------
        aligned_sequence: AlignedSeq
            the aligned sequence to call the building blocks from
        section: BarcodeSection
            the section to call
        caller: BarcodeCaller[OBJ]
            the barcode caller to use for calling the section

        Returns
        -------
        tuple[ValidCall[OBJ], int, int] or FailedDecodeAttempt
            the called object and the start/stop indices if successful,
            else a FailedDecodeAttempt explaining why it failed
        """
        section_codon = aligned_sequence.get_section_barcode(section.section_name)
        best_start_idx = start = aligned_sequence.section_spans[section.section_name][0]
        best_stop_idx = stop = aligned_sequence.section_spans[section.section_name][1]

        call: ValidCall[OBJ] | FailedBarcodeLookup = caller.decode_barcode(section_codon)

        if isinstance(call, FailedBarcodeLookup) and self.wiggle:
            true_section_length = len(section.section_tag)
            best_score = 100.0

            # default to failed lookup
            best_call: ValidCall[OBJ] | FailedBarcodeLookup = FailedBarcodeLookup("placeholder")

            adjusted_section_codon = aligned_sequence.sequence.sequence[start - 1 : stop + 1]

            for length in range(true_section_length, true_section_length + 1, true_section_length - 1):
                wiggle_start = 0
                wiggle_end = wiggle_start + length

                # wiggle out to +2 bases longer than true length
                while wiggle_end <= true_section_length + 2:
                    wiggle_codon = adjusted_section_codon[wiggle_start:wiggle_end]
                    wiggle_call = caller.decode_barcode(wiggle_codon)

                    if not isinstance(wiggle_call, FailedBarcodeLookup):
                        # TODO could have an option to first match instead of best match
                        if wiggle_call.score == 0:  # found perfect match quit now
                            call = wiggle_call
                            best_start_idx = start - true_section_length + wiggle_start
                            best_stop_idx = stop - true_section_length + wiggle_end
                            break
                        elif wiggle_call.score < best_score:  # better match found
                            best_score = wiggle_call.score
                            best_call = wiggle_call
                            best_start_idx = start - true_section_length + wiggle_start
                            best_stop_idx = stop - true_section_length + wiggle_end
                        elif wiggle_call.score == best_score:  # ambiguous match check
                            # if the barcodes are different but have the same score, it is ambiguous
                            if (not isinstance(best_call, FailedBarcodeLookup)) and (wiggle_call.obj != best_call.obj):
                                best_call = AmbiguousBarcodeCall(wiggle_codon)

                    # move the wiggle window
                    wiggle_start += 1
                    wiggle_end = wiggle_start + length

                if not isinstance(call, FailedBarcodeLookup):
                    break  # found perfect match quit now, else keep wiggling

            # wiggle found a non-perfect call
            if best_call is not None:
                call = best_call

        # handle errors
        if isinstance(call, FailedBarcodeLookup):
            if isinstance(call, AmbiguousBarcodeCall):
                return AmbiguousBuildingBlockBarcode(aligned_sequence.sequence, section_codon, section.section_name)
            return FailedBuildingBlockCall(aligned_sequence.sequence, section_codon, section.section_name)
        else:
            return call, best_start_idx, best_stop_idx

    def call_umi(self, aligned_sequence: AlignedSeq) -> UMI | None | FailedDecodeAttempt:
        """
        Call the UMI from the aligned sequence

        Parameters
        ----------
        aligned_sequence: AlignedSeq
            the aligned sequence to call the UMI from

        Returns
        -------
        UMI or None or FailedDecodeAttempt
            the UMI if successful,
            None if no UMI section in library,
            else a FailedDecodeAttempt explaining why it failed
        """
        if self._umi_length == -1:
            return None
        elif self._closest_sec_start_idx == -1:
            raise BarcodeDecodingError("Cannot call UMI before building blocks have been called")
        else:
            if self._direction_from_closest_sec_to_umi == -1:
                umi_start = self._closest_sec_start_idx + self._dist_from_closest_sec_to_umi
            else:
                umi_start = self._closest_sec_stop_idx + self._dist_from_closest_sec_to_umi
            umi_stop = umi_start + self._umi_length
            umi_codon = aligned_sequence.sequence.sequence[umi_start:umi_stop]
            aligned_umi_codon = aligned_sequence.get_section_barcode("umi")
            if ("N" in umi_codon) or ("N" in aligned_umi_codon):
                return UMIContainsAmbiguity(aligned_sequence.sequence, umi_codon)
            if umi_codon != aligned_umi_codon:
                return UMI([umi_codon, aligned_umi_codon])
            else:
                return UMI([umi_codon])


class DELibraryDecoder(LibraryDecoder[DELibrary]):
    """
    Decodes the barcodes for variable sections in a DEL

    Currently only support building block sections

    Parameters
    ----------
    library: DELibrary
        the library to decode from
    wiggle: bool, default = False
        If true, allow for wiggling the tag to find the best match
    default_error_correction_mode_str: str, default = "levenshtein_dist:1,asymmetrical"
        If a decodable section does not have an error correction mode defined,
        use this one by default.
    """

    def __init__(
        self,
        library: DELibrary,
        wiggle: bool = False,
        default_error_correction_mode_str: str = "levenshtein_dist:1,asymmetrical",
    ):
        self._has_doped = len(library.tool_compounds) > 0
        # make "fake" building blocks for the doped compounds
        # ids are equal to the tool compound ids
        self._doped_building_blocks: list[list[TaggedBuildingBlock]] = list()
        self._doped_id_map: dict[str, DopedToolCompound] = {
            tool_comp.compound_id: tool_comp for tool_comp in library.tool_compounds
        }
        for i in range(library.num_cycles):
            cycle_bbs: list[TaggedBuildingBlock] = list()
            for tool_comp in library.tool_compounds:
                cycle_bbs.append(
                    TaggedBuildingBlock(
                        bb_id=tool_comp.compound_id,
                        tag=tool_comp.bb_tags[i],
                    )
                )
            self._doped_building_blocks.append(cycle_bbs)
        super().__init__(library, wiggle, default_error_correction_mode_str)

    def _get_callers(self) -> dict[BarcodeSection, BarcodeCaller]:
        callers: dict[BarcodeSection, BarcodeCaller[TaggedBuildingBlock]] = {}
        for i, (bb_sec, bb_set) in enumerate(self.library.iter_bb_barcode_sections_and_sets()):
            error_correction_mode_str = getattr(bb_sec, "error_correction_mode", self.default_error_correction_mode_str)
            callers[bb_sec] = get_barcode_caller(
                # add in the doped building blocks to the caller
                tag_map={
                    tag: bb for bb in list(bb_set.building_blocks) + self._doped_building_blocks[i] for tag in bb.tags
                },
                error_correction_mode_str=error_correction_mode_str,
            )
        return callers

    def decode_sequence(
        self, library_call: ValidCall[DELibrary], aligned_sequence: AlignedSeq
    ) -> DecodedCompound | FailedDecodeAttempt:
        """
        Given a library call, decode the barcode

        Parameters
        ----------
        library_call: ValidCall[DELibrary]
            the called DEL
        aligned_sequence: AlignedSeq
            the aligned sequence to decode the building blocks from

        Returns
        -------
        DecodedCompound or FailedDecodeAttempt
            the decoded compound if successful,
            else a FailedDecodeAttempt explaining why it failed
        """
        bb_calls: list[ValidCall[TaggedBuildingBlock]] = []

        for i, (bb_section, bb_caller) in enumerate(self._callers.items()):
            bb_call_or_fail = self._call_decodable_sections(aligned_sequence, bb_section, bb_caller)
            if isinstance(bb_call_or_fail, FailedDecodeAttempt):
                return bb_call_or_fail
            else:
                bb_call, start_idx, stop_idx = bb_call_or_fail
                bb_calls.append(bb_call)
                if i == self._closest_sec_idx:
                    self._closest_sec_start_idx = start_idx
                    self._closest_sec_stop_idx = stop_idx

        umi_call = self.call_umi(aligned_sequence)
        if isinstance(umi_call, FailedDecodeAttempt):
            return umi_call
        else:
            if self._has_doped and any([bb_call.obj.bb_id in self._doped_id_map.keys() for bb_call in bb_calls]):
                bb_id_set = {bb_call.obj.bb_id for bb_call in bb_calls}
                if len(bb_id_set) > 1:
                    return ImpossibleBuildingBlockBarcode(
                        aligned_sequence.sequence,
                        list(bb_id_set),
                        library_call.obj.library_id,
                    )
                else:
                    try:
                        doped_compound = self._doped_id_map[list(bb_id_set)[0]]
                    except Exception as e:
                        raise BarcodeDecodingError(
                            f"cannot find doped compound for bb_id '{bb_id_set}'; something went wrong"
                        ) from e
                    tool_call: ValidCall[DopedToolCompound] = ValidCall(
                        doped_compound, sum([bb_call.score for bb_call in bb_calls])
                    )
                    return DecodedToolCompound(
                        tool_library_call=library_call,
                        tool_compound_call=tool_call,
                    )
            return DecodedDELCompound(
                library_call=library_call,
                building_block_calls=bb_calls,
                umi=umi_call,
            )


class ToolCompoundDecoder(LibraryDecoder[TaggedToolCompoundLibrary]):
    """
    Decodes the barcode adn UMI for a tool compound

    Parameters
    ----------
    library: TaggedToolCompoundLibrary
        the tool compound library to decode from
    wiggle: bool, default = False
        If true, allow for wiggling the tag to find the best match
    default_error_correction_mode_str: str, default = "levenshtein_dist:1,asymmetrical"
        If a decodable section does not have an error correction mode defined,
        use this one by default.

    Attributes
    ----------
    tool_caller: BarcodeCaller[TaggedToolCompound]
        the barcode caller for the tool compound
    """

    def __init__(
        self,
        library: TaggedToolCompoundLibrary,
        wiggle: bool = False,
        default_error_correction_mode_str: str = "levenshtein_dist:1,asymmetrical",
    ):
        super().__init__(library, wiggle, default_error_correction_mode_str)

        self._tool_caller = list(self._callers.values())[0]
        self._barcode_section = self.library.barcode_schema.get_section("compound_tag")

    def _get_callers(self) -> dict[BarcodeSection, BarcodeCaller]:
        tool_comp_section = self.library.barcode_schema.get_section("compound_tag")
        error_correction_mode_str = getattr(
            tool_comp_section,
            "error_correction_mode",
            self.default_error_correction_mode_str,
        )
        tool_caller: BarcodeCaller[TaggedToolCompound] = get_barcode_caller(
            tag_map={tool_comp.tag: tool_comp for tool_comp in self.library.compounds},
            error_correction_mode_str=error_correction_mode_str,
        )
        return {tool_comp_section: tool_caller}

    def decode_sequence(
        self, library_call: ValidCall[TaggedToolCompoundLibrary], aligned_sequence: AlignedSeq
    ) -> DecodedCompound | FailedDecodeAttempt:
        """
        Given a library call, decode the barcode

        Parameters
        ----------
        library_call: ValidCall[TaggedToolCompoundLibrary]
            the called tool compound library
        aligned_sequence: AlignedSeq
            the aligned sequence to decode the building blocks from

        Returns
        -------
        DecodedCompound or FailedDecodeAttempt
            the decoded compound if successful,
            else a FailedDecodeAttempt explaining why it failed
        """
        tool_comp_call_or_fail = self._call_decodable_sections(
            aligned_sequence, self._barcode_section, self._tool_caller
        )
        if isinstance(tool_comp_call_or_fail, FailedDecodeAttempt):
            return tool_comp_call_or_fail
        else:
            tool_comp_call, start_idx, stop_idx = tool_comp_call_or_fail
            self._closest_sec_start_idx = start_idx
            self._closest_sec_stop_idx = stop_idx

        umi_call = self.call_umi(aligned_sequence)
        if isinstance(umi_call, FailedDecodeAttempt):
            return umi_call
        else:
            return DecodedToolCompound(tool_library_call=library_call, tool_compound_call=tool_comp_call)


@dataclasses.dataclass(frozen=True)
class DecodingSettings:
    """
    Define parameters for decoding experiments

    More details about the exact effect of these settings can
    be found in the "Decoding" docs

    Notes
    -----
    Only parameters relating to the algorithm should be here
    Setting relating to IO should be handled outside this context

    Parameters
    ----------
    ignore_tool_compounds: bool, default = False
        if true, will ignore any tool compounds during decoding
    demultiplexer_algorithm: Literal["cutadapt", "regex", "full"], default = "regex"
        The demultiplexing algorithm to use.
        - "cutadapt": use a cutadapt to locate sections
        - "regex": use a regular expression based demultiplexer
        - "full": use a full alignment based demultiplexer
    demultiplexer_mode: Literal["library", "single", "flanking"], default = "flanking"
        The demultiplexing section strategy to use.
        - "library": demultiplex by matching just the library tag
        - "single": demultiplex by matching a single static barcode section
        - "flanking": demultiplex by matching barcode sections that flank the library tag
        (flanking means one before and one after the tag)
    realign: bool, default = False
        if true, will perform a local realignment of the read to the
        libraries barcode schema *after* demultiplexing determine the library.
        This could help recover reads that have complex alignments due multiple indels
    library_error_tolerance: int, default = 1
        The number of errors you are willing to tolerate in any given barcode
        section during library demultiplexing. Will apply to each section
        independently. For example, a flanking demultiplexer will allow for
        1 error in *each* of the flanking sections.
    library_error_correction_mode_str: str, default = "levenshtein_dist:2,asymmetrical"
        The error correction mode string to use for library barcode
        calling during demultiplexing.
    min_library_overlap: int , default = 8
        if using a cutadapt style demultiplexer, this is the minimum number of bases
        that must align to the expected barcode section for a match to be called.
        See the cutadapt documentation for more details on this parameter.
    wiggle: bool, default = False
        if true, will extend aligned sections by 1 bp on each side
        and then and scan all possible chunks of the expected barcode length,
        1 smaller and 1 larger than expected length (in that order). If not
        using a local realignment post demultiplexing this can help recover
        reads lost to indels in the barcode region.
    revcomp: bool, default = False
        If true, search the reverse compliment as well.
        In most cases it is faster to use an external tools
        to align and reverse compliment reads before decoding
    max_read_length: int or None, default = None
        maximum length of a read to be considered for decoding
        if above the max, decoding will fail
        if `None` will default to 5x the min_read_length
    min_read_length: int or None, default = None
        minimum length of a read to be considered for decoding
        if below the min, decoding will fail
        if `None` will default to the smallest min match length of
        any library in the collection considered for decoding
        with 10bp of buffer
    default_error_correction_mode_str: str, default = "levenshtein_dist:1,asymmetrical"
        The default error correction mode string to use for decoding.
        If a barcode section lacks a specified error correction mode,
        this mode will be used.
        See the documentation for more details on the format of this string.
    """

    ignore_tool_compounds: bool = False,
    demultiplexer_algorithm: Literal["cutadapt", "regex", "full"] = "regex",
    demultiplexer_mode: Literal["library", "single", "flanking"] = "single",
    realign: bool = False,
    library_error_tolerance: int = 1,
    library_error_correction_mode_str: str = "levenshtein_dist:2,asymmetrical",
    min_library_overlap: int = 8,
    revcomp: bool = False,
    wiggle: bool = False,
    max_read_length: Optional[int] = None,
    min_read_length: Optional[int] = None,
    default_error_correction_mode_str: str = "levenshtein_dist:1,asymmetrical",

    def to_file(self, path: str):
        """
        Save settings to a YAML file

        Parameters
        ----------
        path : str
            path to save settings to
        """
        import yaml
        yaml.dump(asdict(self), open(path, "w"))

    @classmethod
    def from_file(cls, path: str) -> "DecodingSettings":
        """
        Load settings from a YAML file

        Will first check if there is a "decode_settings" key
        and load settings from that sub dict.
        Otherwise, will load from the YAML file keys

        Parameters
        ----------
        path : str
            Path to YAML file

        Returns
        -------
        DecodingSettings

        Raises
        ------
        RuntimeError
            if valid decode settings cannot be loaded from the passed YAML file
        """
        import yaml
        _data = yaml.safe_load(open(path, "r"))
        if "decode_settings" not in _data:
            try:
                return cls(**yaml.safe_load(open(path, "r")))
            except Exception as e:
                raise RuntimeError(f"Failed to load decode settings from {path}") from e
        else:
            try:
                return cls(**_data["decode_settings"])
            except Exception as e:
                raise RuntimeError(f"Failed to load decode settings from {path}") from e


class SelectionDecoder:
    """
    Decode reads that can come from a collection of libraries

    This means library demultiplexing is first performed

    For most users, this in the main entry point for doing decoding.

    Notes
    -----
    The runner will write logs tracking progress
    and other warnings.
    Logger will be named after the runner PID

    Parameters
    ----------
    selection: DELSelection
        the selection with the DEL libraries to decode
    decode_settings: DecodingSettings | None, default = None
        the settings to use for decoding
        if `None`, will use the default settings
    """

    def __init__(
            self,
            selection: DELSelection,
            decode_settings: DecodingSettings | None = None,
    ):
        self.selection: DELSelection = selection
        self.decode_settings: DecodingSettings = decode_settings if decode_settings is not None else DecodingSettings()
        self.decode_stats: DecodeStatistics = DecodeStatistics()

        # parse the demultiplexer settings
        demultiplex_algorithm = self.decode_settings.demultiplexer_algorithm
        demultiplex_mode = self.decode_settings.demultiplexer_mode

        demultiplexer: LibraryDemultiplexer = (get_library_demultiplexer_type(
            demultiplex_mode=demultiplex_mode,
            demultiplex_algorithm=demultiplex_algorithm)
            (
                libraries=self.selection.library_collection,
                tool_compounds=list(self.selection.tool_compounds),
                **asdict(self.decode_settings)
            )
        )

        # initialize all the decoding object required
        self.decoder = DELCollectionDecoder(
            library_demultiplexer=demultiplexer,
            decode_statistics=self.decode_stats,
            wiggle=self.decode_settings.wiggle,
            max_read_length=self.decode_settings.max_read_length,
            min_read_length=self.decode_settings.min_read_length,
            default_error_correction_mode_str=self.decode_settings.default_error_correction_mode_str,
        )

    def decode_read(self, read: SequenceRecord) -> DecodedCompound | FailedDecodeAttempt:
        """
        Decode a single read

        Parameters
        ----------
        read: SequenceReader
            the read to decode

        Returns
        -------
        DecodedDELCompound or FailedDecodeAttempt
            The decoded compound if successful,
            else a `FailedDecodeAttempt` explaining why it failed.
        """
        self.decode_stats.num_seqs_read += 1
        decoded_compound = self.decoder.decode_read(read)
        if isinstance(decoded_compound, DecodedCompound):
            self.decode_stats.num_seqs_decoded_per_lib[decoded_compound.get_library_id()] += 1
        return decoded_compound

    def decode_file(self, fastq_file: os.PathLike, use_tqdm: bool = True) -> Iterator[DecodedCompound | FailedDecodeAttempt]:
        """
        Decode all reads from a fastq file

        Parameters
        ----------
        fastq_file: os.PathLike
            path to the fastq file to decode reads from
        use_tqdm: bool, default = True
            if true, will wrap the iterator in a tqdm progress bar

        Yields
        ------
        DecodedDELCompound or FailedDecodeAttempt
            The decoded compound if successful,
            else a `FailedDecodeAttempt` explaining why it failed.
        """
        from deli.dna.io import SingleFileSequenceReader

        reader = SingleFileSequenceReader(fastq_file)
        for sequence in tqdm(reader.iter_seqs(), disable=not use_tqdm, desc=f"Decoding file: {fastq_file}"):
            yield self.decode_read(sequence)
