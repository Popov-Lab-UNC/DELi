"""code for calling DNA barcodes"""

import json
import os
from collections import defaultdict
from typing import Optional

from dnaio import SequenceRecord

from deli.dels.barcode import BarcodeSection
from deli.dels.building_block import TaggedBuildingBlock
from deli.dels.compound import DELCompound
from deli.dels.library import DELibrary

from ._base import FailedDecodeAttempt
from .barcode_calling import AmbiguousBarcodeCall, BarcodeCaller, FailedBarcodeLookup, ValidCall, get_barcode_caller
from .library_demultiplex import AlignedSeq, LibraryDemultiplexer
from .umi import UMI


MAX_RETRIES = 50  # number of alignments to attempt before giving up


class DecodeStatistics:
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
    seq_lengths: dict[int, int]
        the distribution of read lengths observed
    num_seqs_decoded_per_lib: dict[str, int]
        the number of sequences decoded per library
    num_seqs_degen_per_lib: dict[str, int]
        the number of sequences degened per library
    num_failed_too_short: int
        the number of decoded failed because barcode read was too short
    num_failed_too_long: int
        the number of decoded failed because barcode read was too long
    num_failed_library_call: int
        the number of decoded failed because the library was not called
    num_failed_building_block_call: int
        the number of decoded failed because a building block was not called
    num_failed_alignment: int
        the number of decoded failed because alignment failed
    """

    def __init__(self):
        """Initialize a DecodeStatistics object"""
        self.num_seqs_read: int = 0
        self.seq_lengths: defaultdict[int, int] = defaultdict(int)
        self.num_seqs_decoded_per_lib: defaultdict[str, int] = defaultdict(int)
        self.num_seqs_degen_per_lib: defaultdict[str, int] = defaultdict(int)

        # track the unique failures
        self.num_failed_too_short: int = 0
        self.num_failed_too_long: int = 0
        self.num_failed_library_call: int = 0
        self.num_failed_building_block_call: int = 0
        self.num_failed_alignment: int = 0

    def __add__(self, other) -> "DecodeStatistics":
        """Add two DecodeStatistics objects together"""
        if not isinstance(other, DecodeStatistics):
            raise TypeError(f"unsupported operand type(s) for +: {type(self)} and {type(other)}")
        result = DecodeStatistics()
        result.num_seqs_read = self.num_seqs_read + other.num_seqs_read
        result.num_failed_too_short = self.num_failed_too_short + other.num_failed_too_short
        result.num_failed_too_long = self.num_failed_too_long + other.num_failed_too_long
        result.num_failed_library_call = self.num_failed_library_call + other.num_failed_library_call
        result.num_failed_building_block_call = (
            self.num_failed_building_block_call + other.num_failed_building_block_call
        )
        result.num_failed_alignment = self.num_failed_alignment + other.num_failed_alignment

        # Merge defaultdicts
        result.seq_lengths = defaultdict(
            int,
            {k: self.seq_lengths[k] + other.seq_lengths[k] for k in set(self.seq_lengths) | set(other.seq_lengths)},
        )
        result.num_seqs_decoded_per_lib = defaultdict(
            int,
            {
                k: self.num_seqs_decoded_per_lib[k] + other.num_seqs_decoded_per_lib[k]
                for k in set(self.num_seqs_decoded_per_lib) | set(other.num_seqs_decoded_per_lib)
            },
        )
        result.num_seqs_degen_per_lib = defaultdict(
            int,
            {
                k: self.num_seqs_degen_per_lib[k] + other.num_seqs_degen_per_lib[k]
                for k in set(self.num_seqs_degen_per_lib) | set(other.num_seqs_degen_per_lib)
            },
        )
        return result

    def __str__(self) -> str:
        """Convert the statistic object to a string (new line seperated)"""
        return "\n".join([f"{key}={val}\n" for key, val in self.__dict__.items()])

    def __repr__(self) -> str:
        """Represent the statistic object as string ('; ' seperated)"""
        return "; ".join([f"{key}={val}\n" for key, val in self.__dict__.items()])

    @property
    def num_seqs_decoded(self) -> int:
        """Number of sequences decoded successfully in total"""
        return sum(self.num_seqs_decoded_per_lib.values())

    @property
    def num_seqs_degen(self) -> int:
        """Number of degen sequences observed in total"""
        return sum(self.num_seqs_degen_per_lib.values())

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

        result.seq_lengths = defaultdict(int, {int(key): int(val) for key, val in data.get("seq_lengths", {}).items()})
        result.num_seqs_decoded_per_lib = defaultdict(
            int,
            {str(key): int(val) for key, val in data.get("num_seqs_decoded_per_lib", {}).items()},
        )
        result.num_seqs_degen_per_lib = defaultdict(
            int,
            {str(key): int(val) for key, val in data.get("num_seqs_degen_per_lib", {}).items()},
        )
        # collect error info
        result.num_failed_too_short = data.get("num_failed_too_short", 0)
        result.num_failed_too_long = data.get("num_failed_too_long", 0)
        result.num_failed_library_call = data.get("num_failed_library_call", 0)
        result.num_failed_building_block_call = data.get("num_failed_building_block_call", 0)
        result.num_failed_alignment = data.get("num_failed_alignment", 0)
        return result


class DecodedDELCompound(DELCompound):
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


class AmbigousBuildingBlockBarcode(FailedDecodeAttempt):
    """returned if building block barcode is ambiguous"""

    def __init__(self, sequence, barcode: str, bb_section_name: str):
        super().__init__(
            sequence=sequence,
            reason=f"ambiguous barcode for building block barcode '{barcode}' in section '{bb_section_name}'",
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
    library_decoders: dict[DELibrary, LibraryDecoder]
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
        self.library_decoders: dict[DELibrary, LibraryDecoder] = {
            library: LibraryDecoder(
                library=library, wiggle=wiggle, default_error_correction_mode_str=default_error_correction_mode_str
            )
            for library in self.library_demultiplexer.libraries
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

    def decode_read(self, sequence: SequenceRecord) -> DecodedDELCompound | FailedDecodeAttempt:
        """
        Given a raw read, decode it into a DecodedDELCompound

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
        DecodedDELCompound or FailedDecodeAttempt
        """
        # check lengths
        if len(sequence) > self._max_read_length:
            self.decode_statistics.num_failed_too_long += 1
            return ReadTooLong(sequence)
        if len(sequence) < self._min_read_length:
            self.decode_statistics.num_failed_too_short += 1
            return ReadTooShort(sequence)

        self.decode_statistics.seq_lengths[len(sequence)] += 1

        _call = self.library_demultiplexer.demultiplex(sequence)
        if isinstance(_call, FailedDecodeAttempt):
            self.decode_statistics.num_failed_library_call += 1
            return _call
        else:
            library_call, alignment_iter = _call
            library_decoder = self.library_decoders[library_call.obj]

            bb_calls: FailedDecodeAttempt | list[ValidCall[TaggedBuildingBlock]] = FailedBuildingBlockCall(
                sequence, "Null", "all"
            )  # a generic default, should never be returned
            umi_call: UMI | None = None

            # TODO could speed up by just failing if too many alignments
            for attempt_count, (alignment, _) in enumerate(alignment_iter):
                if attempt_count > MAX_RETRIES:
                    return AlignmentFailed(sequence)

                # call building blocks
                bb_calls = library_decoder.call_building_blocks(alignment)
                if not isinstance(bb_calls, FailedDecodeAttempt):
                    umi_call = library_decoder.call_umi(alignment)
                    break
                else:
                    continue

            if not isinstance(bb_calls, FailedDecodeAttempt):
                return DecodedDELCompound(
                    library_call=library_call,
                    building_block_calls=bb_calls,
                    umi=umi_call,
                )
            else:
                self.decode_statistics.num_failed_building_block_call += 1
                return bb_calls


class BarcodeDecodingError(Exception):
    """raised when there is an error during barcode decoding"""

    pass


class LibraryDecoder:
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

    Attributes
    ----------
    bb_callers: dict[BarcodeSection, BarcodeCaller]
        the barcode callers for each building block section in the library
    """

    def __init__(
        self,
        library: DELibrary,
        wiggle: bool = False,
        default_error_correction_mode_str: str = "levenshtein_dist:1,asymmetrical",
    ):
        self.library = library
        self.wiggle = wiggle

        self.bb_callers: dict[BarcodeSection, BarcodeCaller] = {}
        for bb_sec, bb_set in self.library.iter_bb_barcode_sections_and_sets():
            error_correction_mode_str = getattr(bb_set, "error_correction_mode_str", default_error_correction_mode_str)
            self.bb_callers[bb_sec] = get_barcode_caller(
                tag_map={tag: bb for bb in bb_set.building_blocks for tag in bb.tags},
                error_correction_mode_str=error_correction_mode_str,
            )

        self._dist_from_final_bb_to_umi: int = 0
        self._umi_length: int = -1
        if self.library.barcode_schema.has_umi():
            self._dist_from_final_bb_to_umi = self.library.barcode_schema.get_length_between_sections(
                self.library.barcode_schema.building_block_sections[-1].section_name,
                "umi",
                include_direction=False,  # UMI should always be after BBs
            ) - (
                len(self.library.barcode_schema.building_block_sections[-1].section_overhang)
                if self.library.barcode_schema.building_block_sections[-1].section_overhang
                else 0
            )
            self._umi_length = len(self.library.barcode_schema.get_section("umi"))
        self._last_bb_idx: int = -1

    def call_building_blocks(
        self, aligned_sequence: AlignedSeq
    ) -> list[ValidCall[TaggedBuildingBlock]] | FailedDecodeAttempt:
        """
        Given a library call, decode the barcode

        Parameters
        ----------
        aligned_sequence: AlignedSeq
            the aligned sequence to decode the building blocks from

        Returns
        -------
        list[ValidCall[TaggedBuildingBlock]] or FailedDecodeAttempt
            the building block calls if successful,
            else a FailedDecodeAttempt explaining why it failed
        """
        bb_calls: list[ValidCall[TaggedBuildingBlock]] = []

        _last_bb_idx = -1
        for bb_sec, bb_caller in self.bb_callers.items():
            section_codon = aligned_sequence.get_section_barcode(bb_sec.section_name)
            _last_bb_idx = aligned_sequence.section_spans[bb_sec.section_name][1]
            bb_call = bb_caller.decode_barcode(section_codon)

            if isinstance(bb_call, FailedBarcodeLookup) and self.wiggle:
                true_section_length = len(bb_sec.section_tag)
                best_score = 100.0
                # default to failed lookup
                best_call: ValidCall[DELibrary] | FailedBarcodeLookup = FailedBarcodeLookup("placeholder")

                start, stop = aligned_sequence.section_spans[bb_sec.section_name]
                adjusted_section_codon = aligned_sequence.sequence.sequence[start - 1 : stop + 1]

                for length in range(true_section_length, true_section_length + 1, true_section_length - 1):
                    wiggle_start = 0
                    wiggle_end = wiggle_start + length
                    while wiggle_end <= true_section_length + 2:
                        wiggle_codon = adjusted_section_codon[wiggle_start:wiggle_end]
                        wiggle_call = bb_caller.decode_barcode(wiggle_codon)
                        if not isinstance(wiggle_call, FailedBarcodeLookup):
                            # TODO could have an option to first match instead of best match
                            if wiggle_call.score == 0:  # found perfect match quit now
                                bb_call = wiggle_call
                                _last_bb_idx = stop - true_section_length + wiggle_end
                                break
                            elif wiggle_call.score < best_score:
                                best_score = wiggle_call.score
                                best_call = wiggle_call
                                _last_bb_idx = stop - true_section_length + wiggle_end
                            elif wiggle_call.score == best_score:
                                # if the barcodes are different but have the same score, we cannot decide
                                if (not isinstance(best_call, FailedBarcodeLookup)) and (
                                    wiggle_call.obj != best_call.obj
                                ):
                                    best_call = AmbiguousBarcodeCall(wiggle_codon)
                        # move the wiggle window
                        wiggle_start += 1
                        wiggle_end = wiggle_start + length

                    if not isinstance(bb_call, FailedBarcodeLookup):
                        break  # found perfect match quit now

                if best_call is not None:
                    bb_call = best_call

            if isinstance(bb_call, FailedBarcodeLookup):
                if isinstance(bb_call, AmbiguousBarcodeCall):
                    return AmbigousBuildingBlockBarcode(aligned_sequence.sequence, section_codon, bb_sec.section_name)
                return FailedBuildingBlockCall(aligned_sequence.sequence, section_codon, bb_sec.section_name)
            else:
                bb_calls.append(bb_call)
        self._last_bb_idx = _last_bb_idx
        return bb_calls

    def call_umi(self, aligned_sequence: AlignedSeq) -> UMI | None:
        """
        Given a library call, decode the UMI if present

        Will return None if no UMI is defined in the library schema

        Parameters
        ----------
        aligned_sequence: AlignedSeq
            the aligned sequence to decode the UMI from
        """
        if self._umi_length == -1:
            return None
        elif self._last_bb_idx == -1:
            raise BarcodeDecodingError("Cannot call UMI before building blocks have been called")
        else:
            umi_start = self._last_bb_idx + self._dist_from_final_bb_to_umi
            umi_stop = umi_start + self._umi_length
            umi_codon = aligned_sequence.sequence.sequence[umi_start:umi_stop]
            aligned_umi_codon = aligned_sequence.get_section_barcode("umi")
            if umi_codon != aligned_umi_codon:
                return UMI([umi_codon, aligned_umi_codon])
            else:
                return UMI([umi_codon])
