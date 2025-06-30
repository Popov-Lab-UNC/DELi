"""code for calling DNA barcodes"""

import abc
import json
import os
import sys
import warnings
from collections import defaultdict
from typing import Literal, no_type_check

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from dnaio import SequenceRecord
from numba import njit

from deli.dels import BuildingBlock, DELCompound, DELibrary, DELibraryCollection
from deli.dna import Aligner, HybridSemiGlobalAligner, SemiGlobalAligner

from .bb_calling import BuildingBlockSetTagCaller, FailedBuildingBlockCall, ValidBuildingBlockCall
from .lib_calling import LibraryCaller, SingleReadLibraryCaller, ValidLibraryCall
from .umi import UMI


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
    num_seqs_decoded_per_lib: defaultdict[str, int]
        the number of sequences decoded per library
    num_seqs_degen_per_lib: defaultdict[str, int]
        the number of sequences degened per library
    num_failed_too_short: int
        the number of decoded failed because observed_seq read was too short
    num_failed_too_long: int
        the number of decoded failed because observed_seq read was too long
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
        self.seq_lengths: defaultdict[int, int] = defaultdict(int)
        self.num_seqs_decoded_per_lib: defaultdict[str, int] = defaultdict(int)
        self.num_seqs_degen_per_lib: defaultdict[str, int] = defaultdict(int)

        # track the unique failures
        self.num_failed_too_short: int = 0
        self.num_failed_too_long: int = 0
        self.num_failed_library_call: int = 0
        self.num_failed_library_match_too_short: int = 0
        self.num_failed_building_block_call: int = 0
        self.num_failed_alignment: int = 0
        self.num_failed_umi_match_too_short: int = 0

    def __add__(self, other) -> "DecodeStatistics":
        """Add two DecodeStatistics objects together"""
        if not isinstance(other, DecodeStatistics):
            raise TypeError(f"unsupported operand type(s) for +: {type(self)} and {type(other)}")
        result = DecodeStatistics()
        result.num_seqs_read = self.num_seqs_read + other.num_seqs_read
        result.num_failed_too_short = self.num_failed_too_short + other.num_failed_too_short
        result.num_failed_too_long = self.num_failed_too_long + other.num_failed_too_long
        result.num_failed_library_call = (
            self.num_failed_library_call + other.num_failed_library_call
        )
        result.num_failed_library_match_too_short = (
            self.num_failed_library_match_too_short + other.num_failed_library_match_too_short
        )
        result.num_failed_building_block_call = (
            self.num_failed_building_block_call + other.num_failed_building_block_call
        )
        result.num_failed_alignment = self.num_failed_alignment + other.num_failed_alignment
        result.num_failed_umi_match_too_short = (
            self.num_failed_umi_match_too_short + other.num_failed_umi_match_too_short
        )

        # Merge defaultdicts
        result.seq_lengths = defaultdict(
            int,
            {
                k: self.seq_lengths[k] + other.seq_lengths[k]
                for k in set(self.seq_lengths) | set(other.seq_lengths)
            },
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

        result.seq_lengths = defaultdict(
            int, {int(key): int(val) for key, val in data.get("seq_lengths", {}).items()}
        )
        result.num_seqs_decoded_per_lib = defaultdict(
            int,
            {str(key): int(val) for key, val in data.get("num_seqs_decoded_per_lib", {}).items()},
        )
        result.num_seqs_degen_per_lib = defaultdict(
            int,
            {str(key): int(val) for key, val in data.get("num_seqs_degen_per_lib", {}).items()},
        )
        result.num_failed_too_short = data.get("num_failed_too_short", 0)
        result.num_failed_too_long = data.get("num_failed_too_long", 0)
        result.num_failed_library_call = data.get("num_failed_library_call", 0)
        result.num_failed_library_match_too_short = data.get(
            "num_failed_library_match_too_short", 0
        )
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
    """

    def __init__(
        self,
        library: DELibrary,
        building_blocks: list[BuildingBlock],
        umi: UMI | None = None,
    ):
        """
        Initialize the DecodedDELCompound object

        Parameters
        ----------
        library : LibraryCall
            the called library
        building_blocks : list[BuildingBlock]
            the building block calls
            if no calls were attempted (because library failed) pass a `None`
        umi: UMI or `None`
            the UMI for the read
            if not using umi or no umi in the barcode, use a `None`
        """
        super().__init__(library=library, building_blocks=building_blocks)
        self.umi = umi


class FailedDecode:
    """
    Base class for all failed decodes

    Failed decodes are really just placeholders
    to help track why a decode failed.
    """

    pass


class ReadTooShort(FailedDecode):
    """returned if read to decode was too short (based on settings)"""

    pass


class ReadTooLong(FailedDecode):
    """returned if read to decode was too long (based on settings)"""

    pass


class LibraryLookupFailed(FailedDecode):
    """returned if library lookup was not successful"""

    pass


class LibraryMatchTooShort(FailedDecode):
    """returned if library lookup resulted in a match that was too short to call"""

    pass


class BuildingBlockLookupFailed(FailedDecode):
    """
    Returned if building block lookup was not successful

    This will be true if any of the building blocks failed
    Even if you got 99/100, it is still a failure
    """

    pass


class AlignmentFailed(FailedDecode):
    """Returned if the alignment is not successful during calling"""

    pass


class UMIMatchTooShort(FailedDecode):
    """Returned if UMI match is too short to call post decode"""

    pass


class DELCollectionDecoder:
    """
    A decoder used to decode barcodes from selection using a collection of libraries

    For most use cases, this is the main entry point for decoding.
    Use one of the LibraryDecoders if you only have one (1) DEL in your selection
    As this will be far faster
    """

    def __init__(
        self,
        library_collection: DELibraryCollection,
        library_error_tolerance: float = 0.1,
        min_library_overlap: int | None = None,
        revcomp: bool = False,
        read_type: Literal["single", "paired"] = "single",
        alignment_algorithm: Literal["semi", "hybrid"] = "semi",
        bb_calling_approach: Literal["alignment", "bio"] = "bio",
        max_read_length: int | None = None,
        min_read_length: int | None = None,
        disable_error_correction: bool = False,
        decode_statistics: DecodeStatistics | None = None,
    ):
        """
        Initialize the DELCollectionDecoder

        Parameters
        ----------
        library_collection: DELibraryCollection
            the library collection used for the selection to decode
        library_error_tolerance: float, default = 0.2
            the percent error to be tolerated in the library section
            this will be converted to number of errors based on tag size
            and down to the nearest whole number.
            for example, a library with 14 nucleotides would tolerate
            1, 2, and 4 errors for an error tolerance of 0.1, 0.2 and 0.3 respectively
        min_library_overlap: int or None, default = 7
            the minimum number of nucleotides required to match
            the library tag
            This is because the demultiplexing will accept truncated matches
            at the front/back of the tag. For example, a tag of AGCTGGTTC
            could match a read of GTTC if the min overlap was <=4
            If `None`, will default to the exact length of the tag, meaning
            the whole tag is expected.
            The recommended value is greater than 8, as the odds of a match this strong
            to be accidental are low
        alignment_algorithm: Literal["semi", "hybrid"], default = "semi"
            the algorithm to use for alignment
            only used if bb_calling_approach is "alignment"
        read_type: Literal["single", "paired"], default = "single"
            the type of read
            paired are for paired reads
            all other read types are single
        revcomp: bool, default = False
            If true, search the reverse compliment as well
        max_read_length: int or None, default = None
            maximum length of a read to be considered for decoding
            if above the max, decoding will fail
            if `None` will default to 5x the min_read_length
        min_read_length: int or None, default = None
            minimum length of a read to be considered for decoding
            if below the min, decoding will fail
            if `None` will default to smallest min match length of
            any library in the collection considered for decoding
        bb_calling_approach: Literal["alignment"], default = "alignment"
            the algorithm to use for bb_calling
            right now only "alignment" mode is supported
        disable_error_correction: bool, default = False
            disable error correction for any barcode section
            capable of it
        decode_statistics: DecodeStatistics or None, default = None
            the statistic tracker for the decoding run
            if None, will initialize a new, empty statistic object
        """
        self.decode_statistics: DecodeStatistics
        if decode_statistics is None:
            self.decode_statistics = DecodeStatistics()
        else:
            self.decode_statistics = decode_statistics

        self.library_caller: LibraryCaller
        if read_type == "single":
            self.library_caller = SingleReadLibraryCaller(
                library_collection,
                error_rate=library_error_tolerance,
                min_overlap=min_library_overlap,
                revcomp=revcomp,
            )
        elif read_type == "paired":
            warnings.warn(
                "DELi only support single reads currently; if a paired mode would be nice"
                "please raise an issue to request it",
                stacklevel=1,
            )
            self.library_caller = SingleReadLibraryCaller(
                library_collection,
                error_rate=library_error_tolerance,
                min_overlap=min_library_overlap,
                revcomp=revcomp,
            )
        else:
            raise ValueError(f"Unknown read_type '{read_type}'")

        # determine library caller class
        self.library_decoders: dict[str, LibraryDecoder]
        if bb_calling_approach == "alignment":
            self.library_decoders = {
                library.library_id: DynamicAlignmentLibraryDecoder(
                    library,
                    decode_statistics=decode_statistics,
                    disable_error_correction=disable_error_correction,
                    alignment_algorithm=alignment_algorithm,
                )
                for library in library_collection.libraries
            }
        elif bb_calling_approach == "bio":
            self.library_decoders = {
                library.library_id: BioAlignmentLibraryDecoder(
                    library,
                    decode_statistics=self.decode_statistics,
                    disable_error_correction=disable_error_correction,
                )
                for library in library_collection.libraries
            }
        else:
            raise ValueError(f"unrecognized bb_calling_approach {bb_calling_approach}")

        # set the min/max lengths
        self._min_read_length: int
        if isinstance(min_read_length, int):
            self._min_read_length = min_read_length
        else:
            self._min_read_length = min(
                [lib.barcode_schema.min_length for lib in library_collection.libraries]
            )

        self._max_read_length: int
        if isinstance(max_read_length, int):
            self._max_read_length = max_read_length
        else:
            self._max_read_length = 5 * self._min_read_length

    def decode_read(self, sequence: SequenceRecord) -> DecodedDELCompound | FailedDecode:
        """
        Given a observed_seq read, decode its barcode

        """
        # check lengths
        if len(sequence) > self._max_read_length:
            self.decode_statistics.num_failed_too_long += 1
            return ReadTooLong()
        if len(sequence) < self._min_read_length:
            self.decode_statistics.num_failed_too_short += 1
            return ReadTooShort()

        self.decode_statistics.seq_lengths[len(sequence)] += 1
        library_call = self.library_caller.call_library(sequence)
        if isinstance(library_call, ValidLibraryCall):
            return self.library_decoders[library_call.library.library_id].call_barcode(
                library_call
            )
        else:
            # if library calling fails cannot continue decoding
            self.decode_statistics.num_failed_library_call += 1
            return LibraryLookupFailed()


class LibraryDecoder(abc.ABC):
    """
    Abstract class for all library decoders

    A library decoder assumes the library the read is from is already known.
    These decoders are only concerned with calling or locating all other sections.
    Currently, that is the building block sections and UMI.
    If you would like support for calling another section too, raise an issue
    to add that feature
    """

    def __init__(
        self,
        library: DELibrary,
        decode_statistics: DecodeStatistics | None = None,
        disable_error_correction: bool = False,
    ):
        """
        Initialize the LibraryDecoder

        Parameters
        ----------
        library: DELibrary
            the library to decode from
        decode_statistics: DecodeStatistics or `None`, default = `None`
            the statistic tracker for the decoding run
            if `None` will initialize a new, empty statistic object
        disable_error_correction: bool, default False
            disable error correction for any barcode section
            capable of it
        """
        self.disable_error_correction = disable_error_correction
        self.library = library

        self.decode_statistics: DecodeStatistics
        if decode_statistics is None:
            self.decode_statistics = DecodeStatistics()
        else:
            self.decode_statistics = decode_statistics

    @abc.abstractmethod
    def call_barcode(self, library_call: ValidLibraryCall) -> DecodedDELCompound | FailedDecode:
        """Given a library call, decode the barcode"""
        raise NotImplementedError()


class DynamicAlignmentLibraryDecoder(LibraryDecoder):
    """
    Uses a dynamic programing alignment algorithm to call the barcode

    Will align to the library tag (and other known static regions) and use
    that alignment to locate the required sections to call
    """

    def __init__(
        self,
        library: DELibrary,
        decode_statistics: DecodeStatistics | None = None,
        alignment_algorithm: Literal["semi", "hybrid"] = "semi",
        disable_error_correction: bool = False,
    ):
        """
        Initialize the DynamicAlignmentLibraryDecoder

        Parameters
        ----------
        library: DELibrary
            the library to decode from
        decode_statistics: DecodeStatistics or `None`, default = `None`
            the statistic tracker for the decoding run
            if `None` will initialize a new, empty statistic object
        alignment_algorithm: Literal["semi", "hybrid"], default "semi"
            the type alignment algorithm to use
            semi is a semi global and hybrid is a hybrid semi global alignment
            see aligning docs for more details
        disable_error_correction: bool, default False
            disable error correction for any barcode section
            capable of it
        """
        super().__init__(
            library=library,
            decode_statistics=decode_statistics,
            disable_error_correction=disable_error_correction,
        )

        # assign the aligner
        self.aligner: Aligner
        if alignment_algorithm == "semi":
            self.aligner = SemiGlobalAligner()
        elif alignment_algorithm == "hybrid":
            self.aligner = HybridSemiGlobalAligner()
        else:
            raise ValueError(
                f"unknown building block alignment algorithm: '{alignment_algorithm}';"
                f"currently supported values are ['semi', 'hybrid']"
            )

        self._bb_callers: dict[str, BuildingBlockSetTagCaller] = {
            bb_section.section_name: BuildingBlockSetTagCaller(
                building_block_tag_section=bb_section,
                building_block_set=bb_set,
                disable_error_correction=self.disable_error_correction,
            )
            for bb_section, bb_set in self.library.iter_bb_barcode_sections_and_sets()
        }

        self._barcode_reference = self.library.barcode_schema.get_full_barcode()

        _barcode_section_spans = self.library.barcode_schema.get_section_spans(
            exclude_overhangs=True
        )
        self._bb_sections_spans_to_search_for: dict[str, slice] = {
            section.section_name: _barcode_section_spans[section.section_name]
            for section in self.library.barcode_schema.building_block_sections
        }
        self._umi_span: slice | None = _barcode_section_spans.get("umi", None)

        # check the alignment of callers and sections:
        _missing_sections = set(self._bb_callers.keys()) - set(
            self._bb_sections_spans_to_search_for.keys()
        )
        if len(_missing_sections) != 0:
            raise ValueError(
                f"building block caller for section(s) '{_missing_sections}'"
                f"are missing from the library barcode"
            )

    def call_barcode(self, library_call: ValidLibraryCall) -> DecodedDELCompound | FailedDecode:
        """
        Given a library call, decode the read

        Parameters
        ----------
        library_call: ValidLibraryCall
            the library call to decode

        Returns
        -------
        DecodedDELCompound
        """
        _alignment = self.aligner.align(self._barcode_reference, library_call.sequence.sequence)

        _barcode_section_alignment: dict[str, str] = dict()
        for section_name, section_span in self._bb_sections_spans_to_search_for.items():
            _aligned_span = _alignment.map_seq1_span_to_seq2_span(
                section_span.start, section_span.stop
            )
            _barcode_section_alignment[section_name] = library_call.sequence.sequence[
                slice(*_aligned_span)
            ]

        # call the building blocks
        bb_calls: dict[str, ValidBuildingBlockCall] = dict()
        for bb_section_name, bb_caller in self._bb_callers.items():
            bb_call = bb_caller.call_building_block(_barcode_section_alignment[bb_section_name])
            if isinstance(bb_call, ValidBuildingBlockCall):
                bb_calls[bb_section_name] = bb_call
            else:
                self.decode_statistics.num_failed_building_block_call += 1
                return BuildingBlockLookupFailed()

        # extract the UMI section if there is one
        _umi: UMI | None = None
        if self._umi_span is not None:
            _umi = UMI(
                library_call.sequence.sequence[
                    slice(
                        *_alignment.map_seq1_span_to_seq2_span(
                            self._umi_span.start, self._umi_span.stop
                        )
                    )
                ]
            )

        return DecodedDELCompound(
            library=library_call.library,
            building_blocks=[val.building_block for val in bb_calls.values()],
            umi=_umi,
        )


class BioAlignmentLibraryDecoder(LibraryDecoder):
    """
    Uses a dynamic programing alignment algorithm to call the barcode

    Will align to the library tag (and other known static regions) and use
    that alignment to locate the required sections to call
    """

    def __init__(
        self,
        library: DELibrary,
        decode_statistics: DecodeStatistics | None = None,
        disable_error_correction: bool = False,
    ):
        """
        Initialize the DynamicAlignmentLibraryDecoder

        Parameters
        ----------
        library: DELibrary
            the library to decode from
        decode_statistics: DecodeStatistics or `None`, default = `None`
            the statistic tracker for the decoding run
            if `None` will initialize a new, empty statistic object
        disable_error_correction: bool, default False
            disable error correction for any barcode section
            capable of it
        """
        super().__init__(
            library=library,
            decode_statistics=decode_statistics,
            disable_error_correction=disable_error_correction,
        )

        # define a custom substitution matrix
        # when aligning to N we don't want to penalize for mismatches
        # since N represents an unknown base we need to decode
        alphabet = "ACGTN"
        matrix = substitution_matrices.Array(alphabet, dims=2)
        for base1 in alphabet:
            for base2 in alphabet:
                if base1 == base2:
                    matrix[base1, base2] = 1  # matches score 1
                elif "N" in (base1, base2):
                    matrix[base1, base2] = 0  # mismatches with 'N' score 0
                else:
                    matrix[base1, base2] = -1  # other mismatches score -1

        # assign the aligner
        self.aligner = PairwiseAligner(
            match_score=1,
            mismatch_score=-1,
            open_gap_score=-2,
            extend_gap_score=-1,
            mode="global",
            end_open_gap_score=0,
            end_extend_gap_score=0,
            substitution_matrix=matrix,
        )

        self._bb_callers: dict[str, BuildingBlockSetTagCaller] = {
            bb_section.section_name: BuildingBlockSetTagCaller(
                building_block_tag_section=bb_section,
                building_block_set=bb_set,
                disable_error_correction=self.disable_error_correction,
            )
            for bb_section, bb_set in self.library.iter_bb_barcode_sections_and_sets()
        }
        self._num_bb_to_call = len(self._bb_callers)

        self._barcode_reference = self.library.barcode_schema.get_full_barcode()

        _barcode_section_spans = self.library.barcode_schema.get_section_spans(
            exclude_overhangs=True
        )
        self._bb_sections_spans_to_search_for: dict[str, slice] = {
            section.section_name: _barcode_section_spans[section.section_name]
            for section in self.library.barcode_schema.building_block_sections
        }
        self._umi_span: slice | None = _barcode_section_spans.get("umi", None)

        # check the alignment of callers and sections:
        _missing_sections = set(self._bb_callers.keys()) - set(
            self._bb_sections_spans_to_search_for.keys()
        )
        if len(_missing_sections) != 0:
            raise ValueError(
                f"building block caller for section(s) '{_missing_sections}'"
                f"are missing from the library barcode"
            )

        # pre compile _inverse_indices with dummy data
        # this could make the initialization of the decoder take longer
        # if this function has not been compiled yet
        _tmp_alignment = self.aligner.align("AGCT", "AGCT")[0]
        _inverse_indices(_tmp_alignment.sequences, _tmp_alignment.coordinates)

    def call_barcode(self, library_call: ValidLibraryCall) -> DecodedDELCompound | FailedDecode:
        """
        Given a library call, decode the read

        Parameters
        ----------
        library_call: ValidLibraryCall
            the library call to decode

        Returns
        -------
        DecodedDELCompound
        """
        _alignments = self.aligner.align(library_call.sequence.sequence, self._barcode_reference)

        # if biopython finds too many top scoring alignments, break out
        #  it means that no good alignment exists, so calling is too risky
        if len(_alignments) > 20:
            self.decode_statistics.num_failed_alignment += 1
            return AlignmentFailed()

        # loop through all scoring alignments, skipping if at any point a bb lookup
        # fails, and exiting and returning the first perfect match
        for _alignment in _alignments:
            # _inverse_indices = _alignment.inverse_indices
            inverse_indices = _inverse_indices(_alignment.sequences, _alignment.coordinates)
            bb_calls: dict[str, ValidBuildingBlockCall] = dict()
            for section_name, section_span in self._bb_sections_spans_to_search_for.items():
                _aligned_bb_codon = library_call.sequence.sequence[
                    inverse_indices[1][section_span.start] : inverse_indices[1][
                        section_span.stop - 1
                    ]
                    + 1
                ]
                bb_call = self._bb_callers[section_name].call_building_block(_aligned_bb_codon)
                if isinstance(bb_call, FailedBuildingBlockCall):
                    break
                else:
                    bb_calls[section_name] = bb_call

            if len(bb_calls) != self._num_bb_to_call:
                continue

            # extract the UMI section if there is one
            _umi: UMI | None = None
            if self._umi_span is not None:
                _umi = UMI(
                    library_call.sequence.sequence[
                        inverse_indices[1][self._umi_span.start] : inverse_indices[1][
                            self._umi_span.stop - 1
                        ]
                        + 1
                    ]
                )

            # check umi length
            if (_umi is not None) and (
                len(_umi) < self.library.barcode_schema.get_section_length("umi")
            ):
                warnings.warn(
                    f"UMI length {len(_umi)} is shorter than expected "
                    f"{self.library.barcode_schema.get_section_length('umi')}; "
                    f"This is likely the result of a major gap in an otherwise "
                    f"perfect (decodable) the alignment. "
                    f"This may indicate a problem with missing or extra barcode "
                    f"sections in your library definition.",
                    stacklevel=1,
                )
                return UMIMatchTooShort()

            return DecodedDELCompound(
                library=library_call.library,
                building_blocks=[val.building_block for val in bb_calls.values()],
                umi=_umi,
            )
        self.decode_statistics.num_failed_building_block_call += 1
        return BuildingBlockLookupFailed()


# TODO for some reason, this is needed to prevent crash when running as a CMD,
#  but not in a python interpreter... not sure why?
sys.setrecursionlimit(3000)


@no_type_check
@njit
def _inverse_indices(sequences, coordinates):
    """
    Numba accelerated inverse_indices calculation

    This is lifted from biopython/Bio/Align/__init__.py:inverse_indices:2961
    The inverse_indices attribute is actually a property for alignment objects
    in Biopython, meaning it is just run every time you try and access the
    inverse_indices attribute.
    It is the most expensive part of the alignment, although
    The code can be numba compiled to speed up the frequent calls
    to the function we have to make.
    This function does just that.
    A few changes to the exact functions were made, but the algorithm is the
    same, just compatible with njit now.

    Using njit makes this function about 5 times faster, and since it accounts for
    ~25% of the runtime during decoding, it offers tangible speed ups

    WARNING:
    Users should never touch this function, it is only meant
    to be used by the BioAlignmentLibraryDecoder to make things
    go a bit faster

    Parameters
    ----------
    sequences
        the sequences from the alignment object
    coordinates
        the coordinates from the alignment object

    Returns
    -------
    reverse indices
    """
    a = [np.full(len(sequence), -1) for sequence in sequences]
    n, m = coordinates.shape
    steps = np.diff(coordinates, 1)
    steps_bool = steps != 0
    aligned = np.sum(steps_bool, 0) > 1
    # True for steps in which at least two sequences align, False if a gap
    steps = steps[:, aligned]
    rcs = np.zeros(n)
    for i, row in enumerate(steps):
        if (row >= 0).all():
            rcs[i] = False
        elif (row <= 0).all():
            rcs[i] = True
        else:
            raise ValueError(f"Inconsistent steps in row {i}")
    i = 0
    j = 0
    for k in range(m - 1):
        starts = coordinates[:, k]
        ends = coordinates[:, k + 1]
        for row, start, end, rc in zip(a, starts, ends, rcs):
            if rc == False and start < end:  # noqa: E712
                j = i + end - start
                row[start:end] = np.arange(i, j)
            elif rc == True and start > end:  # noqa: E712
                j = i + start - end
                if end > 0:
                    row[start - 1 : end - 1 : -1] = np.arange(i, j)
                elif start > 0:
                    row[start - 1 :: -1] = np.arange(i, j)
        i = j
    return a


#####
# Calling by regex alignment has be removed, but these functions
# will stick around incase we want to roll that back
#
#
# def _between(
#         val: int,
#         start: int,
#         stop: int,
#         right_inclusive: bool = True,
#         left_inclusive: bool = False):
#     """
#     Helper function for determining if some value is between two values
#
#     Parameters
#     ----------
#     val: int
#         value to check if in between
#     start: int
#         start of the between range
#     stop: int
#         stop of the between range
#     right_inclusive: bool = True
#         include the `start` value in the range
#     left_inclusive: bool = True
#         include the `stop` value in the range
#
#     Returns
#     -------
#     bool
#         if `val` is between `start` and `stop`
#     """
#     return (start - right_inclusive) < val < (stop + left_inclusive)
#
#
# def _match_call_mode(
#         match: FullBarcodeMatch,
#         barcode_schema: BarcodeSchema
# ) -> BarcodeAlignment:
#     """
#     Call mode function used to generate section alignment in match mode
#
#     Notes
#     -----
#     Needs FullBarcodeMatches (generated with the `full_match=True` flag
#     during the matching)
#
#     Parameters
#     ----------
#     match: FullBarcodeMatch
#         the match to generate the alignment for
#     barcode_schema:
#         the schema defining the barcode to align to
#
#     Returns
#     -------
#     BarcodeAlignment
#     """
#     _tmp_alignment = barcode_schema.barcode_spans.copy()
#
#     _errors: List[Tuple[int, str]] =
#     list(zip(*match.get_sorted_match_sequence_error_idx(exclude_substitute=True)))
#     if len(_errors) == 0:
#         return BarcodeAlignment(_tmp_alignment)
#
#     (_err_idx, _err_type) = _errors.pop(0)
#     _buff = 0
#     _spans = barcode_schema.barcode_spans.copy()
#     for _key, (_start, _stop) in _tmp_alignment.items():
#         _start = _start + _buff
#         _stop = _stop + _buff
#         _end = _stop
#         while _err_idx < _end:
#             if _between(_err_idx, _start, _stop):
#                 if _err_type == "insert":
#                     _stop += 1
#                 elif _err_type == "delete":
#                     _stop -= 1
#                 else:
#                     raise ValueError(f"unrecognized error type {_err_type}")
#
#                 try:
#                     (_err_idx, _err_type) = _errors.pop(0)
#                 except IndexError:
#                     (_err_idx, _err_type) = (1e9, "no-more-errors")
#         _spans[_key] = (_start, _stop)
#         _buff += (_stop - _end)
#     return BarcodeAlignment(_spans)
####
