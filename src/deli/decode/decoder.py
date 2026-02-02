"""code for calling DNA barcodes"""

import abc
import dataclasses
import json
import os
from collections import defaultdict
from dataclasses import asdict
from typing import Generic, Iterator, Literal, MutableSequence, Optional, TypeVar

from dnaio import SequenceRecord
from Levenshtein import distance as levenshtein_distance

from deli.dels.barcode import BarcodeSection
from deli.dels.building_block import TaggedBuildingBlock
from deli.dels.combinatorial import DELibrary
from deli.dels.compound import Compound, DELCompound
from deli.dels.library import Library
from deli.dels.tool_compound import TaggedToolCompound, ToolCompound
from deli.enumeration.enumerator import EnumerationRunError
from deli.selection import DELSelection

from .barcode_calling import AmbiguousBarcodeCall, BarcodeCaller, FailedBarcodeLookup, ValidCall, get_barcode_caller
from .base import FailedDecodeAttempt
from .library_demultiplex import (
    AlignedSeq,
    FailedStaticAlignment,
    LibraryDemultiplexer,
    ToolCompoundLibrary_,
    get_library_demultiplexer_type,
)


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


C = TypeVar("C", bound="Compound")


class DecodedCompound(abc.ABC, Generic[C]):
    """Base class for compounds that were decoded from a DNA read"""

    def __init__(self, umi: str):
        self.umi = umi

    @abc.abstractmethod
    def to_compound(self) -> C:
        """
        Convert the decoded compound to a Compound object

        Returns
        -------
        C
        """
        raise NotImplementedError()

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
        - "bb##_id":
            the id for the building block from cycle ##
        - "bb##_smiles":
            the smiles for the building block from cycle ##

        When writing to a cube, if any key is missing it will be filed with "null"
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def get_decode_res_info(self) -> dict[str, str]:
        """
        Collect the info needed for writing to a decode results file

        Dictionary keys will be:
        - "bb_ids":
            the building block IDs for the compound, comma separated and ordered by cycle
        - "bb_scores":
            the building block call score for each call, comma separated and ordered by cycle
        - "umi":
            the UMI sequence

        These keys can be left out if not applicable, but no additional keys should be added
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

    @abc.abstractmethod
    def get_score(self) -> float:
        """
        Get the decode score of the decoded compound

        Returns
        -------
        float
        """
        raise NotImplementedError()


class DecodedDELCompound(DecodedCompound[DELCompound]):
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
    library : DELibrary
        the called library
    building_block_calls : list[ValidCall[TaggedBuildingBlock]]
        the building block calls *in order of the library's building block sections*
    umi: str or `None`
        the UMI for the read
        if not using umi or no umi in the barcode, use a `None`
    """

    def __init__(
        self,
        library: DELibrary,
        building_block_calls: list[ValidCall[TaggedBuildingBlock]],
        umi: str,
    ):
        self.library = library
        self.building_block_calls = building_block_calls
        super().__init__(umi=umi)

    def to_compound(self) -> DELCompound:
        """
        Convert the decoded DEL compound to a Compound object

        Returns
        -------
        DELCompound
        """
        return DELCompound(library=self.library, building_blocks=[bb_call.obj for bb_call in self.building_block_calls])

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
        - "bb##":
            the id for the building block from cycle ##
        - "bb##_smiles":
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
        compound = self.to_compound()

        row_dict: dict[str, str] = {
            "library_id": compound.library.library_id,
            "compound_id": compound.compound_id,
        }

        if enumerate_compounds:
            try:
                row_dict["smiles"] = self.library.enumerator.enumerate_by_bbs(compound.building_blocks)
            except EnumerationRunError:
                row_dict["smiles"] = "ENUMERATION_FAILED"

        for idx, bb in enumerate(compound.building_blocks, start=1):
            row_dict[f"bb{idx}"] = bb.bb_id
            row_dict[f"bb{idx}_smiles"] = bb.smi if bb.has_smiles() else "null"
        return row_dict

    def get_decode_res_info(self) -> dict[str, str]:
        """
        Collect the info needed for writing to a decode results file

        Dictionary keys will be:
        - "bb_ids":
            the building block IDs for the compound, comma separated and ordered by cycle
        - "bb_scores":
            the building block call score for each call, comma separated and ordered by cycle
        - "umi":
            the UMI sequence

        Returns
        -------
        dict[str, str]
        """
        return {
            "bb_ids": ",".join(bb_call.obj.bb_id for bb_call in self.building_block_calls),
            "bb_scores": ",".join(str(bb_call.score) for bb_call in self.building_block_calls),
            "umi": self.umi,
        }

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
        return self.to_compound().enumerate().smi

    def get_library_id(self) -> str:
        """
        Get the library ID of the decoded compound

        Returns
        -------
        str
            the library ID
        """
        return self.library.library_id

    def get_score(self) -> float:
        """
        Get the decode score of the decoded compound

        Returns
        -------
        float
        """
        return sum(bb_call.score for bb_call in self.building_block_calls)


class DecodedToolCompound(DecodedCompound[ToolCompound]):
    """
    Holds information about a decoded barcode for a Tool compound

    Tool compounds are single compounds added to the DEL selection that exist outside the DELs

    Unlike Decoded DEL which has building block sections to call, a tool compound does not, since we know
    exactly what compound it is (there is only one option). Instead, it includes a score for how well the
    read aligned to the compound's tag, which is used to help determine if the call is correct

    Parameters
    ----------
    tool_compound: ToolCompound
        the decoded tool compound library
    alignment_score: float
        the score for how well the read aligned to the
    """

    def __init__(
        self,
        tool_compound: ToolCompound,
        alignment_score: float,
        umi: str,
    ):
        self.alignment_score = alignment_score
        self.tool_compound: ToolCompound = tool_compound
        super().__init__(umi=umi)

    def to_compound(self) -> ToolCompound:
        """
        Convert the decoded tool compound to a ToolCompound object

        Returns
        -------
        ToolCompound
        """
        return self.tool_compound

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

    def get_decode_res_info(self) -> dict[str, str]:
        """
        Collect the info needed for writing to a decode results file

        Since Tool compounds have no decodable sections, only will have the UMI key

        Dictionary keys will be:
        - "umi":
            the UMI sequence

        Returns
        -------
        dict[str, str]
        """
        return {"umi": self.umi}

    def get_smiles(self) -> str:
        """
        Get the SMILES of the decoded tool compound

        Returns
        -------
        str
            the SMILES of the compound

        Raises
        ------
        RuntimeError
            if the tool compound has no SMILES
        """
        if self.tool_compound.has_smiles():
            return self.tool_compound.smi
        else:
            raise RuntimeError(f"Single compound {self.tool_compound.compound_id} has no SMILES")

    def get_library_id(self) -> str:
        """
        Get the library ID of the decoded compound

        Returns
        -------
        str
            the library ID
        """
        return self.tool_compound.compound_id

    def get_score(self) -> float:
        """
        Get the decode score of the decoded compound

        Returns
        -------
        float
        """
        return self.alignment_score


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


class UMIContainsAmbiguity(FailedDecodeAttempt):
    """Returned if UMI contains ambiguous bases ('N')"""

    def __init__(self, sequence: SequenceRecord, umi_sequence: str):
        super().__init__(sequence=sequence, reason=f"UMI contains ambiguous bases: {umi_sequence}")


class UMIContainsINDEL(FailedDecodeAttempt):
    """Returned if UMI appears to have an INDEL"""

    def __init__(self, sequence: SequenceRecord, umi_sequence: str):
        super().__init__(sequence=sequence, reason=f"INDEL Detected in UMI: {umi_sequence}")


class BarcodeDecodingError(Exception):
    """raised when there is an error during barcode decoding"""

    pass


def _decode_section_first_greedy(
    codon: str, caller: BarcodeCaller, true_length: int
) -> tuple[ValidCall | FailedBarcodeLookup, int, int]:
    """
    Decode a section codon by greedily taking the first match found no matter the score

    This means a perfect match could be missed if the first match is imperfect but within tolerance.
    While faster on average, the gains from this are minimal compared to the best greedy approach,
    and the risk of missing a perfect match is high, so this method is not recommended in most cases.

    Notes
    -----
    This algorithm works by scanning the codon with various windows, one for each possible error type:
    - substitution / no error (true length window)
    - deletion (true length - 1 window)
    - insertion (true length + 1 window)

    All valid windows for the codon are generated in the order listed above and checked against the barcode caller.
    The caller will return either a ValidCall or a FailedBarcodeLookup. ValidCalls can be imperfect depending on
    the error tolerance of the caller. In this method, the first ValidCall found is returned immediately.

    Parameters
    ----------
    codon: str
        the codon sequence to decode
    caller: BarcodeCaller
        the barcode caller to use for decoding the codon windows
    true_length: int
        the expected length of the codon without any errors.
        This is used to determine the window sizes for scanning.

    Returns
    -------
    tuple[ValidCall | FailedBarcodeLookup, int, int]
        the decoded call (or failed lookup) and the start/stop indices of the matched window with respect to the codon
        if the lookup failed, start/stop indices are both -1
    """
    for scan_length in (true_length, true_length - 1, true_length + 1):  # sub, del, insert
        scan_start = 0
        scan_stop = scan_length
        while scan_stop <= len(codon):
            scan_codon = codon[scan_start:scan_stop]
            call = caller.decode_barcode(scan_codon)
            if not isinstance(call, FailedBarcodeLookup):
                return call, scan_start, scan_stop
            scan_start += 1
            scan_stop = scan_start + scan_length

    return FailedBarcodeLookup(codon), -1, -1


def _decode_section_best_greedy(
    codon: str, caller: BarcodeCaller, true_length: int, fail_on_conflict: bool = True
) -> tuple[ValidCall | FailedBarcodeLookup, int, int]:
    """
    Decode a section codon by greedily finding the best scoring call

    This approach will greedily search for the best scoring call across all possible windows.
    If a perfect match is found, it will return immediately as there is not a better possible call.
    Otherwise, it will return the first call with the best score across all possible windows.
    This means if more than one call has the same best score, the first one found will be returned, not
    both nor neither.

    Notes
    -----
    This algorithm works by scanning the codon with various windows, one for each possible error type:
    - substitution / no error (true length window)
    - deletion (true length - 1 window)
    - insertion (true length + 1 window)

    All valid windows for the codon are generated in the order listed above and checked against the barcode caller.
    The caller will return either a ValidCall or a FailedBarcodeLookup. ValidCalls can be imperfect depending on
    the error tolerance of the caller.

    Parameters
    ----------
    codon: str
        the codon sequence to decode
    caller: BarcodeCaller
        the barcode caller to use for decoding the codon windows
    true_length: int
        the expected length of the codon without any errors.
        This is used to determine the window sizes for scanning.

    Returns
    -------
    tuple[ValidCall | FailedBarcodeLookup, int, int]
        the decoded call (or failed lookup) and the start/stop indices of the matched window with respect to the codon
        if the lookup failed, start/stop indices are returned but meaningless
    """
    best_call: ValidCall | FailedBarcodeLookup = FailedBarcodeLookup(codon)
    best_score: float = 100.0  # large number
    best_start: int = -1
    best_stop: int = -1

    for scan_length in (true_length, true_length - 1, true_length + 1):  # sub, insert, del
        scan_start = 0
        scan_stop = scan_length
        while scan_stop <= len(codon):
            scan_codon = codon[scan_start:scan_stop]
            call = caller.decode_barcode(scan_codon)
            if not isinstance(call, FailedBarcodeLookup):
                if call.score == 0:  # perfect match found; can stop searching
                    return call, scan_start, scan_stop

                if fail_on_conflict and (best_score == call.score):
                    best_call = AmbiguousBarcodeCall(codon)
                else:
                    best_call = call
                    best_score = call.score
                    best_start = scan_start
                    best_stop = scan_stop
            scan_start += 1
            scan_stop = scan_start + scan_length

    # fail if no call found, else return the best one found
    return best_call, best_start, best_stop


def _decode_section_return_all(
    codon: str, caller: BarcodeCaller, true_length: int
) -> dict[int, list[tuple[ValidCall, int, int]]] | FailedBarcodeLookup:
    """
    Decode a section codon, return all possible valid calls sorted by score.

    This approach will scan all possible windows and return all valid calls found.
    Calls are split by their score into a dictionary where the key is the score.
    This is not a greedy approach, as it returns all possible calls, but it also is the slowest.
    Unless complex logic is used later to handle possible call conflicts between other sections,
    this method is not recommended.


    Notes
    -----
    This algorithm works by scanning the codon with various windows, one for each possible error type:
    - substitution / no error (true length window)
    - deletion (true length - 1 window)
    - insertion (true length + 1 window)

    All valid windows for the codon are generated in the order listed above and checked against the barcode caller.
    The caller will return either a ValidCall or a FailedBarcodeLookup. ValidCalls can be imperfect depending on
    the error tolerance of the caller.

    Parameters
    ----------
    codon: str
        the codon sequence to decode
    caller: BarcodeCaller
        the barcode caller to use for decoding the codon windows
    true_length: int
        the expected length of the codon without any errors.
        This is used to determine the window sizes for scanning.

    Returns
    -------
    dict[int, list[tuple[ValidCall, int, int]]] | FailedBarcodeLookup
        The valid calls found keyed by their score, with the start/stop indices of the matched window with
        respect to the codon.
        If the lookup failed, returns a FailedBarcodeLookup (not a dict)
    """
    calls: defaultdict[int, list[tuple[ValidCall, int, int]]] = defaultdict(list)

    for scan_length in (true_length, true_length - 1, true_length + 1):  # sub, insert, del
        scan_start = 0
        scan_stop = scan_length
        while scan_stop <= len(codon):
            scan_codon = codon[scan_start:scan_stop]
            call = caller.decode_barcode(scan_codon)
            if not isinstance(call, FailedBarcodeLookup):
                calls[int(call.score)].append((call, scan_start, scan_stop))
            scan_start += 1
            scan_stop = scan_start + scan_length

    # fail if no call found
    if len(calls) == 0:
        return FailedBarcodeLookup(codon)
    else:
        return calls


class _SequenceDecoder(abc.ABC):
    """Base interface for all sequence decoders"""

    @abc.abstractmethod
    def decode_sequence(self, aligned_sequence: AlignedSeq) -> DecodedCompound | FailedDecodeAttempt:
        """
        Given a sequence Alignment, decode it to determine the compound ID and UMI

        Parameters
        ----------
        aligned_sequence: AlignedSeq
            the aligned sequence to decode the building blocks from

        Returns
        -------
        DecodedCompound or FailedDecodeAttempt
            the decoded compound if successful,
            else a FailedDecodeAttempt explaining why it failed
        """
        raise NotImplementedError()


class DELibraryDecoder(_SequenceDecoder):
    """Base class for a library decoder"""

    def __init__(
        self,
        library: DELibrary,
        wiggle: bool = False,
        global_adjustment: bool = False,
        default_error_correction_mode_str="disable",
    ):
        self.library = library
        self.wiggle = wiggle
        self.global_adjustment = global_adjustment
        self.default_error_correction_mode_str = default_error_correction_mode_str

        self._callers: dict[BarcodeSection, BarcodeCaller[TaggedBuildingBlock]] = {
            bb_sec: get_barcode_caller(bb_set.tag_map, bb_sec.error_correction_mode)
            for bb_sec, bb_set in library.iter_bb_barcode_sections_and_sets()
        }

        _required_sections: MutableSequence[BarcodeSection] = list()

        # all building block sections + umi if present
        for sec in library.barcode_schema.building_block_sections + [library.barcode_schema.umi_section]:
            if sec is not None:
                _required_sections.append(sec)
        self._section_order = [
            bb_sec for bb_sec in library.barcode_schema.barcode_sections if bb_sec in _required_sections
        ]  # order of the sections

        # UMI location could end up moving depending on how the decoded sections are called. If they are determined
        # to have INDELs in them, the UMI could shift in either direction. To account for this, we need to know which
        # decodable section are closest to the UMI on either side and how far apart these sections are. Then at
        # decode time, we can adjust the UMI location based on these distances and the determined section spans.
        # This means the aligned UMI position is not actually used, rather the UMI position is completely determined
        # by the final locations of the decoded barcoded sections post successful decode.

        # distance between the start of the UMI section and the end of the closest decodable section that is before it
        _umi_sec_idx = self._section_order.index(library.barcode_schema.umi_section)
        self._section_before_umi: None | BarcodeSection = (
            None if _umi_sec_idx == 0 else self._section_order[_umi_sec_idx - 1]
        )
        self._section_dist_before_umi = (
            0
            if self._section_before_umi is None
            else library.barcode_schema.get_length_between_sections(
                library.barcode_schema.umi_section.section_name, self._section_before_umi.section_name
            )
        )
        # distance between the end of the UMI section and the start of the closest decodable section that is after it
        self._section_after_umi: None | BarcodeSection = (
            None if _umi_sec_idx == (len(self._section_order) - 1) else self._section_order[_umi_sec_idx + 1]
        )
        self._section_dist_after_umi = (
            0
            if self._section_after_umi is None
            else library.barcode_schema.get_length_between_sections(
                library.barcode_schema.umi_section.section_name, self._section_after_umi.section_name
            )
        )

    def decode_sequence(self, aligned_sequence: AlignedSeq) -> DecodedCompound | FailedDecodeAttempt:
        """
        Given a library call, decode the barcode sections

        Parameters
        ----------
        aligned_sequence: AlignedSeq
            the aligned sequence to decode the building blocks from

        Returns
        -------
        DecodedCompound or FailedDecodeAttempt
            the decoded compound if successful,
            else a FailedDecodeAttempt explaining why it failed
        """
        cur_position = 0  # track position to avoid overlaps
        global_adj = 0  # track global adjustment due to indels in previous sections

        called_sec_spans: dict[BarcodeSection, tuple[int, int]] = dict()
        called_secs: dict[BarcodeSection, ValidCall[TaggedBuildingBlock]] = dict()

        for section, caller in self._callers.items():
            start, stop = aligned_sequence.section_spans[section.section_name]

            # if using global adjustment, adjust start/stop positions
            if self.global_adjustment:
                start += global_adj
                stop += global_adj

            # adjust span in case of overlap with last section
            start = max(start, cur_position)  # move start back if overlapping
            stop = max(start + len(section.section_tag), stop)  # ensure at least full length of section tag

            if self.wiggle:  # wiggle increase the span by 1 on the right side
                stop += 1
            codon = aligned_sequence.sequence.sequence[start:stop]
            call, call_start, call_stop = _decode_section_best_greedy(codon, caller, len(section.section_tag))

            if isinstance(call, FailedBarcodeLookup):
                if isinstance(call, AmbiguousBarcodeCall):
                    return AmbiguousBuildingBlockBarcode(
                        sequence=aligned_sequence.sequence,
                        barcode=codon,
                        bb_section_name=section.section_name,
                    )
                else:
                    return FailedBuildingBlockCall(
                        sequence=aligned_sequence.sequence,
                        barcode=codon,
                        bb_section_name=section.section_name,
                    )

            # save the calls
            called_sec_spans[section] = (start + call_start, start + call_stop)
            called_secs[section] = call

            # update current position and global adj
            cur_position = start + call_stop
            # account for wiggle that added 1 by subbing it of the stop again
            global_adj += (call_stop - call_start) - (stop - start - self.wiggle)

        # call the UMI
        if (self._section_after_umi is not None) and (self._section_before_umi is not None):
            umi_start = called_sec_spans[self._section_before_umi][1] + self._section_dist_before_umi
            umi_stop = called_sec_spans[self._section_before_umi][0] - self._section_dist_after_umi
            if (umi_stop - umi_start) != len(self.library.barcode_schema.umi_section.section_tag):
                return UMIContainsINDEL(
                    sequence=aligned_sequence.sequence,
                    umi_sequence=aligned_sequence.sequence.sequence[umi_start:umi_stop],
                )
            umi_call = aligned_sequence.sequence.sequence[umi_start:umi_stop]
        elif self._section_before_umi is not None:
            umi_start = called_sec_spans[self._section_before_umi][1] + self._section_dist_before_umi
            umi_stop = umi_start + len(self.library.barcode_schema.umi_section.section_tag)
            umi_call = aligned_sequence.sequence.sequence[umi_start:umi_stop]
        elif self._section_after_umi is not None:
            umi_stop = called_sec_spans[self._section_after_umi][0] - self._section_dist_after_umi
            umi_start = umi_stop - len(self.library.barcode_schema.umi_section.section_tag)
            umi_call = aligned_sequence.sequence.sequence[umi_start:umi_stop]
        else:
            return FailedDecodeAttempt(aligned_sequence.sequence, "No sections to determine UMI location")

        if "N" in umi_call:
            return UMIContainsAmbiguity(
                sequence=aligned_sequence.sequence,
                umi_sequence=umi_call,
            )

        if all(bb_call.obj.is_real() for bb_call in called_secs.values()):  # real DEL compound
            return DecodedDELCompound(
                library=self.library,
                building_block_calls=[called_secs[sec] for sec in self.library.barcode_schema.building_block_sections],
                umi=umi_call,
            )
        elif all(bb_call.obj.is_real() for bb_call in called_secs.values()):  # mix of real and fake, impossible
            return ImpossibleBuildingBlockBarcode(
                sequence=aligned_sequence.sequence,
                bb_ids=[bb_call.obj.bb_id for bb_call in called_secs.values()],
                library_id=self.library.library_id,
            )
        else:  # all tool compound block calls
            tool_ids = [bb_call.obj.bb_id for bb_call in called_secs.values()]
            if len(set(tool_ids)) > 1:  # more than one tool compound found, impossible
                return ImpossibleBuildingBlockBarcode(
                    sequence=aligned_sequence.sequence,
                    bb_ids=[bb_call.obj.bb_id for bb_call in called_secs.values()],
                    library_id=self.library.library_id,
                )
            return DecodedToolCompound(
                tool_compound=self.library.get_tool_compound(tool_ids[0]),
                alignment_score=sum(bb_call.score for bb_call in called_secs.values()),
                umi=umi_call,
            )


class ToolCompoundDecoder(_SequenceDecoder):
    """
    Decoder for Tool Compounds

    Unlike libraries, tool compounds don't have decodable sections. Instead, they have a reference tag
    that is aligned to the read. If the alignment is good enough (based on a distance cutoff), the tool compound
    is considered called.

    Parameters
    ----------
    tagged_tool_compound: TaggedToolCompound
        the tagged tool compound to decode
    dist_cutoff: int, default = 2
        the maximum Levenshtein distance allowed between the read and the tool compound reference sequence
    """

    def __init__(self, tagged_tool_compound: TaggedToolCompound, dist_cutoff: int = 2):
        self.tagged_tool_compound = tagged_tool_compound
        self.dist_cutoff = dist_cutoff

    def decode_sequence(self, aligned_sequence: AlignedSeq) -> DecodedCompound | FailedDecodeAttempt:
        """
        Given a library call, decode the barcode sections

        Unlike other library decoders, tool compounds don't always have decodable sections. To validate that it is
        a correct call, they will look for a

        Parameters
        ----------
        aligned_sequence: AlignedSeq
            the aligned sequence to decode the building blocks from

        Returns
        -------
        DecodedCompound or FailedDecodeAttempt
            the decoded compound if successful,
            else a FailedDecodeAttempt explaining why it failed
        """
        # align the read to the tool compound tag
        ref_codon_start, ref_codon_stop = aligned_sequence.section_spans[
            self.tagged_tool_compound.barcode_schema.tool_compound_ref_section.section_name
        ]
        ref_codon = aligned_sequence.sequence.sequence[ref_codon_start:ref_codon_stop]
        real_codon = self.tagged_tool_compound.barcode_schema.tool_compound_ref_section.section_tag
        distance = levenshtein_distance(ref_codon, real_codon)

        if distance > self.dist_cutoff:
            return AlignmentFailed(aligned_sequence.sequence)
        else:
            return DecodedToolCompound(
                tool_compound=self.tagged_tool_compound,
                alignment_score=float(distance),
                umi=aligned_sequence.sequence.sequence[
                    aligned_sequence.section_spans["umi"][0] : len(
                        self.tagged_tool_compound.barcode_schema.umi_section.section_tag
                    )
                ],
            )


class DecodedCollectionCompound:
    """Holds the results of decoding a sequence from a DEL collection"""

    def __init__(self, library_call: ValidCall[DELibrary | ToolCompoundLibrary_], decoded_compound: DecodedCompound):
        self.library_call = library_call
        self.decoded_compound = decoded_compound

    def get_overall_score(self) -> float:
        """
        Get the overall score of the decoded collection compound

        This is the sum of the library call score and the decoded compound score

        Returns
        -------
        float
        """
        return self.library_call.score + self.decoded_compound.get_score()

    def to_decode_res_row_dict(self) -> dict[str, str]:
        """
        Convert the decoded collection compound to a dictionary of values for writing to a decode results file

        Dictionary keys will be:
        - "library_id":
            the library ID
        - "library_score":
            the score of the library call
        - "overall_score":
            the overall score of the decoded collection compound
        - plus all keys from the decoded compound's decode res row dict

        Returns
        -------
        dict[str, str]
        """
        row_dict: dict[str, str] = {
            "library_id": self.library_call.obj.library_id,
            "library_score": str(self.library_call.score),
            "overall_score": str(self.get_overall_score()),
        }
        row_dict.update(self.decoded_compound.get_decode_res_info())
        return row_dict

    def get_library_id(self) -> str:
        """
        Get the library ID of the decoded compound

        Returns
        -------
        str
            the library ID
        """
        return self.library_call.obj.library_id


class DELCollectionDecoder:
    """
    Decodes reads into DEL compounds

    Parameters
    ----------
    library_demultiplexer: LibraryDemultiplexer
        the library demultiplexer to use for calling libraries
    wiggle: bool, default = False
        If true, allow for wiggling the tag to find the best match
    global_adjustment: bool, default = False
        If true, adjust the positions of later sections based on indels found in earlier sections
        when decoding a read.
        Note: Should only be used when the alignment method is not flanking, since flanking alignments
        already account for indels in the barcode region.
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
        global_adjustment: bool = False,
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
        self.library_decoders: dict[Library, _SequenceDecoder] = dict()
        for library in self.library_demultiplexer.all_libraries:
            if isinstance(library, DELibrary):
                self.library_decoders[library] = DELibraryDecoder(
                    library=library,
                    wiggle=wiggle,
                    default_error_correction_mode_str=default_error_correction_mode_str,
                    global_adjustment=global_adjustment,
                )
            elif isinstance(library, ToolCompoundLibrary_):
                self.library_decoders[library] = ToolCompoundDecoder(
                    tagged_tool_compound=library.compound,
                )

        # set the min/max lengths
        self._min_read_length: int
        if isinstance(min_read_length, int):
            self._min_read_length = min_read_length
        else:
            self._min_read_length = (
                min([lib.barcode_schema.min_length for lib in library_demultiplexer.all_libraries]) - 10
            )  # allow some slack for seq errors

        self._max_read_length: int
        if isinstance(max_read_length, int):
            self._max_read_length = max_read_length
        else:
            self._max_read_length = 5 * self._min_read_length

    def decode_read(self, sequence: SequenceRecord) -> DecodedCollectionCompound | FailedDecodeAttempt:
        """
        Given a raw read, decode it to determine the compound it encodes

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
        DecodedCollectionCompound or FailedDecodeAttempt
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

            # loop through possible alignments until one decodes successfully, they run out, or max attempts reached
            decoded_compound: None | DecodedCompound | FailedDecodeAttempt = None
            for attempt, (seq_alignment, _) in enumerate(alignment_iter):
                if attempt >= MAX_RETRIES:
                    break
                decoded_compound = library_decoder.decode_sequence(seq_alignment)
                if isinstance(decoded_compound, DecodedCompound):  # quit once an alignment is successful
                    break

            if decoded_compound is None:  # no alignments were found, alignment failed
                decoded_compound = AlignmentFailed(sequence.sequence)

            # track the number of errors
            if isinstance(decoded_compound, FailedDecodeAttempt):
                if isinstance(decoded_compound, FailedBuildingBlockCall):
                    self.decode_statistics.num_failed_building_block_call += 1
                elif isinstance(decoded_compound, AmbiguousBuildingBlockBarcode):
                    self.decode_statistics.num_failed_ambiguous_building_block_call += 1
                elif isinstance(decoded_compound, (UMIContainsAmbiguity, UMIContainsINDEL)):
                    self.decode_statistics.num_failed_umi += 1
                elif isinstance(decoded_compound, ImpossibleBuildingBlockBarcode):
                    self.decode_statistics.num_failed_building_block_call += 1
                elif isinstance(decoded_compound, AlignmentFailed):
                    self.decode_statistics.num_failed_alignment += 1
                return decoded_compound
            else:
                self.decode_statistics.num_seqs_decoded_per_lib[decoded_compound.get_library_id()] += 1
                return DecodedCollectionCompound(
                    library_call=library_call,
                    decoded_compound=decoded_compound,
                )


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

    ignore_tool_compounds: bool = False
    demultiplexer_algorithm: Literal["cutadapt", "regex", "full"] = "regex"
    demultiplexer_mode: Literal["library", "single", "flanking"] = "single"
    realign: bool = False
    library_error_tolerance: int = 1
    library_error_correction_mode_str: str = "levenshtein_dist:2,asymmetrical"
    min_library_overlap: int = 8
    revcomp: bool = False
    wiggle: bool = False
    max_read_length: Optional[int] = None
    min_read_length: Optional[int] = None
    default_error_correction_mode_str: str = "levenshtein_dist:1,asymmetrical"

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

        demultiplexer: LibraryDemultiplexer = get_library_demultiplexer_type(
            demultiplex_mode=demultiplex_mode, demultiplex_algorithm=demultiplex_algorithm
        )(
            libraries=self.selection.library_collection,
            tool_compounds=list(self.selection.tool_compounds),
            **asdict(self.decode_settings),
        )

        # initialize all the decoding object required
        self.decoder = DELCollectionDecoder(
            library_demultiplexer=demultiplexer,
            decode_statistics=self.decode_stats,
            wiggle=self.decode_settings.wiggle,
            global_adjustment=self.decode_settings.demultiplexer_mode == "single"
            and (not self.decode_settings.realign),
            max_read_length=self.decode_settings.max_read_length,
            min_read_length=self.decode_settings.min_read_length,
            default_error_correction_mode_str=self.decode_settings.default_error_correction_mode_str,
        )

    def decode_read(self, read: SequenceRecord) -> DecodedCollectionCompound | FailedDecodeAttempt:
        """
        Decode a single read

        Parameters
        ----------
        read: SequenceRecord
            the read to decode

        Returns
        -------
        DecodedCollectionCompound or FailedDecodeAttempt
            The decoded compound if successful,
            else a `FailedDecodeAttempt` explaining why it failed.
        """
        self.decode_stats.num_seqs_read += 1
        decoded_compound = self.decoder.decode_read(read)
        if isinstance(decoded_compound, DecodedCompound):
            self.decode_stats.num_seqs_decoded_per_lib[decoded_compound.get_library_id()] += 1
        return decoded_compound

    def decode_file(
        self, fastq_file: os.PathLike, use_tqdm: bool = True
    ) -> Iterator[DecodedCollectionCompound | FailedDecodeAttempt]:
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
        from tqdm import tqdm

        from deli.dna.io import SingleFileSequenceReader

        reader = SingleFileSequenceReader(fastq_file)
        for sequence in tqdm(reader.iter_seqs(), disable=not use_tqdm, desc=f"Decoding file: {fastq_file}"):
            yield self.decode_read(sequence)
