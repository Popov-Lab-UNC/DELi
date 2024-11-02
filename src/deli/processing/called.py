"""code for calling del matches to del compounds"""

import enum
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Union

from Levenshtein import distance as levenshtein_distance

from deli.dels import (
    BarcodeSchema,
    BuildingBlock,
    DELibrary,
    DELibrarySchemaGroup,
    Index,
    IndexSet,
    Umi,
    get_min_index_distance,
    get_min_library_tag_distance,
)
from deli.sequence import SemiGlobalAlignment

from ..dels.barcode import BarcodeSection
from .match import BarcodeMatch


class BarcodeCallingError(Exception):
    """Error for issues that arise due to calling setting conflicts"""

    pass


class BarcodeAlignment(dict):
    """a wrapper of dict for alignments to help with readability"""

    pass


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


def _align_call_mode(match: BarcodeMatch, barcode_schema: BarcodeSchema) -> BarcodeAlignment:
    """
    Call mode function used to generate section alignment in align mode

    Notes
    -----
    Can use Fast or Full matches; just needs a truncated match region

    Parameters
    ----------
    match: FullBarcodeMatch
        the match to generate the alignment for
    barcode_schema: BarcodeSchema
        the schema defining the barcode to align to

    Returns
    -------
    BarcodeAlignment
    """
    seq_alignment = SemiGlobalAlignment.globally_align(
        seq1=barcode_schema.full_barcode,
        seq2=match.match_sequence.sequence,
    )

    _spans = {}
    for _section_name, (
        _section_start,
        _section_stop,
    ) in barcode_schema.barcode_spans.copy().items():
        _spans[_section_name] = (
            seq_alignment.alignment_lookup[_section_start],
            seq_alignment.alignment_lookup[_section_stop],
        )

    return BarcodeAlignment(_spans)


class CallModes(enum.Enum):
    """Enum for calling mode settings"""

    # functions need to be manually added as members to the enum
    ALIGN = enum.member(_align_call_mode)

    def __call__(self, *args, **kwargs):
        """Makes enum values called (since Enum only holds functions)"""
        return self.value(*args, **kwargs)


@dataclass
class SectionCall:
    """base class for all calls"""

    score: float = field(default=0.0, kw_only=True)


@dataclass
class FailedCall(SectionCall):
    """class for a failed call"""

    section_name: str = ""
    score: float = field(default=float("-inf"), kw_only=True)

    def __str__(self):
        """Always returns 'FAILED'"""
        return "FAILED"

    def to_dict(self):
        """Convert into a csv row dict"""
        return {self.section_name: "FAILED"}


@dataclass
class IndexCall(SectionCall):
    """Class to hold info of successful index call"""

    called_index: Index

    def __str__(self):
        """Return the called indexes id as a str"""
        return self.called_index.index_id

    def to_dict(self):
        """Convert into a csv row dict"""
        return {"index": self.called_index.index_id}


@dataclass
class LibraryCall(SectionCall):
    """Class to hold info of successful library call"""

    called_library: DELibrary

    def __str__(self):
        """Returns the called libraries id as a str"""
        return self.called_library.library_id

    def to_dict(self):
        """Convert into a csv row dict"""
        return {"library": self.called_library.library_id}


@dataclass
class UmiCall(SectionCall):
    """Class to hold info of successful umi call"""

    called_umi: Umi

    def __str__(self):
        """Returns the called umi as a str"""
        return self.called_umi.umi

    def to_dict(self):
        """Convert into a csv row dict"""
        return {"umi": self.called_umi.umi}


@dataclass
class BBCall(SectionCall):
    """Class to hold info of successful building block cycle calls"""

    bb_cycle: str
    bb_call: BuildingBlock

    def __str__(self):
        """Returns the called blocks id as a str"""
        return self.bb_call.bb_id

    def to_dict(self):
        """Convert into a csv row dict"""
        return {self.bb_cycle: self.bb_call.bb_id}


class CalledBarcode:
    """class for containing info on a called barcode"""

    def __init__(
        self,
        parent_match: BarcodeMatch,
        library_call: Union[LibraryCall, FailedCall],
        bb_calls: List[Union[BBCall, FailedCall]],
        umi_call: Union[UmiCall, FailedCall],
        index_call: Optional[Union[IndexCall, FailedCall]] = None,
    ):
        """
        Initialize a CalledBarcode object

        Parameters
        ----------
        parent_match: BarcodeMatch
            the match the call come from
        library_call: Union[LibraryCall, FailedCall]
            the called library (can be failed call)
        bb_calls: List[Union[BBCall, FailedCall]]
            the building block calls as a list (can be failed)
        umi_call: Union[UmiCall, FailedCall]
            the UMI call (can be failed)
        index_call: Optional[Union[IndexCall, FailedCall]]
            the index call
            leave as `None` if not de-multiplexing
        """
        self.parent_match = parent_match
        self.library_call = library_call
        self.bb_calls = bb_calls
        self.umi_call = umi_call
        self.index_call = index_call

    def get_all_calls(self) -> List[SectionCall]:
        """Return all the call object from the called barcode"""
        return (
            [self.library_call]
            + self.bb_calls
            + [self.umi_call]
            + ([self.index_call] if self.index_call else [])
        )

    def called_successfully(self) -> bool:
        """If all barcode sections called successfully return true"""
        return 0 <= self.get_overall_score()

    def get_overall_score(self) -> float:
        """Return the overall score of all compounds"""
        return sum([_call.score for _call in self.get_all_calls()])

    def to_row_dict(self) -> Dict[str, str]:
        """Convert the called barcode to csv row for easy saving"""
        _csv_dict = {
            "from_read": self.parent_match.match_sequence.read_id,
            "from_match": self.parent_match.match_id,
            "match_seq": self.parent_match.sequence,
        }

        _csv_dict.update(self.library_call.to_dict())

        # add the bb calls
        for bb_call in self.bb_calls:
            _csv_dict.update(bb_call.to_dict())

        # add DEL_ID if the call was good, else set it to failed
        if self.called_successfully():
            _csv_dict["DEL_ID"] = (
                str(self.library_call)
                + "_"
                + "_".join([str(bb_call) for bb_call in self.bb_calls])
            )
        else:
            _csv_dict["DEL_ID"] = "FAILED"

        # if an index is called add the index
        if self.index_call is not None:
            _csv_dict.update(self.index_call.to_dict())

        return _csv_dict


class BarcodeCaller:
    """Used to call matches into decoded DELs"""

    def __init__(
        self,
        libraries: DELibrarySchemaGroup,
        indexes: Optional[Union[List[Index], IndexSet]] = None,
        call_mode: CallModes = CallModes.ALIGN,
    ):
        """
        Initialize a DEL BarcodeCaller object

        Parameters
        ----------
        libraries: Union[List[DELibrary], MegaDELibrary]
            the possible libraries that can be called from
            can be a list of libraries or a mega library
            must all use teh same barcode schema as the passed
             `barcode_schema`
        indexes: Optional[Union[List[Index], IndexSet]] = IndexSet([])
            the possible indexes that can be called from
            can be a list of indexes or an index set
            if not de-multiplexing can be ignored
        call_mode: CallModes, default = CallModes.ALIGN
            the methods used for calling
            must be a valid CallModes function
        """
        # handle no indexes
        if indexes is None:
            indexes = IndexSet([])

        self.indexes = indexes if isinstance(indexes, IndexSet) else IndexSet(indexes)
        self.libraries = libraries
        self.call_mode = call_mode

        # precompute calling distance cutoffs
        self._max_library_dist: float = get_min_library_tag_distance(self.libraries) / 2
        self._max_index_dist: float = get_min_index_distance(self.indexes) / 2
        self._skip_calling_index: bool = len(self.indexes) <= 1
        self._skip_calling_lib: bool = not self.libraries.requires_multistep_calling

    @staticmethod
    def _make_call(query: str, refs: List[str], dist_cutoff: float) -> Tuple[int, float]:
        """Help func for repeated calls"""
        # return -1 if the call fails
        called: int = -1
        called_dist: float = -1.0
        for i, ref in enumerate(refs):
            dist = levenshtein_distance(ref, query)
            if dist <= dist_cutoff:
                if called != -1:
                    return -1, -1.0
                called_dist = dist
                called = i
        return called, called_dist

    @staticmethod
    def _decode_seq(seq: str, section: BarcodeSection) -> Tuple[str, bool]:
        if section.decoder:
            return section.decoder.decode_sequence(seq)
        return seq, True

    def _call_index(
        self, match: BarcodeMatch, alignment: BarcodeAlignment
    ) -> Union[FailedCall, IndexCall]:
        """Call the index"""
        # if skip index call the first index always
        if self._skip_calling_index:
            return IndexCall(self.indexes[0])

        aligned_index_seq = match.match_sequence.sequence[slice(*alignment["index"])]
        aligned_index_seq, _ = self._decode_seq(
            aligned_index_seq, self.libraries[0].barcode_schema["index"]
        )

        best_index_idx, best_index_dist = self._make_call(
            query=aligned_index_seq,
            refs=[index.dna_tag for index in self.indexes],
            dist_cutoff=self._max_index_dist,
        )

        if best_index_idx == -1:
            return FailedCall(section_name="index")
        else:
            return IndexCall(self.indexes[best_index_idx], score=best_index_dist)

    def _call_library(
        self, match: BarcodeMatch, alignment: BarcodeAlignment
    ) -> Union[FailedCall, LibraryCall]:
        """Call the library"""
        # if skip library call the first library always
        if self._skip_calling_lib:
            return LibraryCall(self.libraries[0])

        aligned_lib_seq = match.match_sequence.sequence[slice(*alignment["library"])]
        aligned_lib_seq, _ = self._decode_seq(
            aligned_lib_seq, self.libraries[0].barcode_schema["library"]
        )

        best_lib_idx, best_lib_dist = self._make_call(
            query=aligned_lib_seq,
            refs=[library.library_tag for library in self.libraries],
            dist_cutoff=self._max_library_dist,
        )

        if best_lib_idx == -1:
            return FailedCall(section_name="library")
        else:
            return LibraryCall(self.libraries[best_lib_idx], score=best_lib_dist)

    @staticmethod
    def _call_umi(match: BarcodeMatch, alignment: BarcodeAlignment) -> Union[UmiCall, FailedCall]:
        """Call the umi"""
        if "umi" in alignment.keys():
            return UmiCall(Umi(match.match_sequence.sequence[slice(*alignment["umi"])]))
        else:
            return FailedCall(section_name="umi")

    def _call_bb(
        self,
        match: BarcodeMatch,
        alignment: BarcodeAlignment,
        called_lib: Union[LibraryCall, FailedCall],
    ) -> List[Union[FailedCall, BBCall]]:
        """Call the building blocks"""
        if isinstance(called_lib, FailedCall):
            return [FailedCall(section_name="bb1")]

        _lib = called_lib.called_library
        bb_regions = _lib.barcode_schema.get_bb_regions()

        bb_calls: List[Union[BBCall, FailedCall]] = []

        for bb_region, bb_set in zip(bb_regions, _lib.iter_bb_sets()):
            aligned_bb_seq = match.match_sequence.sequence[slice(*alignment[bb_region])]
            aligned_bb_seq, _passed_decode = self._decode_seq(
                aligned_bb_seq, _lib.barcode_schema[bb_region]
            )

            if not _passed_decode:
                bb_calls.append(FailedCall(section_name=bb_region))

            bb_call = bb_set.search_tags(aligned_bb_seq)
            if isinstance(bb_call, BuildingBlock):
                bb_calls.append(BBCall(bb_region, bb_call))
            else:
                bb_calls.append(FailedCall(section_name=bb_region))

        return bb_calls

    def call_tag(self, match: BarcodeMatch) -> CalledBarcode:
        """
        Call/decode a match into the DEL compound the DNA encodes

        Notes
        -----
        Even if some parts of the calling fails (e.i. no match for bb bases)
        a CalledBarcode is still returned.
        CalledBarcodes have functions that can assess if the call
        was successful

        Parameters
        ----------
        match: BarcodeMatch
            the barcode match to call

        Returns
        -------
        CalledBarcode
        """
        # call index/library first if multiple libraries
        if self._skip_calling_lib:
            # get global alignment
            alignment = self.call_mode(match, self.libraries[0].barcode_schema)

            # call the library bases
            library_call = self._call_library(match, alignment)
        else:
            lib_alignment = self.call_mode(match, self.libraries.library_call_schema)
            # call the library bases
            library_call = self._call_library(match, lib_alignment)

            if isinstance(library_call, FailedCall):
                return CalledBarcode(
                    parent_match=match,
                    index_call=FailedCall(),
                    library_call=library_call,
                    umi_call=FailedCall(),
                    bb_calls=[FailedCall(section_name="bb1")],
                )

            alignment = self.call_mode(match, library_call.called_library.barcode_schema)

        # call index if it is needed
        if len(self.indexes) != 0:
            index_call = self._call_index(match, alignment)
        else:
            index_call = None

        bb_calls = self._call_bb(match, alignment, library_call)
        umi_call = self._call_umi(match, alignment)

        return CalledBarcode(
            parent_match=match,
            index_call=index_call,
            library_call=library_call,
            umi_call=umi_call,
            bb_calls=bb_calls,
        )

    def call_tags(self, matches: List[BarcodeMatch]) -> List[CalledBarcode]:
        """
        calls a list of barcode matches

        Notes
        -----
        See `BarcodeMatch` for more details on calling

        See Also
        --------
        BarcodeMatch

        Parameters
        ----------
        matches: List[BarcodeMatch]
            matches to call

        Returns
        -------
        List[CalledBarcode]
        """
        return [self.call_tag(match) for match in matches]
