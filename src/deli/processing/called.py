import enum
import warnings
from typing import Literal, Union, List, Optional

from Levenshtein import distance as levenshtein_distance

from deli.dels import BarcodeSchema, Index, DELibrary, get_min_index_distance
from deli.sequence import SemiGlobalAlignment

from .match import FastBarcodeMatch, FullBarcodeMatch


class BarcodeCallingError(Exception):
    pass


class BarcodeAlignment(dict):
    pass


def _between(
        val: int,
        start: int,
        stop: int,
        right_inclusive: bool = True,
        left_inclusive: bool = False):
    """
    Helper function for determining if some value is between two values

    Parameters
    ----------
    val: int
        value to check if in between
    start: int
        start of the between range
    stop: int
        stop of the between range
    right_inclusive: bool = True
        include the `start` value in the range
    left_inclusive: bool = True
        include the `stop` value in the range

    Returns
    -------
    bool
        if `val` is between `start` and `stop`
    """
    return (start - right_inclusive) < val < (stop + left_inclusive)


def _match_call_mode(
        match: FullBarcodeMatch,
        barcode_schema: BarcodeSchema
) -> BarcodeAlignment:
    """
    Call mode function used to generate section alignment in match mode

    Notes
    -----
    Needs FullBarcodeMatches (generated with the `full_match=True` flag
    during the matching)

    Parameters
    ----------
    match: FullBarcodeMatch
        the match to generate the alignment for
    barcode_schema:
        the schema defining the barcode to align to

    Returns
    -------
    BarcodeAlignment
    """
    _tmp_alignment = barcode_schema.barcode_spans.copy()

    _errors = list(zip(*match.get_sorted_match_sequence_error_idx(exclude_substitute=True)))
    if len(_errors) == 0:
        return BarcodeAlignment(_tmp_alignment)

    (_err_idx, _err_type) = _errors.pop(0)
    _buff = 0
    _spans = barcode_schema.barcode_spans.copy()
    for _key, (_start, _stop) in _tmp_alignment.items():
        _start = _start + _buff
        _stop = _stop + _buff
        _end = _stop
        while _err_idx < _end:
            if _between(_err_idx, _start, _stop):
                if _err_type == "insert":
                    _stop += 1
                elif _err_type == "delete":
                    _stop -= 1
                else:
                    raise ValueError(f"unrecognized error type {_err_type}")

                try:
                    (_err_idx, _err_type) = _errors.pop(0)
                except IndexError:
                    (_err_idx, _err_type) = (1e9, "no-more-errors")
        _spans[_key] = (_start, _stop)
        _buff += (_stop - _end)
    return BarcodeAlignment(_spans)


def _align_call_mode(
        match: Union[FastBarcodeMatch, FullBarcodeMatch],
        barcode_schema: BarcodeSchema
) -> BarcodeAlignment:
    """
    Call mode function used to generate section alignment in align mode

    Notes
    -----
    Can use Fast or Full matches; just needs a truncated match region

    Parameters
    ----------
    match: FullBarcodeMatch
        the match to generate the alignment for
    barcode_schema:
        the schema defining the barcode to align to

    Returns
    -------
    BarcodeAlignment
    """
    seq_alignment = SemiGlobalAlignment.globally_align(
        seq1=barcode_schema.full_barcode,
        seq2=match.match_sequence.sequence,
    )

    # start only after the first non-gap in the barcode
    _start = next(i for i, item in enumerate(seq_alignment) if item[0] is not None)
    _match_start = seq_alignment.alignment[_start][1]

    barcode_sections = list(barcode_schema.barcode_spans.copy().items())
    _section_name, (_section_start, _section_stop) = barcode_sections.pop(0)

    _spans = {}
    _last_non_null = None
    for barcode_idx, match_idx in seq_alignment.alignment[_start:]:
        if barcode_idx == _section_stop:
            _spans[_section_name] = (_match_start, _last_non_null + 1)
            _match_start = _last_non_null + 1
            try:
                _section_name, (_section_start, _section_stop) = barcode_sections.pop(0)
            except IndexError:
                break
        if match_idx is not None:
            _last_non_null = match_idx
    return BarcodeAlignment(_spans)


class CallModes(enum.Enum):
    """Enum for calling mode settings"""
    MATCH = _match_call_mode
    ALIGN = _align_call_mode


class CallOutcome:
    def __init__(
            self,
    ):
        pass


class BarcodeCaller:
    def __init__(
            self,
            barcode_schema: BarcodeSchema,
            libraries: List[DELibrary],
            indexes: Optional[List[Index]] = None,
            call_mode: CallModes = CallModes.ALIGN,
    ):
        self.barcode_schema = barcode_schema
        self.indexes = indexes
        self.libraries = libraries
        self.call_mode = call_mode

        if self.indexes:
            self.max_index_dist = get_min_index_distance(self.indexes)
        else:
            self.max_index_dist = None

        self.max_index_dist = get_min_index_distance(self.indexes)

    def _align(self, match: Union[FastBarcodeMatch, FullBarcodeMatch]) -> BarcodeAlignment:
        """wraps the align function from the call mode enum"""
        return self.call_mode.value(match, self.barcode_schema)

    def _call_index(self, match: Union[FastBarcodeMatch, FullBarcodeMatch], alignment: BarcodeAlignment) -> Index:
        # this error should never be seen by a user; here for development
        if self.indexes is None:
            raise RuntimeError("attempted to call index when not index was passed")

        if self.indexes and ("index" not in alignment.keys()):
            raise ValueError("barcode schema lacks index section; index cannot be called")

        matched_index = match.match_sequence.sequence[alignment["index"]]




# class CalledSequence:
#     """
#     This class is used to take a SequenceMatch and attempt to assign it a INDEX_ID, LIB_ID and DEL_ID for
#     downstream analysis and processing
#
#     Parameters
#     ----------
#     sequence_match : SequenceMatch
#         the SequenceMatch object to attempt DEL calling on
#     make: DELMake
#         the DELMake object defining the make up of the barcode scheme used for this read
#     mode: Literal["low-quality", "normal", "high-quality"], default = "high-quality"
#         which level of processing to use when trying to make the call
#     strict: bool, default = False
#         use strict mode or not
#
#     Attributes
#     ----------
#     barcode_pat : dict
#         dictionary of the (inclusive) start and (non-inclusive) stop index of each section: `dict[key] = [start, stop)`
#     match_risk: float
#         the numerical risk associated with the SequenceMatch
#     seq: str
#         the raw DNA sequence of the SequenceMatch
#     index: str
#         the raw DNA sequence of the index section in the barcode
#     primer: str
#         the raw DNA sequence of the primer section in the barcode
#     lib_barcode: str
#         the raw DNA sequence of the library barcode section in the barcode
#     lib_barcode_toe: str
#         the raw DNA sequence of the library barcode toe section in the barcode
#     bba: str
#         the raw DNA sequence of the bba section in the barcode
#     bba_toe: str
#         the raw DNA sequence of the bba toe section in the barcode
#     bbb: str
#         the raw DNA sequence of the bbb section in the barcode
#     bbb_toe: str
#         the raw DNA sequence of the bbb toe section in the barcode
#     bbc: str
#         the raw DNA sequence of the bbc section in the barcode
#     bbc_toe: str
#         the raw DNA sequence of the bbc toe section in the barcode
#     closing_primer: str
#         the raw DNA sequence of the closing primer section in the barcode
#     umi: str
#         the raw DNA sequence of the umi section in the barcode
#     del_id: str
#         the called DEL ID for SquenceMatch
#         will be "FAILED" if a DEL id could not be called (either BB lookup failed of library failed to be called)
#     lib_toe_risk: float
#         numerical risk associated with the library toe region
#     bba_call: str
#         the BB ID for the BB in cycle A
#         if `called_lib` is "FAILED", will be "FAILED"
#         if a substitution mutation is detected as the reason the lookup failed will be "failed_snp"
#         if a insertion mutation is detected as the reason the lookup failed will be "failed_insert"
#         if a deletion mutation is detected as the reason the lookup failed will be "failed_delete"
#     bbb_call: str
#         the BB ID for the BB in cycle B
#         if `called_lib` is "FAILED", will be "FAILED"
#         if a substitution mutation is detected as the reason the lookup failed will be "failed_snp"
#         if a insertion mutation is detected as the reason the lookup failed will be "failed_insert"
#         if a deletion mutation is detected as the reason the lookup failed will be "failed_delete"
#     bbc_call: str
#         the BB ID for the BB in cycle C
#         if `called_lib` is "FAILED", will be "FAILED"
#         if a substitution mutation is detected as the reason the lookup failed will be "failed_snp"
#         if a insertion mutation is detected as the reason the lookup failed will be "failed_insert"
#         if a deletion mutation is detected as the reason the lookup failed will be "failed_delete"
#     bba_toe_risk: float
#         numerical risk associated with the bba toe region
#     bbb_toe_risk: float
#         numerical risk associated with the bbb toe region
#     bbc_toe_risk: float
#         numerical risk associated with the bbc toe region
#     called_index_risk: float
#         numerical risk associated with the index call that was made
#     called_index: str
#         the index call that was made
#         will be "FAILED" if an index could not be called
#     called_lib_risk: float
#         numerical risk associated with the library call that was made
#     called_lib: str
#         the library call that was made
#         will be "FAILED" if an library could not be called
#
#     Notes
#     -----
#     It is almost always recommended to use 'high-quality' mode.
#     Only use the others if you know exactly what they do and why you want to run them
#
#     'normal' and 'low-quality' mode will only work if the matching mode used to generate the SequenceMatch was
#     'full' or 'minimum', otherwise it will fail and possibly produce lots of wrong calls
#     """
#     def __init__(self, sequence_match: SequenceMatch, make: DELMake,
#                  mode: Literal["low-quality", "normal", "high-quality"] = "high-quality",
#                  strict: bool = False):
#         self.sequence_match = sequence_match
#         self.mode = mode
#         self.make = make
#
#         if strict:
#             self._alignment_seq = make["BARCODE_FULL_SEQUENCE"]
#         else:
#             self._alignment_seq = make["BARCODE_SEQUENCE"]
#
#         # get the index mapping to each region
#         self.barcode_pat = {}
#         self._get_alignment()
#
#         # attributes used to assess call quality
#         self.match_risk = self.sequence_match.total_errors()
#         self.seq = self.sequence_match.match_sequence
#         self.index = self.seq[self.barcode_pat["index"][0]:self.barcode_pat["index"][1]]
#         self.primer = self.seq[self.barcode_pat["primer"][0]:self.barcode_pat["primer"][1]]
#         self.lib_barcode = self.seq[self.barcode_pat["lib_barcode"][0]:self.barcode_pat["lib_barcode"][1]]
#         self.lib_barcode_toe = self.seq[self.barcode_pat["lib_barcode_toe"][0]:self.barcode_pat["lib_barcode_toe"][1]]
#         self.bba = self.seq[self.barcode_pat["bb1"][0]:self.barcode_pat["bb1"][1]]
#         self.bba_toe = self.seq[self.barcode_pat["bb1_toe"][0]:self.barcode_pat["bb1_toe"][1]]
#         self.bbb = self.seq[self.barcode_pat["bb2"][0]:self.barcode_pat["bb2"][1]]
#         self.bbb_toe = self.seq[self.barcode_pat["bb2_toe"][0]:self.barcode_pat["bb2_toe"][1]]
#         self.bbc = self.seq[self.barcode_pat["bb3"][0]:self.barcode_pat["bb3"][1]]
#         self.bbc_toe = self.seq[self.barcode_pat["bb3_toe"][0]:self.barcode_pat["bb3_toe"][1]]
#         self.closing_primer = self.seq[self.barcode_pat["closing_primer"][0]:self.barcode_pat["closing_primer"][1]]
#         self.umi = self.seq[self.barcode_pat["umi"][0]:self.barcode_pat["umi"][1]]
#
#         # risk and call attributes
#         self.del_id = None
#         self.lib_toe_risk = None
#         self.bba_call = None
#         self.bbb_call = None
#         self.bbc_call = None
#         self.bba_toe_risk = None
#         self.bbb_toe_risk = None
#         self.bbc_toe_risk = None
#
#         self.called_index_risk = None
#         self.called_index = None
#         self.called_lib_risk = None
#         self.called_lib = None
#
#
#     def call_bbs(self):
#         """
#         Attempts to call the BB IDs for the DEL
#
#         In order to function, needs three json files, one for each A-B-C cycle that defines a map from sequence to
#         BB_ID.
#
#         The above files should be stored in `./data/bb_codes/<CALLED_LIB>/<CALLED_LIB_BB<A,B,C>.json` where CALLED_LIB
#         is the name of the library that was called.
#
#         If the called library lacks existing files in the above locations, the called BB IDs will all be "FAILED"
#
#         If the lookup fails, will use the alignment to detect (with moderate accuracy) the sequencing error that
#         occurred resulting the failed lookup.
#         If `called_lib` is "FAILED" will be "FAILED"
#         If a substitution mutation is detected as the reason the lookup failed will be "failed_snp"
#         If an insertion mutation is detected as the reason the lookup failed will be "failed_insert"
#         If a deletion mutation is detected as the reason the lookup failed will be "failed_delete"
#
#         Notes
#         -----
#         Must be called after the library is called with `self.call_lib()`
#         """
#         self.bba_call = "FAILED"
#         self.bbb_call = "FAILED"
#         self.bbc_call = "FAILED"
#
#         if self.called_lib is None:
#             raise ValueError("cannot call bbs without calling lib first")
#         if self.called_lib == "FAILED" or self.called_lib not in self.make["BB_CODES_NO_TOE"].keys():
#             self.del_id = "FAILED"
#             self.lib_toe_risk = 100
#             self.bba_toe_risk = 100
#             self.bbb_toe_risk = 100
#             self.bbc_toe_risk = 100
#             return
#
#         bba_length = len(self.bba)
#         bbb_length = len(self.bbb)
#         bbc_length = len(self.bbc)
#
#         if bba_length > self.make["BB_LENGTHS"][self.called_lib]:
#             self.bba_call = "failed_insert"
#         if bbb_length > self.make["BB_LENGTHS"][self.called_lib]:
#             self.bbb_call = "failed_insert"
#         if bbc_length > self.make["BB_LENGTHS"][self.called_lib]:
#             self.bbc_call = "failed_insert"
#
#         if bba_length < self.make["BB_LENGTHS"][self.called_lib]:
#             self.bba_call = "failed_delete"
#         if bbb_length < self.make["BB_LENGTHS"][self.called_lib]:
#             self.bbb_call = "failed_delete"
#         if bbc_length < self.make["BB_LENGTHS"][self.called_lib]:
#             self.bbc_call = "failed_delete"
#
#         _a = self.make["BB_CODES_NO_TOE"][self.called_lib][self.make["CYCLE_A_NAME"]].get(self.bba)
#         _b = self.make["BB_CODES_NO_TOE"][self.called_lib][self.make["CYCLE_B_NAME"]].get(self.bbb)
#         _c = self.make["BB_CODES_NO_TOE"][self.called_lib][self.make["CYCLE_C_NAME"]].get(self.bbc)
#
#         if self.bba_call == "FAILED":
#             if _a is None:
#                 self.bba_call = "failed_snp"
#             else:
#                 self.bba_call = _a
#
#         if self.bbb_call == "FAILED":
#             if _b is None:
#                 self.bbb_call = "failed_snp"
#             else:
#                 self.bbb_call = _b
#
#         if self.bbc_call == "FAILED":
#             if _c is None:
#                 self.bbc_call = "failed_snp"
#             else:
#                 self.bbc_call = _c
#
#         if all([_a, _b, _c]):
#             self.del_id = "-".join([_a, _b, _c])
#         else:
#             self.del_id = "FAILED"
#
#         self.lib_toe_risk = levenshtein_distance(self.make["DEL_BARCODE_SECTIONS"]["lib_barcode_toe"], self.lib_barcode_toe)
#         self.bba_toe_risk = levenshtein_distance(self.make["DEL_BARCODE_SECTIONS"]["bb1_toe"], self.bba_toe)
#         self.bbb_toe_risk = levenshtein_distance(self.make["DEL_BARCODE_SECTIONS"]["bb2_toe"], self.bbb_toe)
#         self.bbc_toe_risk = levenshtein_distance(self.make["DEL_BARCODE_SECTIONS"]["bb3_toe"], self.bbc_toe)
#
#     def call_index(self, included_index: list[str] = None, min_dist: int = None):
#         """
#         Will attempt to use the alignment to call the index (experiment ID) of the SequenceMatch.
#         Will only look for index that are passed in `included_index`.
#         The index sequences are stored in the `/data/experiment_index.json` file.
#
#         When mode is "high-quality", the index will be called successfully if the Levenshtein distance between the
#         index DNA sequence from the match and a given index sequence is less than the minimum distance between all
#         indexes in `included_index`
#
#         When mode is "normal" the index will be called successfully if the smallest Levenshtein distance between the
#         index DNA sequence from the match and a given index is less than `INDEX_RISK_NORMAL_THRESHOLD` compared to
#         the next smallest distance to another index
#
#         When mode is "low-quality" the index will be called based on the index sequence with the lowest Levenshtein
#         distance to the matched index DNA sequence
#
#         `called_index` will be "FAILED" if an index cannot be assigned using the respective method
#
#         Parameters
#         ----------
#         included_index: list[str], default = `None`
#             the list of indexes names that should be searched for
#             index names should exist with a corresponding sequence in the `/data/experiment_index.json` file.
#         min_dist: int, default = None
#             the min Levenshtein distance between the passed index
#
#         Notes
#         -----
#         "low-quality" mode can never fail to call an index
#         If left as `None` included_index will default to all indexes in the config files
#         """
#         distances = get_dist_from_index(self.index, included_index)
#         self.called_index = "FAILED"
#         self.called_index_risk = 100
#
#         if self.mode == "high-quality":
#             if min_dist is None:
#                 _dist_thresh = get_min_index_distance(included_index) / 2
#             else:
#                 _dist_thresh = min_dist / 2
#             for _index_key, dist in distances.items():
#                 if dist < _dist_thresh:
#                     self.called_index = _index_key
#                     self.called_index_risk = dist
#                     break
#                 if dist < self.called_index_risk:
#                     self.called_index_risk = dist
#         elif self.mode == "normal":
#             _curr_dist = 100
#             for _index_key, dist in distances.items():
#                 if dist < _curr_dist:
#                     if (_curr_dist - dist) < INDEX_RISK_NORMAL_THRESHOLD:
#                         self.called_index = "FAILED"
#                     else:
#                         self.called_index = _index_key
#                     _curr_dist = dist
#             self.called_index_risk = _curr_dist
#         else:
#             _curr_dist = 100
#             _curr_fail = False
#             for _index_key, dist in distances.items():
#                 if dist < _curr_dist:
#                     self.called_index = _index_key
#                     _curr_dist = dist
#             self.called_index_risk = _curr_dist
#
#     def call_lib(self, included_libs: list[str] = None, min_dist: int = None):
#         """
#         Will attempt to use the alignment to call the library ID of the SequenceMatch.
#         Will only look for index that are passed in `included_libs`.
#         The index sequences are stored in the `/data/libs.json` file.
#
#         When mode is "high-quality", the library will be called successfully if the Levenshtein distance between the
#         library DNA sequence from the match and a given library sequence is less than the minimum distance between all
#         libraries in `included_libs`
#
#         When mode is "normal" the library will be called successfully if the smallest Levenshtein distance between the
#         library DNA sequence from the match and a given library is less than `LIBRARY_RISK_NORMAL_THRESHOLD` compared to
#         the next smallest distance to another library
#
#         When mode is "low-quality" the library will be called based on the library sequence with the lowest Levenshtein
#         distance to the matched library DNA sequence
#
#         `called_lib` will be "FAILED" if a library cannot be assigned using the respective method
#
#         Parameters
#         ----------
#         included_libs: list[str], default = None
#             the list of library names that should be searched for
#             library names should exist with a corresponding sequence in the `/data/libs.json` file.
#         min_dist: int, default = None
#             the min Levenshtein distance between the passed libraries
#
#         Notes
#         -----
#         "low-quality" mode can never fail to call a library
#         if left as `None`, included_libs will default to all libs in the config files
#         """
#         distances = get_dist_from_lib(self.lib_barcode, included_libs)
#         self.called_lib = "FAILED"
#         self.called_lib_risk = 100
#
#         if self.mode == "high-quality":
#             if min_dist is None:
#                 _dist_thresh = get_min_index_distance(included_libs) / 2
#             else:
#                 _dist_thresh = min_dist / 2
#             for _lib_key, dist in distances.items():
#                 if dist < _dist_thresh:
#                     self.called_lib = _lib_key
#                     self.called_lib_risk = dist
#                     break
#                 if dist < self.called_lib_risk:
#                     self.called_lib_risk = dist
#         elif self.mode == "normal":
#             _curr_dist = 100
#             for _lib_key, dist in distances.items():
#                 if dist < _curr_dist:
#                     if (_curr_dist - dist) < LIBRARY_RISK_NORMAL_THRESHOLD:
#                         self.called_lib = "FAILED"
#                     else:
#                         self.called_lib = _lib_key
#                     _curr_dist = dist
#             self.called_lib_risk = _curr_dist
#         else:
#             _curr_dist = 100
#             _curr_fail = False
#             for _lib_key, dist in distances.items():
#                 if dist < _curr_dist:
#                     self.called_lib = _lib_key
#                     _curr_dist = dist
#             self.called_lib_risk = _curr_dist
#
#     def to_row(self) -> str:
#         """
#         Convert the CalledSequence to a string representation that can be written to a csv file
#
#         Returns
#         -------
#         csv_row: str
#             the CalledSequence as a comma seperated string and newline terminal chr
#         """
#         # REMOVE ONCE WE SWAP TO HAMMING ENCODED LIBRARIES
#         hamming_fixable = all([
#             self.bba_call not in ["failed_insert", "failed_delete", "FAILED"],
#             self.bbb_call not in ["failed_insert", "failed_delete", "FAILED"],
#             self.bbc_call not in ["failed_insert", "failed_delete", "FAILED"],
#             self.del_id == "FAILED",
#             self.called_index != "FAILED"
#         ])
#
#         return (
#             f'{self.sequence_match.match_id},{self.sequence_match.read_id},'
#             f'{self.del_id != "FAILED" and self.called_index != "FAILED" and self.called_lib != "FAILED"},'
#             f'{hamming_fixable},'
#             f'{self.called_lib_risk+self.called_index_risk+self.bba_toe_risk+self.bbb_toe_risk+self.bbc_toe_risk},'
#             f'{self.match_risk},{self.called_lib},{self.called_lib_risk},{self.called_index},{self.called_index_risk},'
#             f'{self.del_id},{self.bba_call},{self.bba_toe_risk},{self.bbb_call},{self.bbb_toe_risk},{self.bbc_call},'
#             f'{self.bbc_toe_risk},{self.umi},{self.seq}\n'
#         )
#
#     def __str__(self):
#         return self.del_id
#
#     def __repr__(self):
#         return str(self)
#
#
# def call_hit(hit: SequenceMatch, libs: list[str], indexes: list[str], make: DELMake, strict: bool = False,
#              min_lib_dist: int = None, min_idx_dist: int = None) -> CalledSequence or None:
#     """
#     A helper function the get the multiprocessing module to run without failing in `sequence_caller.py`
#
#     Parameters
#     ----------
#     hit: SequenceMatch
#         the sequence match to call
#     libs: list[str]
#         the list of libraries to search for when calling the library
#     indexes: list[str]
#         the list of indexes to search for when calling the index
#     make: DELMake
#         the DELMake object defining the make of the barcode for the SequenceMatch
#     strict: boolean, default = False
#         whether to use the strict alignment mode
#     min_lib_dist: int, default = None
#         the minimum distance between any pair of library sequences
#         if left as `None` will generate this number on the fly using `libs`
#     min_idx_dist: int, default = None
#         the minimum distance between any pair of index sequences
#         if left as `None` will generate this number on the fly using `indexes`
#
#     Notes
#     -----
#     it is expensive to do the on-the-fly min_dist calculations.
#     For the most part, nearly every call during a DELi run will have the same set of indexes and libraries.
#     So it is advised to precalculate these and pass them to avoid this slow down
#
#     Returns
#     -------
#     called_hit: CalledSequence or None
#         the CalledSequence object or None if the calling sequence raised an error
#         will log the error in the log file
#     """
#
#     if min_lib_dist is None:
#         min_lib_dist = get_min_lib_distance(libs)
#     if min_idx_dist is None:
#         min_idx_dist = get_min_index_distance(indexes)
#
#     try:
#         _called_hit = CalledSequence(hit, mode="high-quality", strict=strict, make=make)
#         _called_hit.call_lib(included_libs=libs, min_dist=min_lib_dist)
#         _called_hit.call_index(included_index=indexes, min_dist=min_idx_dist)
#         _called_hit.call_bbs()
#         return _called_hit
#     except Exception:
#         return None
