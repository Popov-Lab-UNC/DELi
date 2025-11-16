"""code for demultiplexing in DELi"""
import abc
import itertools
import sys
import warnings
from collections import defaultdict
from typing import Optional, Iterator, no_type_check, TypeVar, Generic

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from cutadapt.adapters import SingleMatch, SingleAdapter, AnywhereAdapter, FrontAdapter, \
    RightmostBackAdapter, LinkedAdapter, LinkedMatch
from cutadapt.adapters import Match as OGCutadaptMatch
from dnaio import SequenceRecord
from numba import njit
from regex import Pattern, Match, compile, BESTMATCH

from deli.dels.library import DELibraryCollection, DELibrary
from deli.dels.barcode import BarcodeSchema

from .barcode_calling import ValidCall, BarcodeCaller, get_barcode_caller
from ._base import FailedDecodeAttempt


# good old type hinting stuff
Q_co = TypeVar("Q_co", bound="_Query", covariant=True)
M_co = TypeVar("M_co", bound="_Match", covariant=True)
Q = TypeVar("Q", bound="_Query")


class AlignedSeq:
    """A sequence that has been aligned to a reference"""

    def __init__(self, sequence: SequenceRecord, section_spans: dict[str, tuple[int, int]]) -> None:
        self.sequence = sequence
        self.section_spans = section_spans

    def get_span(self, section_name: str) -> str:
        """Get the library sequence from the aligned sequence"""
        start, stop = self.section_spans[section_name]
        return self.sequence.sequence[start:stop]


class SectionSequenceAligner:
    def __init__(self, barcode_schema: BarcodeSchema, alignment_sections: tuple[str, ...]):
        self.barcode_schema = barcode_schema
        self.alignment_sections = alignment_sections

        _required_sections = self.barcode_schema.get_required_sections()
        self._spans = self.barcode_schema.get_section_spans(exclude_overhangs=False)
        self._span_map = {sec_nam: (0, 100000) for sec_nam in _required_sections}

        self._adjustment_table: dict[str, dict[str, tuple[int, int]]] = dict()
        for alignment_section in self.alignment_sections:
            # alignment_section_span = self._spans[alignment_section]
            self._adjustment_table[alignment_section] = dict()
            for required_section in _required_sections:
                required_barcode_section = self.barcode_schema.get_section(required_section)
                dist = self.barcode_schema.get_length_between_sections(alignment_section, required_section, include_direction=True)
                if dist < 0:
                    adjusted_stop = dist - 1 - (len(required_barcode_section.section_overhang) if required_barcode_section.section_overhang else 0)
                    adjusted_start = adjusted_stop - len(required_barcode_section.section_tag)
                else:
                    adjusted_start = dist + 1
                    adjusted_stop = adjusted_start + len(required_barcode_section.section_tag)
                self._adjustment_table[alignment_section][required_section] = (adjusted_start, adjusted_stop)

    def align_sequence(self, sequence: SequenceRecord, alignment_section_spans: dict[str, tuple[int, int]]) -> AlignedSeq:
        _span_map = self._span_map.copy()
        for alignment_section, (align_start, align_stop) in alignment_section_spans.items():
            for required_section, (adj_start, adj_stop) in self._adjustment_table[alignment_section].items():
                cur_start, cur_stop = _span_map[required_section]
                _span_map[required_section] = (min(align_start + adj_start, 0, cur_start), max(align_stop + adj_stop, len(sequence), cur_stop))

        return AlignedSeq(sequence=sequence, section_spans=_span_map)

    def get_library_seq(self, sequence: SequenceRecord, alignment_section_spans: dict[str, tuple[int, int]]) -> str:
        _library_span = (0, 100000)
        for alignment_section, (align_start, align_stop) in alignment_section_spans.items():
            adj_start, adj_stop = self._adjustment_table[alignment_section]["library"]
            _library_span = (min(align_start + adj_start, 0, _library_span[0]), max(align_stop + adj_stop, len(sequence), _library_span[1]))
        return sequence.sequence[_library_span[0]:_library_span[1]]

    def iter_possible_library_seqs(self, sequence: SequenceRecord, alignment_section_spans: dict[str, tuple[int, int]]) -> Iterator[str]:
        library_span= self.get_library_seq(sequence, alignment_section_spans)
        library_tag_length = len(self.barcode_schema.library_section.section_tag)
        for length in [library_tag_length, library_tag_length - 1]:
            start = 0
            end = start + length
            while end <= len(library_span):
                yield library_span[start: end]
                start += 1
                end = start + length


class LibraryAligner:
    def __init__(self, barcode_schema: BarcodeSchema, include_library_section: bool = True):
        self.barcode_schema = barcode_schema

        # Custom substitution matrix with 0 score for 'N' mismatches
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
        self.barcode_ref = barcode_schema.get_full_barcode()
        self.spans = barcode_schema.get_section_spans(exclude_overhangs=True)
        self._required_sections = barcode_schema.get_required_sections()
        if not include_library_section:
            self._required_sections.remove("library")

    def align(self, sequence: SequenceRecord) -> Iterator[tuple[AlignedSeq, float]]:
        aligned_sec_spans: dict[str, tuple[int, int]] = dict()
        alignments = self.aligner.align(sequence.sequence, self.barcode_ref)
        for alignment in alignments:
            _reverse_idx = _query_to_ref_map(alignment.coordinates)
            for barcode_sec in self._required_sections:
                sec_span = self.spans[barcode_sec]
                aligned_sec_spans[barcode_sec] = (int(_reverse_idx[sec_span.start]), int(_reverse_idx[sec_span.stop - 1] + 1))
            yield AlignedSeq(sequence, aligned_sec_spans), alignment.score


class FailedStaticAlignment(FailedDecodeAttempt):
    """A sequence that failed to align to a static reference"""
    def __init__(self, sequence: SequenceRecord):
        super().__init__(sequence, "Failed to locate static observed_barcode regions")


class FailedLibraryBarcodeLookup(FailedDecodeAttempt):
    """A sequence that failed to align to a static reference"""
    def __init__(self, sequence: SequenceRecord):
        super().__init__(sequence, "Failed to locate static observed_barcode regions")


class _Query(Generic[M_co], abc.ABC):
    def __init__(self, libraries: list[DELibrary], section_names: tuple[str, ...]):
        self.libraries = libraries
        self.section_names = section_names

        self._sec_name2seq: dict[str, str] = dict()
        self._dist_between_section: dict[tuple[str, str], list[int]] = dict()

        for library in self.libraries:
            _last_idx = -1
            _prev_sec_name: str | None = None
            for section_name in section_names:
                try:
                    sec = library.barcode_schema.get_section(section_name)
                except Exception as e:
                    raise LibraryDemultiplexerError(
                        f"Barcode section name {section_name} not found in query "
                        f"covered library {library.library_id}"
                    ) from e
                if section_name not in self._sec_name2seq:
                    self._sec_name2seq[sec.section_name] = sec.get_dna_sequence()
                elif self._sec_name2seq[sec.section_name] != sec.get_dna_sequence():
                    raise LibraryDemultiplexerError(
                        f"Barcode section {section_name} for library {library.library_id} does not match "
                        f"sequences across libraries covered by the query; "
                        f"expected {self._sec_name2seq[sec.section_name]}, found {sec.get_dna_sequence()}"
                    )
                cur_idx = library.barcode_schema.barcode_sections.index(sec)
                if cur_idx <= _last_idx:
                    raise LibraryDemultiplexerError(
                        f"Barcode section {section_name} for library {library.library_id} is out of order; "
                        f"expected order {section_names}"
                    )
                else:
                    _last_idx = cur_idx

                # this will collect the min and max distances between the sections
                if _prev_sec_name is not None:
                    dist = library.barcode_schema.get_length_between_sections(_prev_sec_name, section_name)
                    key = (_prev_sec_name, section_name)
                    if key not in self._dist_between_section:
                        self._dist_between_section[key] = [dist, dist]
                    elif dist > self._dist_between_section[key][1]:
                        self._dist_between_section[key][1] = dist
                    elif dist < self._dist_between_section[key][0]:
                        self._dist_between_section[key][0] = dist
                _prev_sec_name = section_name

        # get trim adjustments to required regions for each library
        self.trim_adjustments: dict[DELibrary, tuple[int, int]] = dict()
        for library in self.libraries:
            _min_sec_start = 100000
            _max_sec_end = 0
            for section_name in section_names:
                static_section_span = library.barcode_schema.section_spans[section_name]
                _min_sec_start = min(_min_sec_start, static_section_span.start)
                _max_sec_end = max(_max_sec_end, static_section_span.stop)
            start_adj = min(library.barcode_schema.required_start - _min_sec_start, 0)
            end_adj = max(library.barcode_schema.required_end - _max_sec_end, 0)
            self.trim_adjustments[library] = (start_adj, end_adj)

    @abc.abstractmethod
    def __hash__(self):
        raise NotImplementedError()

    @abc.abstractmethod
    def search(self, sequence: str) -> Optional[M_co]:
        raise NotImplementedError()


class _Match(Generic[Q_co], abc.ABC):
    def __init__(self, query: Q_co, score: float | int):
        self.query = query
        self.score = score

    @abc.abstractmethod
    def get_section_spans(self) -> dict[str, tuple[int, int]]:
        raise NotImplementedError()

    @abc.abstractmethod
    def trim_seq(self, library_call: ValidCall[DELibrary], sequence: SequenceRecord) -> SequenceRecord:
        raise NotImplementedError()


class RegexError(Exception):
    """Exception raised for errors in the regex demultiplexing process."""
    pass


class _RegexQuery(_Query["_RegexMatch"]):
    def __init__(self, libraries: list[DELibrary], section_names: tuple[str, ...], error_tolerance: int = 1):
        super().__init__(libraries, section_names)

        # build the regex string
        self.pattern = ""
        for i ,(section_name, sequence) in enumerate(self._sec_name2seq.items()):
            self.pattern += f"(?:(?P<{section_name}>{sequence})){{e<={error_tolerance}}}"
            if i < len(self._sec_name2seq) - 1:  # flanked by other sections
                key = (section_names[i], section_names[i + 1])
                min_dist, max_dist = self._dist_between_section[key]
                self.pattern += f".{{{min_dist-5},{max_dist+5}}}"
        self.compiled_pattern: Pattern = compile(self.pattern, BESTMATCH)

    def __hash__(self):
        return hash(self.pattern)

    def search(self, sequence: str) -> Optional["_RegexMatch"]:
        """Search for the pattern in the sequence"""
        match = self.compiled_pattern.search(sequence)
        return _RegexMatch(self, sum(match.fuzzy_counts), match) if match is not None else None


class _RegexMatch(_Match[_RegexQuery]):
    def __init__(self, query: _RegexQuery, score: float | int, match: Match):
        super().__init__(query, score)
        self.match = match

    def get_section_spans(self) -> dict[str, tuple[int, int]]:
        """Get the spans of the sections from a regex match"""
        return {
            section_name: span for section_name, span in zip(
                self.query.section_names, self.match.regs[1:], strict=True
            )
        }

    def trim_seq(self, library_call: ValidCall[DELibrary], sequence: SequenceRecord) -> SequenceRecord:
        """Trim a sequence based given the regex pattern used"""
        start_adj, end_adj = self.query.trim_adjustments[library_call.obj]
        # build the trim
        trim = slice(
            max(0, min(self.match.regs[1][0] - end_adj, library_call.obj.barcode_schema.required_start) - 10),
            max(self.match.regs[-1][1] + end_adj, library_call.obj.barcode_schema.required_end) + 10
        )
        return sequence[trim]


class CutAdaptError(Exception):
    """Exception raised for errors in the cutadapt demultiplexing process."""
    pass


class _CutadaptQuery(_Query["_CutadaptMatch"]):
    def __init__(self, libraries: list[DELibrary], section_names: tuple[str, ...], error_tolerance: int = 1, min_overlap: int = 8):
        if len(section_names) > 2:
            raise CutAdaptError(
                "Cutadapt based queries can only handle up to 2 static sections; "
                "found {len(section_names)} sections"
            )
        super().__init__(libraries, section_names)

        self._hash_name = "---".join(self._sec_name2seq.values())

        # determine which type of adapter to build for each section
        adapters: list[SingleAdapter] = []
        for section in section_names:
            distance_from_front_percent: float = 0.0
            distance_from_back_percent: float = 0.0
            for library in self.libraries:
                start, end = library.barcode_schema.section_spans[section]
                distance_from_front_percent = max(distance_from_front_percent, start/len(library.barcode_schema))
                distance_from_back_percent = max(distance_from_back_percent, (len(library.barcode_schema) - end)/len(library.barcode_schema))

            if ((distance_from_front_percent > 0.1) and (distance_from_back_percent > 0.1)) or ((distance_from_front_percent <= 0.1) and (distance_from_back_percent <= 0.1)) :
                warnings.warn(
                    "Cutadapt based queries may not perform well when a given static section "
                    "is located in the middle or on both sides of the barcode for libraries in the query"
                )
                adapters.append(AnywhereAdapter(
                    sequence=self._sec_name2seq[section],
                    max_errors=error_tolerance,
                    min_overlap=min_overlap,
                    name=section,
                ))
            elif distance_from_front_percent <= 0.1 :
                adapters.append(FrontAdapter(
                    seq=self._sec_name2seq[section],
                    max_errors=error_tolerance,
                    min_overlap=min_overlap,
                    name=section,
                ))
            else:
                adapters.append(RightmostBackAdapter(
                    seq=self._sec_name2seq[section],
                    max_errors=error_tolerance,
                    min_overlap=min_overlap,
                    name=section,
                ))

        # if two adapters, link them together
        self.adapter: LinkedAdapter | SingleAdapter
        if len(adapters) == 2:
            self.adapter = LinkedAdapter(
                front_adapter=adapters[0],
                back_adapter=adapters[1],
                front_required=True,
                back_required=True,
                name=None
            )
        else:
            self.adapter = adapters[0]

    def __hash__(self):
        return hash(self._hash_name)

    def search(self, sequence: str) -> Optional["_CutadaptMatch"]:
        """Search for the pattern in the sequence using cutadapt"""
        match: LinkedMatch | SingleMatch | None = self.adapter.match_to(sequence)
        if match is None:
            return None
        elif isinstance(match, LinkedMatch):
            return _LinkedCutadaptMatch(self, match.errors, match)
        else:
            return _SingleCutadaptMatch(self, match.errors, match)


class _CutadaptMatch(_Match[_CutadaptQuery], abc.ABC):
    match: OGCutadaptMatch
    def __init__(self, query: _CutadaptQuery, score: float | int):
        super().__init__(query, score)


class _LinkedCutadaptMatch(_CutadaptMatch):
    def __init__(self, query: _CutadaptQuery, score: float | int, match: LinkedMatch):
        super().__init__(query, score)
        self.match: LinkedMatch = match

    def get_section_spans(self) -> dict[str, tuple[int, int]]:
        """Get the spans of the sections from a cutadapt match"""
        front_match = self.match.front_match
        back_match = self.match.back_match

        return {
            front_match.adapter.name: (front_match.rstart, front_match.rstop),
            back_match.adapter.name: (back_match.rstart, back_match.rstop),
        }

    def trim_seq(self, library_call: ValidCall[DELibrary], sequence: SequenceRecord) -> SequenceRecord:
        """Trim a sequence given the cutadapt match used"""
        start_adj, end_adj = self.query.trim_adjustments[library_call.obj]
        # build the trim
        trim_start = max(0, min(self.match.front_match.rstart - end_adj, library_call.obj.barcode_schema.required_start) - 10)
        trim_end = max(self.match.back_match.rstop + end_adj, library_call.obj.barcode_schema.required_end) + 10
        return sequence[trim_start: trim_end]


class _SingleCutadaptMatch(_CutadaptMatch):
    def __init__(self, query: _CutadaptQuery, score: float | int, match: SingleMatch):
        super().__init__(query, score)
        self.match: SingleMatch = match

    def get_section_spans(self) -> dict[str, tuple[int, int]]:
        """Get the spans of the sections from a cutadapt match"""
        return {
            self.match.adapter.name: (self.match.rstart, self.match.rstop),
        }

    def trim_seq(self, library_call: ValidCall[DELibrary], sequence: SequenceRecord) -> SequenceRecord:
        """Trim a sequence given the cutadapt match used"""
        start_adj, end_adj = self.query.trim_adjustments[library_call.obj]
        trim_start = max(0, min(self.match.rstart - end_adj, library_call.obj.barcode_schema.required_start) - 10)
        trim_end = max(self.match.rstop + end_adj, library_call.obj.barcode_schema.required_end) + 10
        return sequence[trim_start: trim_end]


class LibraryDemultiplexerError(Exception):
    """Exception raised for errors in the library demultiplexing process."""
    pass


class LibraryDemultiplexer(abc.ABC):
    def __init__(self, libraries: DELibraryCollection, revcomp: bool = True):
        self.revcomp = revcomp
        self.libraries = libraries

    @abc.abstractmethod
    def demultiplex(self, sequence: SequenceRecord) -> tuple[ValidCall[DELibrary], Iterator[tuple[AlignedSeq, float]]] | FailedDecodeAttempt:
        raise NotImplementedError()


class QueryBasesLibraryDemultiplexer(LibraryDemultiplexer, Generic[Q], abc.ABC):
    def __init__(self, *args, realign: bool = False, build_dynamic_aligners_: bool = False, build_section_seq_aligners_: bool = False, **kwargs):
        self.realign = realign
        super().__init__(*args, **kwargs)

        self._library_aligners: dict[DELibrary, LibraryAligner] = dict()
        self._section_seq_aligners: dict[DELibrary, SectionSequenceAligner] = dict()
        self._unique_section_seq_aligners: dict[Q, list[SectionSequenceAligner]] = dict()
        self._queries: list[Q] = self._get_queries()

        if build_dynamic_aligners_ or self.realign:  # only initialize these if we need them
            self._library_aligners = self._get_library_aligners()

        if build_section_seq_aligners_:  # only initialize these if we need them
            self._section_seq_aligners = self._get_section_seq_aligners()
            for query in self._queries:
                aligner_group: list[SectionSequenceAligner] = []
                aligner_group_ids: set[int] = set()
                for library in query.libraries:
                    aligner = self._section_seq_aligners[library]
                    if id(aligner) not in aligner_group_ids:
                        aligner_group.append(aligner)
                self._unique_section_seq_aligners[query] = aligner_group

    @abc.abstractmethod
    def _get_queries(self) -> list[Q]:
        raise NotImplementedError()

    def _get_library_aligners(self) -> dict[DELibrary, LibraryAligner]:
        aligners: dict[DELibrary, LibraryAligner] = {}
        for library in self.libraries:
            _found_existing: bool = False
            for other_library, aligner in aligners.items():
                if library.barcode_schema == other_library.barcode_schema:
                    aligners[library] = aligner
                    _found_existing = True
                    break
            if not _found_existing:
                aligners[library] = LibraryAligner(library.barcode_schema)
        return aligners

    def _get_section_seq_aligners(self) -> dict[DELibrary, SectionSequenceAligner]:
        aligners: dict[DELibrary, SectionSequenceAligner] = {}
        for query in self._queries:
            for library in query.libraries:
                _found_existing: bool = False
                for other_library, aligner in aligners.items():
                    if library.barcode_schema == other_library.barcode_schema:
                        aligners[library] = aligner
                        _found_existing = True
                        break
                if not _found_existing:
                    aligners[library] = SectionSequenceAligner(library.barcode_schema, query.section_names)
        return aligners

    def _get_best_match(self, sequence: SequenceRecord) -> _Match | None:
        best_match: _Match | None = None
        best_score: float = float('inf')

        for query in self._queries:
            match = query.search(sequence.sequence)
            if match is None:
                continue
            elif match.score < best_score:
                best_match = match
                best_score = match.score
                if best_score == 0:
                    break  # perfect match found

        if self.revcomp and (best_score > 0):
            best_revcomp_match: _Match | None = None
            best_revcomp_score: float = float('inf')

            # try the reverse complement
            revcomp_seq = sequence.reverse_complement()
            for query in self._queries:
                match = query.search(revcomp_seq.sequence)
                if match is None:
                    continue
                elif match.score < best_score:
                    best_revcomp_match = match
                    best_revcomp_score = match.score
                    if best_revcomp_score == 0:
                        break  # perfect match found

            if best_revcomp_score < best_score:
                best_match = best_revcomp_match
                # update sequence record in place to be revcomp
                sequence.sequence = revcomp_seq.sequence
                sequence.name += " [revcomp]"
                if sequence.qualities is not None:
                    sequence.qualities = revcomp_seq.qualities

        return best_match

    def _make_alignments(self, best_match: _Match, library_call: ValidCall[DELibrary], sequence: SequenceRecord, _static_alignment_spans: Optional[dict[str, tuple[int, int]]] = None)-> Iterator[tuple[AlignedSeq, float]]:
        if self.realign:
            # this will trim the sequence to the region we care about with a 10bp buffer on each side
            trimmed_seq = best_match.trim_seq(library_call, sequence)
            realigner: LibraryAligner = self._library_aligners[library_call.obj]
            return realigner.align(trimmed_seq)
        else:
            _static_alignment_spans = best_match.get_section_spans() if _static_alignment_spans is None else _static_alignment_spans
            return iter([(self._section_seq_aligners[library_call.obj].align_sequence(
                sequence, _static_alignment_spans
            ), 0.0)])

class LibraryTagLibraryDemultiplexer(QueryBasesLibraryDemultiplexer[Q], abc.ABC):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, build_section_seq_aligners_=True, **kwargs)

    def demultiplex(self, sequence: SequenceRecord) -> tuple[ValidCall[DELibrary], Iterator[tuple[AlignedSeq, float]]] | FailedDecodeAttempt:
        # query, score, regex.Match object
        best_match = self._get_best_match(sequence)

        # failed to match any static region query
        if best_match is None:
            return FailedStaticAlignment(sequence)
        else:
            library_call = ValidCall(best_match.query.libraries[0], best_match.score)
            return library_call, self._make_alignments(best_match, library_call, sequence)


class LibraryTagRegexLibraryDemultiplexer(LibraryTagLibraryDemultiplexer[_RegexQuery]):
    def __init__(self, libraries: DELibraryCollection, error_tolerance: int = 1, realign: bool = True, revcomp: bool = True):
        self.error_tolerance = error_tolerance
        super().__init__(libraries=libraries, realign=realign, revcomp=revcomp)

    def _get_queries(self) -> list[_RegexQuery]:
        queries: list[_RegexQuery] = list()
        for library in self.libraries:
            before_lib_section_names = [
                sec.section_name for sec in library.barcode_schema.get_static_sections_before_library()
            ]
            after_lib_section_names = [
                sec.section_name for sec in library.barcode_schema.get_static_sections_after_library()
            ]
            queries.append(
                _RegexQuery(
                    libraries=[library],
                    section_names=tuple(before_lib_section_names + ["library"] + after_lib_section_names),
                    error_tolerance=self.error_tolerance
                )
            )
        return queries


class LibraryTagCutadaptLibraryDemultiplexer(LibraryTagLibraryDemultiplexer[_CutadaptQuery]):
    def __init__(self, libraries: DELibraryCollection, error_tolerance: int = 1, realign: bool = True, min_overlap: int = 8, revcomp: bool = True):
        self.error_tolerance = error_tolerance
        self.min_overlap = min_overlap
        super().__init__(libraries=libraries, realign=realign, revcomp=revcomp)

    def _get_queries(self) -> list[_CutadaptQuery]:
        queries: list[_CutadaptQuery] = list()
        for library in self.libraries:
            queries.append(
                _CutadaptQuery(
                    libraries=[library],
                    section_names=tuple(["library"]),
                    error_tolerance=self.error_tolerance
                )
            )
        return queries


class NoLibraryTagLibraryDemultiplexer(QueryBasesLibraryDemultiplexer[Q], abc.ABC):
    def __init__(self, *args, error_correction_mode_str: str = "levenshtein_dist:1,asymmetrical", **kwargs):
        super().__init__(*args, build_section_seq_aligners_=True, **kwargs)

        # maps queries to barcode callers for the libraries in that query
        self._library_callers: dict[_RegexQuery, BarcodeCaller[DELibrary]] = {
            query: get_barcode_caller(
                {library.library_tag: library for library in query.libraries},
                error_correction_mode_str=error_correction_mode_str
            ) for query in self._queries
        }

    @abc.abstractmethod
    def _get_query_object(self, libraries: list[DELibrary], section_names: tuple[str, ...]) -> type[Q]:
        raise NotImplementedError()

    def demultiplex(self, sequence: SequenceRecord) -> tuple[ValidCall[DELibrary], Iterator[tuple[AlignedSeq, float]]] | FailedDecodeAttempt:
        best_match = self._get_best_match(sequence)

        # failed to match any static region query
        if best_match is None:
            return FailedStaticAlignment(sequence)

        # these are the spans of the static regions from the match
        _static_alignment_spans = best_match.get_section_spans()

        # Now find the best matching library within the best regex query
        library_call: ValidCall[DELibrary] | None = None
        if len(best_match.query.libraries) == 1:
            library_call = ValidCall(best_match.query.libraries[0], 0.0) # only one possible library
        else:
            _best_call_score = 100  # lower scores are better

            # This could probably be optimized by avoiding re-extraction of the same sequences
            # if we know that the static sections are the same distance from the library tags
            # across several aligners. Instead of doing that now, we just track which sequences
            # we've already tried to avoid some redundant work.
            observed_possible_lib_seqs: set[str] = set()

            for init_aligner in self._unique_section_seq_aligners[best_match.query]:
                _best_internal_call = None
                _best_internal_score = 100
                for possible_library_seq in init_aligner.iter_possible_library_seqs(
                    sequence, _static_alignment_spans
                ):
                    # track which sequences we've already tried
                    if possible_library_seq in observed_possible_lib_seqs:
                        continue  # skip if we've already tried this sequence
                    else:
                        observed_possible_lib_seqs.add(possible_library_seq)

                    _library_call = self._library_callers[best_match.query].decode_barcode(possible_library_seq)
                    if _library_call is not None:
                        if _library_call.score < _best_internal_score:
                            _best_internal_call = _library_call
                            _best_call_score = _library_call.score

                if _best_internal_call is not None:
                    if _best_call_score < _best_call_score:
                        library_call = _best_internal_call
                        _best_call_score = _best_call_score

        # at this point we did everything we can to try and call the library
        if library_call is None:
            return FailedLibraryBarcodeLookup(sequence)
        else:
            return library_call, self._make_alignments(best_match, library_call, sequence, _static_alignment_spans=_static_alignment_spans)


class SinglePrimerLibraryDemultiplexer(NoLibraryTagLibraryDemultiplexer[Q], abc.ABC):
    """
    Demultiplexer for single primer regex style demultiplexing
    """
    def __init__(self, *args, error_tolerance: int = 1, **kwargs):
        """
        Initialize the SinglePrimerDemultiplexer

        Parameters
        ----------
        libraries: DELibraryCollection
            the libraries to use for demultiplexing
        realign: bool
            whether to realign the sequences after demultiplexing
            If also will use the regex alignment
            This can be prone to higher error rates
        """
        self.error_tolerance = error_tolerance
        super().__init__(*args, **kwargs)

    def _get_queries(self) -> list[Q]:
        """generate the regex query sets for the single primer demultiplexer"""
        """
        Given a list of DELs, get the set of demultiplex queries required

        Returns a mapping of static sequence to the list of DELs it will
        be able to demultiplex.

        This algorithm greedily picks sequences that cover the most remaining
        libraries until all libraries are covered. It does not attempt to optimize
        for distance between chosen sequences or length of sequences. In practice,
        this should still yield good results since most DEL runs will have use
        libraries that have only 6-8 unique static regions among them. If its
        becomes problematic, a more complex algorithm can be implemented later.
        """

        # map all possible static observed_barcode to the observed_barcode schemas they cover
        static_seq_groups: defaultdict[tuple[str, str], set[int]] = defaultdict(set)
        for idx, library in enumerate(self.libraries):
            static_section_tags = [(sec.get_dna_sequence(), sec.section_name) for sec in library.barcode_schema.static_sections]
            for sec_seq, sec_name in static_section_tags:
                static_seq_groups[(sec_seq, sec_name)].add(idx)

        # greedily pick sequences that cover the most remaining observed_barcode schemas
        candidates = static_seq_groups.copy()
        uncovered = set(range(len(self.libraries)))
        chosen: list[Q] = []
        while uncovered:
            best_seq: tuple[str, str] | None = None
            best_cover: int = 0
            for _sec_seq, grp in candidates.items():
                cover = len(grp & uncovered)
                if cover > best_cover:
                    best_cover = cover
                    best_seq = _sec_seq
            if best_seq is None:
                raise RuntimeError(
                    f"this is impossible, something went wrong. Please raise an issue"
                )
            else:
                new_idxes = set(candidates[best_seq]) & uncovered
                libraries: list[DELibrary] = [self.libraries.libraries[i] for i in new_idxes]

                chosen.append(self._get_query_object(
                    libraries=libraries,
                    section_names=(best_seq[1],)
                ))

                uncovered -= static_seq_groups[best_seq]
                del candidates[best_seq]
        return chosen


class FlankingPrimersLibraryDemultiplexer(NoLibraryTagLibraryDemultiplexer[Q], abc.ABC):
    """
    Demultiplexer for single primer regex style demultiplexing
    """
    def __init__(self, *args, error_tolerance: int = 1, **kwargs):
        """
        Initialize the SinglePrimerDemultiplexer

        Parameters
        ----------
        libraries: DELibraryCollection
            the libraries to use for demultiplexing
        realign: bool
            whether to realign the sequences after demultiplexing
            If also will use the regex alignment
            This can be prone to higher error rates
        """
        self.error_tolerance = error_tolerance
        super().__init__(*args, **kwargs)

    def _get_queries(self) -> list[Q]:
        """generate the regex query sets for the flanking primer demultiplexer"""

        # map all possible static observed_barcode to the observed_barcode schemas they cover
        static_seq_groups: defaultdict[tuple[tuple[str, str], tuple[str, str]], set[int]] = defaultdict(set)
        for idx, library in enumerate(self.libraries):
            before_lib_sections = library.barcode_schema.get_static_sections_before_library()
            after_lib_sections = library.barcode_schema.get_static_sections_after_library()

            if len(before_lib_sections) == 0 or len(after_lib_sections) == 0:
                raise LibraryDemultiplexerError(
                    f"Library {library.library_id} does not have static sections both "
                    f"before and after the library barcode section; "
                    f"cannot use a FlankingPrimersRegexLibraryDemultiplexer"
                )

            before_seqs_and_names = [(sec.section_name, sec.get_dna_sequence()) for sec in before_lib_sections]
            after_seqs_and_names = [(sec.section_name, sec.get_dna_sequence()) for sec in after_lib_sections]

            possible_combos = itertools.product(before_seqs_and_names, after_seqs_and_names)

            for combo in possible_combos:
                static_seq_groups[combo].add(idx)

        # greedily pick sequences that cover the most remaining observed_barcode schemas
        candidates = static_seq_groups.copy()
        uncovered = set(range(len(self.libraries)))
        chosen: list[_RegexQuery] = []
        while uncovered:
            best_seq: tuple[tuple[str, str], tuple[str, str]] | None = None
            best_cover: int = 0
            for _sec_seq, grp in candidates.items():
                cover = len(grp & uncovered)
                if cover > best_cover:
                    best_cover = cover
                    best_seq = _sec_seq
            if best_seq is None:
                raise RuntimeError(
                    f"this is impossible, something went wrong. Please raise an issue"
                )
            else:
                new_idxes = set(candidates[best_seq]) & uncovered
                libraries: list[DELibrary] = [self.libraries.libraries[i] for i in new_idxes]

                chosen.append(self._get_query_object(
                    libraries=libraries,
                    section_names=tuple(list(sec_name for sec_name, _ in best_seq)),
                ))

                uncovered -= static_seq_groups[best_seq]
                del candidates[best_seq]
        return chosen


class SinglePrimerCutadaptLibraryDemultiplexer(SinglePrimerLibraryDemultiplexer[_CutadaptQuery]):
    """
    Demultiplexer for single primer regex style demultiplexing
    """
    def __init__(self, libraries: DELibraryCollection, revcomp: bool = True, realign: bool = False, error_tolerance: int = 1, error_correction_mode_str: str = "levenshtein_dist:1,asymmetrical", min_overlap: int = 8):
        """
        Initialize the SinglePrimerDemultiplexer

        Parameters
        ----------
        libraries: DELibraryCollection
            the libraries to use for demultiplexing
        realign: bool
            whether to realign the sequences after demultiplexing
            If also will use the regex alignment
            This can be prone to higher error rates
        """
        self.min_overlap = min_overlap
        super().__init__(libraries=libraries, realign=realign, error_correction_mode_str=error_correction_mode_str, error_tolerance=error_tolerance, revcomp=revcomp)

    def _get_query_object(self, libraries: list[DELibrary], section_names: tuple[str, ...]) -> _CutadaptQuery:
        return _CutadaptQuery(
            libraries=libraries,
            section_names=section_names,
            error_tolerance=self.error_tolerance,
            min_overlap=self.min_overlap
        )


class FlankingPrimersCutadaptLibraryDemultiplexer(FlankingPrimersLibraryDemultiplexer[_CutadaptQuery]):
    """
    Demultiplexer for single primer regex style demultiplexing
    """
    def __init__(self, libraries: DELibraryCollection, revcomp: bool = True, realign: bool = False, error_tolerance: int = 1,
                 error_correction_mode_str: str = "levenshtein_dist:1,asymmetrical", min_overlap: int = 8):
        """
        Initialize the SinglePrimerDemultiplexer

        Parameters
        ----------
        libraries: DELibraryCollection
            the libraries to use for demultiplexing
        realign: bool
            whether to realign the sequences after demultiplexing
            If also will use the regex alignment
            This can be prone to higher error rates
        """
        self.min_overlap = min_overlap
        super().__init__(libraries=libraries, realign=realign, error_correction_mode_str=error_correction_mode_str, error_tolerance=error_tolerance, revcomp=revcomp)

    def _get_query_object(self, libraries: list[DELibrary], section_names: tuple[str, ...]) -> _CutadaptQuery:
        return _CutadaptQuery(
            libraries=libraries,
            section_names=section_names,
            error_tolerance=self.error_tolerance,
            min_overlap=self.min_overlap
        )


class SinglePrimerRegexLibraryDemultiplexer(SinglePrimerLibraryDemultiplexer[_RegexQuery]):
    """
    Demultiplexer for single primer regex style demultiplexing
    """
    def __init__(self, libraries: DELibraryCollection, realign: bool = False, revcomp: bool = True, error_tolerance: int = 1, error_correction_mode_str: str = "levenshtein_dist:1,asymmetrical"):
        """
        Initialize the SinglePrimerDemultiplexer

        Parameters
        ----------
        libraries: DELibraryCollection
            the libraries to use for demultiplexing
        realign: bool
            whether to realign the sequences after demultiplexing
            If also will use the regex alignment
            This can be prone to higher error rates
        """
        super().__init__(libraries=libraries, realign=realign, error_correction_mode_str=error_correction_mode_str, error_tolerance=error_tolerance, revcomp=revcomp)

    def _get_query_object(self, libraries: list[DELibrary], section_names: tuple[str, ...]) -> _RegexQuery:
        return _RegexQuery(
            libraries=libraries,
            section_names=section_names,
            error_tolerance=self.error_tolerance,
        )


class FlankingPrimersRegexLibraryDemultiplexer(FlankingPrimersLibraryDemultiplexer[_RegexQuery]):
    """
    Demultiplexer for single primer regex style demultiplexing
    """
    def __init__(self, libraries: DELibraryCollection, revcomp: bool = True, realign: bool = False, error_tolerance: int = 1,
                 error_correction_mode_str: str = "levenshtein_dist:1,asymmetrical"):
        """
        Initialize the SinglePrimerDemultiplexer

        Parameters
        ----------
        libraries: DELibraryCollection
            the libraries to use for demultiplexing
        realign: bool
            whether to realign the sequences after demultiplexing
            If also will use the regex alignment
            This can be prone to higher error rates
        """
        super().__init__(libraries=libraries, realign=realign, error_correction_mode_str=error_correction_mode_str, error_tolerance=error_tolerance, revcomp=revcomp)

    def _get_query_object(self, libraries: list[DELibrary], section_names: tuple[str, ...]) -> _RegexQuery:
        return _RegexQuery(
            libraries=libraries,
            section_names=section_names,
            error_tolerance=self.error_tolerance,
        )


class FailedFullBarcodeAlignment(FailedDecodeAttempt):
    """
    Failed to align full barcode to any library
    """
    def __init__(self, sequence: SequenceRecord):
        super().__init__(sequence, "Failed to align full barcode to any library")


class FullSeqAlignmentLibraryDemultiplexer(LibraryDemultiplexer):
    """
    Demultiplexer that aligns the full sequence to each library
    """
    def __init__(self, libraries: DELibraryCollection, revcomp: bool = True):
        """
        Initialize the FullSeqAlignmentLibraryDemultiplexer

        Notes
        -----
        For full sequence alignments there is no heuristic to avoid
        running the revcomp everytime, doubling the amount of work to
        do. It would be better to use a separate external software to
        do standardize all reads into a single orientation before
        demultiplexing with this method.

        Parameters
        ----------
        libraries: DELibraryCollection
            the libraries to use for demultiplexing
        revcomp: bool, default=True
            whether to also consider the reverse complement of the sequence
        """
        super().__init__(libraries=libraries, revcomp=revcomp)

        self.libraries = self.libraries
        self.aligner = [LibraryAligner(library.barcode_schema, include_library_section=False) for library in self.libraries]

    def demultiplex(self, sequence: SequenceRecord) -> tuple[ValidCall[DELibrary], Iterator[tuple[AlignedSeq, float]]] | FailedDecodeAttempt:
        best_alignments: tuple[ValidCall[DELibrary], Iterator[tuple[AlignedSeq, float]]] | None = None
        best_score = 0
        for library, aligner in zip(self.libraries, self.aligner):
            alignments = aligner.align(sequence)
            alignment, score = alignments.__next__()
            if score > best_score:
                best_alignments = (ValidCall(library, score=0), itertools.chain(iter([(alignment, score), ]), alignments))
                best_score = score

        # biopython alignments do not track errors; it won't be possible to know
        # if the match that exists already is perfect with processing the alignment
        # but that processing part is so expensive that it is not worth doing here
        if self.revcomp:
            revcomp_seq = sequence.reverse_complement()
            best_revcomp_alignments: tuple[ValidCall[DELibrary], Iterator[tuple[AlignedSeq, float]]] | None = None
            best_revcomp_score = 0
            for library, aligner in zip(self.libraries, self.aligner):
                alignments = aligner.align(revcomp_seq)
                alignment, score = alignments.__next__()
                if score > best_revcomp_score:
                    best_revcomp_alignments = (ValidCall(library, score=0), itertools.chain(iter([(alignment, score), ]), alignments))
                    best_revcomp_score = score
            if best_revcomp_score > best_score:
                best_alignments = best_revcomp_alignments
                # update sequence record in place to be revcomp
                sequence.sequence = revcomp_seq.sequence
                sequence.name += " [revcomp]"
                if sequence.qualities is not None:
                    sequence.qualities = revcomp_seq.qualities
        return best_alignments if best_alignments is not None else FailedFullBarcodeAlignment(sequence)


# TODO for some reason, this is needed to prevent crash when running as a CMD,
#  but not in a python interpreter... not sure why?
sys.setrecursionlimit(3000)


@no_type_check
@njit
def _query_to_ref_map(coords):
    """
    Generate map from query indices to reference indices

    Map is an array of length equal to the query sequence length
    Each index indicates the reference index that the query index maps to

    Notes
    -----
    Uses the coordinates from a Biopython alignment object
    """
    query_len = coords[1, -1]
    mapping = np.full(query_len, -1, dtype=np.int64)
    n_blocks = coords.shape[1] - 1
    for i in range(n_blocks):
        ref_start, ref_end = coords[0, i], coords[0, i + 1]
        query_start, query_end = coords[1, i], coords[1, i + 1]
        if (ref_end > ref_start) and (query_end > query_start):
            block_len = ref_end - ref_start
            for j in range(block_len):
                mapping[query_start + j] = ref_start + j
    next_val = -1
    for i in range(query_len - 1, -1, -1):
        if mapping[i] == -1:
            if next_val != -1:
                mapping[i] = next_val
        else:
            next_val = mapping[i]
    if next_val == -1:
        mapping[:] = query_len
    else:
        for i in range(query_len):
            if mapping[i] == -1:
                mapping[i] = query_len
    return mapping
