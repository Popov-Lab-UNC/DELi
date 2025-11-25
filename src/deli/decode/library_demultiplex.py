"""code for demultiplexing in DELi"""

import abc
import itertools
import sys
import warnings
from collections import defaultdict
from typing import Generic, Iterator, Optional, TypeAlias, TypeVar, no_type_check

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from cutadapt.adapters import (
    AnywhereAdapter,
    FrontAdapter,
    LinkedAdapter,
    LinkedMatch,
    RightmostBackAdapter,
    SingleAdapter,
    SingleMatch,
)
from cutadapt.adapters import Match as OGCutadaptMatch
from dnaio import SequenceRecord
from numba import njit
from regex import BESTMATCH, Match, Pattern, compile

from deli.dels.barcode import BarcodeSchema
from deli.dels.combinatorial import DELibrary, DELibraryCollection
from deli.dels.tool_compounds import TaggedToolCompoundLibrary

from ._base import FailedDecodeAttempt
from .barcode_calling import BarcodeCaller, FailedBarcodeLookup, ValidCall, get_barcode_caller


# good old type hinting stuff
Q_co = TypeVar("Q_co", bound="_Query", covariant=True)
M_co = TypeVar("M_co", bound="_Match", covariant=True)
Q = TypeVar("Q", bound="_Query")
TaggedLibrary: TypeAlias = DELibrary | TaggedToolCompoundLibrary
ValidLibraryCall: TypeAlias = ValidCall[DELibrary | TaggedToolCompoundLibrary]


class AlignedSeq:
    """
    A sequence that has been aligned to a reference

    Any algorithm could have aligned the sequence, this is just a container
    to make the observed read to that alignment.

    Notes
    -----
    Alignments in DELi are mapping of the start and end indices of a
    barcode section in the read sequence.

    Parameters
    ----------
    sequence: SequenceRecord
        the sequence that has been aligned
    section_spans: dict[str, tuple[int, int]]
        the mapping of section names to their (start, stop) indices in the sequence
        determined by the alignment algorithm
    """

    def __init__(self, sequence: SequenceRecord, section_spans: dict[str, tuple[int, int]]) -> None:
        self.sequence = sequence
        self.section_spans = section_spans

    def get_section_barcode(self, section_name: str) -> str:
        """
        Get the barcode for a given section from the aligned sequence

        Parameters
        ----------
        section_name: str
            the name of the barcode section to get the sequence for

        Returns
        -------
        str
        """
        start, stop = self.section_spans[section_name]
        return self.sequence.sequence[start:stop]


class SectionSequenceMapper(abc.ABC):
    """
    Mapper for barcodes based on known static barcode sections

    Given the known locations of pre-defined static barcode sections,
    a mapper will be able to locate and extract the part of the read
    that corresponds to other barcode sections

    It does this buy pre-compiling a map of adjustments from the static
    sections to the required sections based on the distance and direction.
    This way, at runtime it just needs to do a bit of math to locate the
    sections.

    Parameters
    ----------
    barcode_schema: BarcodeSchema
        the barcode schema to use for alignment
    alignment_sections: tuple[str, ...]
        the names of the static sections to use for alignment
    """

    def __init__(self, barcode_schema: BarcodeSchema, alignment_sections: tuple[str, ...]):
        self.barcode_schema = barcode_schema
        self.alignment_sections = alignment_sections

        self._required_sections = self.barcode_schema.get_required_section_names()
        self._spans = self.barcode_schema.get_section_spans(exclude_overhangs=False)
        self._span_map = {sec_nam: (100000, 0) for sec_nam in self._required_sections}

        self._adjustment_table: dict[str, dict[str, tuple[int, int]]] = dict()
        for alignment_section in self.alignment_sections:
            self._adjustment_table[alignment_section] = dict()
            for required_section in self._required_sections:
                required_barcode_section = self.barcode_schema.get_section(required_section)
                dist = self.barcode_schema.get_length_between_sections(
                    alignment_section, required_section, include_direction=True
                )
                if self.barcode_schema.get_direction_of_sections(alignment_section, required_section) < 0:
                    adjusted_stop = (
                        dist
                        - 1
                        - (
                            len(required_barcode_section.section_overhang)
                            if required_barcode_section.section_overhang
                            else 0
                        )
                    )
                    adjusted_start = adjusted_stop - len(required_barcode_section.section_tag)
                else:
                    adjusted_start = dist + 1
                    adjusted_stop = adjusted_start + len(required_barcode_section.section_tag)
                self._adjustment_table[alignment_section][required_section] = (adjusted_start, adjusted_stop)

    @abc.abstractmethod
    def _get_sections_to_map(self) -> list[str]:
        """Get the names of the sections to map to"""
        raise NotImplementedError()

    def calculate_section_span(
        self, section_name: str, alignment_section_spans: dict[str, tuple[int, int]]
    ) -> tuple[int, int]:
        """
        Given a section name and static section spans, determine the span of the section

        Parameters
        ----------
        section_name: str,
            the name of the section to calculate the span for
        alignment_section_spans: dict[str, tuple[int, int]]
            the mapping of alignment section names to their (start, stop) indices in the sequence

        Returns
        -------
        tuple[int, int]
            the (start, stop) indices of the section in the sequence
        """
        _span = (100000, 0)  # init to extreme values
        for alignment_section, (align_start, align_stop) in alignment_section_spans.items():
            adj_start, adj_stop = self._adjustment_table[alignment_section][section_name]
            if adj_start < 0:
                _span = (
                    max(min(align_start + adj_start + 1, _span[0]), 0),
                    max(align_start + adj_stop + 1, _span[1]),
                )
            else:
                _span = (
                    max(min(align_stop + adj_start - 1, _span[0]), 0),
                    max(align_stop + adj_stop - 1, _span[1]),
                )
        return _span


class LibraryLocator(SectionSequenceMapper):
    """
    Given a read and known static barcode sections, locate the library section

    This is a specialized mapper that only locates the library section
    based on the known static barcode sections. This is useful for
    demultiplexing reads to libraries quickly, without needing to
    align all the required sections.
    """

    def _get_sections_to_map(self) -> list[str]:
        """Only need to map library section for this one"""
        return ["library"]

    def get_library_seq(self, sequence: SequenceRecord, alignment_section_spans: dict[str, tuple[int, int]]) -> str:
        """
        Get the library sequence from the aligned sequence based on the provided alignment section spans

        This will avoid aligning the rest of the required sections and just return the library to save
        some time if that is all that is needed.

        Parameters
        ----------
        sequence: SequenceRecord
            the sequence to get the library from
        alignment_section_spans: dict[str, tuple[int, int]]
            the mapping of alignment section names to their (start, stop) indices in the sequence
        """
        _library_span = self.calculate_section_span("library", alignment_section_spans)
        return sequence.sequence[_library_span[0] : _library_span[1]]

    def iter_possible_library_seqs(
        self, sequence: SequenceRecord, alignment_section_spans: dict[str, tuple[int, int]]
    ) -> Iterator[str]:
        """
        Iterate over possible library sequences based on the provided alignment section spans

        It is possible if using more than one static section that the library
        codon is not perfectly defined due to INDELs in the read. It's also possible
        that an INDEL elsewhere cause a misalignment. This will wiggle the library
        codon by one base in either direction to generate possible library sequences
        of length equal to the library tag length or one less (in that order)

        This is very helpful for calling the library, as it helps account for
        error that occur outside the library codon region

        Parameters
        ----------
        sequence: SequenceRecord
            the sequence to get the library from
        alignment_section_spans: dict[str, tuple[int, int]]
            the mapping of alignment section names to their (start, stop) indices in the sequence

        Yields
        ------
        str
            possible library barcodes
        """
        library_span = self.get_library_seq(sequence, alignment_section_spans)
        library_tag_length = len(self.barcode_schema.library_section.section_tag)
        for length in [library_tag_length, library_tag_length - 1]:
            start = 0
            end = start + length
            while end <= len(library_span):
                yield library_span[start:end]
                start += 1
                end = start + length


class SectionSequenceAligner(SectionSequenceMapper):
    """
    A sequence aligner based on static section locations

    The alignment algorithm uses the positions of known static barcode
    sections to locate the region that other sections should be found
    in the observed sequence. This works by pre-compling a map of
    "adjustments" from the static sections to the required sections
    based on the distance and direction between the sections.
    This enables rapid alignment of sequences based on the static sections.

    The aligner is only based on the distances between sections and their
    names, not the nucleotides in them. Therefore, if all your DELs in a
    given collection have the same sections, with the same lengths you
    only need one aligner object for all of them. You can figure out this
    by using the `is_schema_compatible` method on the `BarcodeSchema`.

    If multiple static sections are used, they might disagree slightly on
    where the section is due to INDELs in the read. In this case, the
    union of their predicted spans is used. For example if two static
    sections predict the same required section to be at (10, 20) and (12, 22),
    the final predicted span will be (10, 22).

    Notes
    -----
    Will only align the sections listed as required by the
    `get_required_sections` method of the `BarcodeSchema`.
    """

    def _get_sections_to_map(self) -> list[str]:
        """Map to only the required sections minus the library tag (already located)"""
        required = self.barcode_schema.get_required_section_names()
        required.remove("library")  # not need for this library
        return required

    def align_sequence(
        self, sequence: SequenceRecord, alignment_section_spans: dict[str, tuple[int, int]]
    ) -> Iterator[tuple[AlignedSeq, float]]:
        """
        Align a sequence based on the provided alignment section spans

        Parameters
        ----------
        sequence: SequenceRecord
            the sequence to align
        alignment_section_spans: dict[str, tuple[int, int]]
            the mapping of alignment section names to their (start, stop) indices in the sequence
            determined by the alignment algorithm

        Returns
        -------
        AlignedSeq
        """
        _span_map = self._span_map.copy()
        for required_section in self._required_sections:
            _span_map[required_section] = self.calculate_section_span(required_section, alignment_section_spans)
        return iter([(AlignedSeq(sequence=sequence, section_spans=_span_map), 0)])


class BarcodeAligner:
    """
    Aligns sequences to the fully barcode reference using pairwise alignment

    This will take a sequence and run a pairwise alignment against a library's
    full barcode schema. This means any static sections (including overhangs)
    will be used to align the sequence. This is better at detecting INDELs,
    recovering more reads compared to a static alignment, but is also more
    expensive.

    Notes
    -----
    DELi uses a semi-global alignment with a custom substitution matrix to
    allow for variable 'N' in the read to match any base with 0 score penalty.
    It is implemented using Biopython's PairwiseAligner.

    Numba is also used to optimize the coordinate mapping from the alignment,
    so the first time a library aligner is used could result in a 1-2 second
    compile overhead (but only once)

    Parameters
    ----------
    barcode_schema: BarcodeSchema
        the barcode schema to use for alignment
    include_library_section: bool
        whether to include the library section in the required sections
        for alignment. Default is True.
    """

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
        self._required_sections = barcode_schema.get_required_section_names()
        if not include_library_section:
            self._required_sections.remove("library")

    def align_sequence(self, sequence: SequenceRecord) -> Iterator[tuple[AlignedSeq, float]]:
        """
        Align a sequence to the barcode reference

        Since there could be many alignment with the same score,
        this will yield all alignments found with their score.

        Notes
        -----
        These alignments are done on the fly as requested, allowing
        for early stopping if a perfect alignment is found early on

        While the alignment code is will optimized (written in C)
        the coordinate mapping is done in raw Python. It is
        accelerated a bit using numba, but it is still a bottleneck.

        Parameters
        ----------
        sequence: SequenceRecord
            the sequence to align

        Yields
        ------
        tuple[AlignedSeq, float]
        """
        aligned_sec_spans: dict[str, tuple[int, int]] = dict()
        alignments = self.aligner.align(sequence.sequence, self.barcode_ref)
        for alignment in alignments:
            _reverse_idx = _query_to_ref_map(alignment.coordinates)
            for barcode_sec in self._required_sections:
                sec_span = self.spans[barcode_sec]
                aligned_sec_spans[barcode_sec] = (
                    int(_reverse_idx[sec_span.start]),
                    int(_reverse_idx[sec_span.stop - 1] + 1),
                )
            yield AlignedSeq(sequence, aligned_sec_spans), alignment.score


class FailedStaticAlignment(FailedDecodeAttempt):
    """A sequence that failed to align to a static reference"""

    def __init__(self, sequence: SequenceRecord):
        super().__init__(sequence, "Failed to locate static barcode regions")


class FailedLibraryBarcodeLookup(FailedDecodeAttempt):
    """A sequence that failed to align to a static reference"""

    def __init__(self, sequence: SequenceRecord):
        super().__init__(sequence, "Failed to lookup the library codon")


class _Query(Generic[M_co], abc.ABC):
    """
    Abstract class for queries used to demultiplex reads based on libraries

    These queries are used to determine which libraries are possible
    for a given read, and then provide a static alignment to locate
    the library barcode (and other required sections) to be called.

    Notes
    -----
    Query objects are covariant with the respective Match object

    Parameters
    ----------
    libraries: list[TaggedLibrary]
        the libraries covered by this query
    section_names: tuple[str, ...]
        the names of the static sections to use for querying
    """

    def __init__(self, libraries: list[TaggedLibrary], section_names: tuple[str, ...]):
        self.libraries: list[TaggedLibrary] = libraries
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
                        f"Barcode section name {section_name} not found in query covered library {library.library_id}"
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
        self._trim_adjustments: dict[TaggedLibrary, tuple[int, int]] = dict()
        for library in self.libraries:
            _min_sec_start = 100000
            _max_sec_end = 0
            for section_name in section_names:
                static_section_span = library.barcode_schema.section_spans[section_name]
                _min_sec_start = min(_min_sec_start, static_section_span.start)
                _max_sec_end = max(_max_sec_end, static_section_span.stop)
            start_adj = min(library.barcode_schema.required_start - _min_sec_start, 0)
            end_adj = max(library.barcode_schema.required_end - _max_sec_end, 0)
            self._trim_adjustments[library] = (start_adj, end_adj)

    @abc.abstractmethod
    def __hash__(self):
        """Queries must be hashable"""
        raise NotImplementedError()

    @abc.abstractmethod
    def search(self, sequence: str) -> Optional[M_co]:
        """
        Given a sequence, search for a match to the query

        Parameters
        ----------
        sequence: str
            the sequence to search

        Returns
        -------
        Optional[M_co]
            the match found, or None if no match
        """
        raise NotImplementedError()


class _Match(Generic[Q_co], abc.ABC):
    """
    Abstract class for matches found by a query

    Matches provide the ability to determine where the
    static sections used in the query to made it are
    located in the sequence read.

    Notes
    -----
    Match objects are covariant with the respective Query object

    Parameters
    ----------
    query: Q_co
        the query that produced this match
    score: float | int
        the score of the match
    """

    def __init__(self, query: Q_co, score: float | int):
        self.query = query
        self.score = score

    @abc.abstractmethod
    def get_section_spans(self) -> dict[str, tuple[int, int]]:
        """
        Get the spans of the sections from the match

        Returns
        -------
        dict[str, tuple[int, int]]
            mapping of section names to their (start, stop) indices in the sequence
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def trim_seq(self, library_call: ValidLibraryCall, sequence: SequenceRecord) -> SequenceRecord:
        """
        Trim a sequence based on the match and library call

        Parameters
        ----------
        library_call: ValidLibraryCall
            the library call to use for trimming

        sequence: SequenceRecord
            the sequence to trim

        Returns
        -------
        SequenceRecord
            the trimmed sequence
        """
        raise NotImplementedError()


class RegexError(Exception):
    """Exception raised for errors in the regex demultiplexing process."""

    pass


class _RegexQuery(_Query["_RegexMatch"]):
    """
    Regex based query for demultiplexing reads based on libraries

    Regex queries will build a fuzzy match group for each static section, and
    a variable length wildcard between them based on the distance between the
    static sections. A buffer of +-5 is used to compensate for possible INDELs.
    Static sections will also have match group names idential to their section names.
    An example looks like:
    `(?:(?P<primer1>AGCTAGCT)){e<=1}.{{10,20}}(?:(?P<primer2>CGTACGTA)){e<=1}`

    Notes
    -----
    By default Regex will always look for the best match, not the first
    viable one

    Parameters
    ----------
    libraries: list[TaggedLibrary]
        the libraries covered by this query
    section_names: tuple[str, ...]
        the names of the static sections to use for querying
    error_tolerance: int
        the number of errors to allow per static section

    Attributes
    ----------
    pattern: str
        the regex pattern used for matching
    compiled_pattern: Pattern
        the compiled regex pattern used for matching
    """

    def __init__(self, libraries: list[TaggedLibrary], section_names: tuple[str, ...], error_tolerance: int = 1):
        super().__init__(libraries, section_names)

        # build the regex string
        self.pattern = ""
        for i, (section_name, sequence) in enumerate(self._sec_name2seq.items()):
            self.pattern += f"(?:(?P<{section_name}>{sequence})){{e<={error_tolerance}}}"
            if i < len(self._sec_name2seq) - 1:  # flanked by other sections
                key = (section_names[i], section_names[i + 1])
                min_dist, max_dist = self._dist_between_section[key]
                self.pattern += f".{{{min_dist - 5},{max_dist + 5}}}"
        self.compiled_pattern: Pattern = compile(self.pattern, BESTMATCH)

    def __hash__(self):
        """A regex query is hashed based on its pattern"""
        return hash(self.pattern)

    def search(self, sequence: str) -> Optional["_RegexMatch"]:
        """
        Search for the pattern in the sequence

        Parameters
        ----------
        sequence: str
            the sequence to search

        Returns
        -------
        Optional[_RegexMatch]
            the match found, or None if no match
        """
        match = self.compiled_pattern.search(sequence)
        return _RegexMatch(self, sum(match.fuzzy_counts), match) if match is not None else None


class _RegexMatch(_Match[_RegexQuery]):
    """
    A match found by a regex query

    Notes
    -----
    Somewhat confusingly there is also a `regex.Match` object that
    is returned after a regex Pattern is searched. This object
    wraps that object to provide additional functionality specific
    to DELi.

    Parameters
    ----------
    query: _RegexQuery
        the query that produced this match
    score: float | int
        the score of the match
    match: Match
        the regex match object
    """

    def __init__(self, query: _RegexQuery, score: float | int, match: Match):
        super().__init__(query, score)
        self.match = match

    def get_section_spans(self) -> dict[str, tuple[int, int]]:
        return {
            section_name: span for section_name, span in zip(self.query.section_names, self.match.regs[1:], strict=True)
        }

    def trim_seq(self, library_call: ValidLibraryCall, sequence: SequenceRecord) -> SequenceRecord:
        start_adj, end_adj = self.query._trim_adjustments[library_call.obj]
        # build the trim
        trim = slice(
            max(0, min(self.match.regs[1][0] - end_adj, library_call.obj.barcode_schema.required_start) - 10),
            max(self.match.regs[-1][1] + end_adj, library_call.obj.barcode_schema.required_end) + 10,
        )
        return sequence[trim]


class CutAdaptError(Exception):
    """Exception raised for errors in the cutadapt demultiplexing process."""

    pass


class _CutadaptQuery(_Query["_CutadaptMatch"]):
    """
    A Cutadapt based query for demultiplexing reads based on libraries

    This will use the cutadapt hybrid alignment algorithm to locate
    static sections in the read. Depending on where the static sections
    are located in the barcode, different adapter types will be used
    to optimize performance. If a static section is located near the
    front (within 10% of the barcode length) a FrontAdapter will be used.
    If it is located near the back (within 10% of the barcode length from
    the end) a RightmostBackAdapter will be used. If it is located in the
    middle (more than 10% from either end) an AnywhereAdapter will be used.
    Cutadapt queries can only support up to 2 static sections. If 2 sections are used,
    they will be linked together using a LinkedAdapter. You can read more
    about this from the cutadapt docs:
    https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types

    Parameters
    ----------
    libraries: list[TaggedLibrary]
        the libraries covered by this query
    section_names: tuple[str, ...]
        the names of the static sections to use for querying
    error_tolerance: int
        the number of errors to allow per static section
    min_overlap: int
        the minimum overlap required for a match
    """

    def __init__(
        self,
        libraries: list[TaggedLibrary],
        section_names: tuple[str, ...],
        error_tolerance: int = 1,
        min_overlap: int = 8,
    ):
        if len(section_names) > 2:
            raise CutAdaptError(
                f"Cutadapt based queries can only handle up to 2 static sections; found {len(section_names)} sections"
            )
        super().__init__(libraries, section_names)

        self._hash_name = "---".join(self._sec_name2seq.values())

        # determine which type of adapter to build for each section
        adapters: list[SingleAdapter] = []
        for section in section_names:
            distance_from_front_percent: float = 0.0
            distance_from_back_percent: float = 0.0
            for library in self.libraries:
                span = library.barcode_schema.section_spans[section]
                start = span.start
                end = span.stop
                distance_from_front_percent = max(distance_from_front_percent, start / len(library.barcode_schema))
                distance_from_back_percent = max(
                    distance_from_back_percent, (len(library.barcode_schema) - end) / len(library.barcode_schema)
                )

            if ((distance_from_front_percent > 0.1) and (distance_from_back_percent > 0.1)) or (
                (distance_from_front_percent <= 0.1) and (distance_from_back_percent <= 0.1)
            ):
                warnings.warn(
                    "Cutadapt based queries may not perform well when a given static section "
                    "is located in the middle or on both sides of the barcode for libraries in the query",
                    stacklevel=1,
                )
                adapters.append(
                    AnywhereAdapter(
                        sequence=self._sec_name2seq[section],
                        max_errors=error_tolerance,
                        min_overlap=min_overlap,
                        name=section,
                    )
                )
            elif distance_from_front_percent <= 0.1:
                adapters.append(
                    FrontAdapter(
                        sequence=self._sec_name2seq[section],
                        max_errors=error_tolerance,
                        min_overlap=min_overlap,
                        name=section,
                    )
                )
            else:
                adapters.append(
                    RightmostBackAdapter(
                        sequence=self._sec_name2seq[section],
                        max_errors=error_tolerance,
                        min_overlap=min_overlap,
                        name=section,
                    )
                )

        # if two adapters, link them together
        self.adapter: LinkedAdapter | SingleAdapter
        if len(adapters) == 2:
            self.adapter = LinkedAdapter(
                front_adapter=adapters[0],
                back_adapter=adapters[1],
                front_required=True,
                back_required=True,
                name=adapters[0].name,
            )
        else:
            self.adapter = adapters[0]

    def __hash__(self):
        """A cutadapt query is hashed based on its adapter sequences"""
        return hash(self._hash_name)

    def search(self, sequence: str) -> "_CutadaptMatch | None":
        """
        Search for the pattern in the sequence using cutadapt

        Parameters
        ----------
        sequence: str

        Returns
        -------
        _CutadaptMatch | None
            the match found, or None if no match
        """
        match: LinkedMatch | SingleMatch | None = self.adapter.match_to(sequence)
        if match is None:
            return None
        elif isinstance(match, LinkedMatch):
            return _LinkedCutadaptMatch(self, match.errors, match)
        else:
            return _SingleCutadaptMatch(self, match.errors, match)


class _CutadaptMatch(_Match[_CutadaptQuery], abc.ABC):
    """
    A base class for cutadapt based matches

    Cutadapt will return different match objects (also called Matches confusingly)
    depending on the type of adapter used. This base class is for any possible
    cutadapt match to help with type hinting.
    """

    match: OGCutadaptMatch

    def __init__(self, query: _CutadaptQuery, score: float | int):
        super().__init__(query, score)


class _LinkedCutadaptMatch(_CutadaptMatch):
    """
    A match found by a cutadapt linked adapter query

    Parameters
    ----------
    query: _CutadaptQuery
        the query that produced this match
    score: float | int
        the score of the match
    match: LinkedMatch
        the cutadapt linked match object
    """

    def __init__(self, query: _CutadaptQuery, score: float | int, match: LinkedMatch):
        super().__init__(query, score)
        self.match: LinkedMatch = match

    def get_section_spans(self) -> dict[str, tuple[int, int]]:
        front_match = self.match.front_match
        back_match = self.match.back_match

        return {
            front_match.adapter.name: (front_match.rstart, front_match.rstop),
            back_match.adapter.name: (back_match.rstart + front_match.rstop, back_match.rstop + front_match.rstop),
        }

    def trim_seq(self, library_call: ValidLibraryCall, sequence: SequenceRecord) -> SequenceRecord:
        start_adj, end_adj = self.query._trim_adjustments[library_call.obj]
        # build the trim
        trim_start = max(
            0, min(self.match.front_match.rstart - end_adj, library_call.obj.barcode_schema.required_start) - 10
        )
        trim_end = max(self.match.back_match.rstop + end_adj, library_call.obj.barcode_schema.required_end) + 10
        return sequence[trim_start:trim_end]


class _SingleCutadaptMatch(_CutadaptMatch):
    """
    A match found by a cutadapt single adapter query

    Parameters
    ----------
    query: _CutadaptQuery
        the query that produced this match
    score: float | int
        the score of the match
    match: SingleMatch
        the cutadapt single match object
    """

    def __init__(self, query: _CutadaptQuery, score: float | int, match: SingleMatch):
        super().__init__(query, score)
        self.match: SingleMatch = match

    def get_section_spans(self) -> dict[str, tuple[int, int]]:
        return {
            self.match.adapter.name: (self.match.rstart, self.match.rstop),
        }

    def trim_seq(self, library_call: ValidLibraryCall, sequence: SequenceRecord) -> SequenceRecord:
        start_adj, end_adj = self.query._trim_adjustments[library_call.obj]
        trim_start = max(0, min(self.match.rstart - end_adj, library_call.obj.barcode_schema.required_start) - 10)
        trim_end = max(self.match.rstop + end_adj, library_call.obj.barcode_schema.required_end) + 10
        return sequence[trim_start:trim_end]


class LibraryDemultiplexerError(Exception):
    """Exception raised for errors in the library demultiplexing process."""

    pass


class LibraryDemultiplexer(abc.ABC):
    """
    Base class for all library demultiplexers objects

    Library demultiplexers take in a collection of libraries
    and provide the ability to both:
    1. Determine which library a given sequence belongs to
    2. Provide an alignment(s) of the sequence to the barcode schema of that library

    In DELi we call these two combined steps "library demultiplexing"
    """

    def __init__(
        self,
        libraries: DELibraryCollection,
        tool_compounds: Optional[list[TaggedToolCompoundLibrary]],
        revcomp: bool = True,
    ):
        self.revcomp = revcomp
        self.tool_compounds = tool_compounds if tool_compounds is not None else []
        self.libraries = libraries

        self.all_libraries: list[TaggedLibrary] = list(libraries.libraries) + self.tool_compounds

    @abc.abstractmethod
    def demultiplex(
        self, sequence: SequenceRecord
    ) -> tuple[ValidLibraryCall, Iterator[tuple[AlignedSeq, float]]] | FailedDecodeAttempt:
        """
        Demultiplex a sequence into a library call and alignments

        This both determine the library the sequence belongs to, as
        generate a set of alignments of the sequence to the library's
        barcode schema, mapping the barcode sections to the read



        Parameters
        ----------
        sequence: SequenceRecord
            the sequence to demultiplex

        Returns
        -------
        tuple[ValidCall[DELibrary | ToolCompoundLibrary], Iterator[tuple[AlignedSeq, float]]] | FailedDecodeAttempt
            the library call and alignments, or a failed decode attempt
            if demultiplexing failed
        """
        raise NotImplementedError()


class QueryBasesLibraryDemultiplexer(LibraryDemultiplexer, Generic[Q], abc.ABC):
    """
    Base class for a query based library demultiplexer

    Query based demultiplexers generate a set of queries
    based on settings and the static sections of the libraries
    provided. These queries are then used to generate an
    alignment and library call for a given sequence.

    Query bases demultiplexers can optionally use dynamic
    aligners to realign the barcode region after the library
    has been called.

    Also, depending on the use case, static section aligners
    and or dynamic aligners might need to be built. This is
    controlled using the `build_dynamic_aligners_` and
    `build_section_seq_aligners_` parameters. They cannot
    be set by the users, and are only here to help with
    class inheritance and construction at the development
    side.

    Attributes
    ----------
    _library_aligners: dict[DELibrary, BarcodeAligner]
        the dynamic barcode aligners for each library.
        will be a unique object for every library
        will be empty if `build_dynamic_aligners_` is False
    _section_seq_aligners: dict[DELibrary, DELSectionSequenceAligner]
        the static section sequence aligners for each library
        may be shared between libraries if they have
        compatible barcode schemas
        will be empty if `build_section_seq_aligners_` is False
    _library_locaters: dict[Q, list[LibraryLocators]]
        the required library locators for each query needed to
        correctly locate the library barcode for all libraries
        the respective query covers
    _queries: list[Q]
        the queries used for demultiplexing
    """

    def __init__(
        self,
        *args,
        realign: bool = False,
        build_dynamic_aligners_: bool = False,
        build_section_seq_aligners_: bool = False,
        build_library_locaters_: bool = False,
        **kwargs,
    ):
        self.realign = realign
        super().__init__(*args, **kwargs)

        self._library_aligners: dict[TaggedLibrary, BarcodeAligner] = dict()
        self._section_seq_aligners: dict[TaggedLibrary, SectionSequenceAligner] = dict()
        self._library_locaters: dict[Q, list[LibraryLocator]] = dict()
        self._queries: list[Q] = self._get_queries()

        if build_dynamic_aligners_ or self.realign:  # only initialize these if we need them
            self._library_aligners = self._get_library_aligners()

        if build_section_seq_aligners_:  # only initialize these if we need them
            self._section_seq_aligners = self._get_section_seq_aligners()

        if build_library_locaters_:
            self._library_locaters = self._get_library_locators()

    @abc.abstractmethod
    def _get_queries(self) -> list[Q]:
        """
        Get the queries used for demultiplexing

        This function should be implemented by subclasses
        to generate the appropriate queries based on
        the libraries and settings provided
        """
        raise NotImplementedError()

    def _get_library_aligners(self) -> dict[TaggedLibrary, BarcodeAligner]:
        """
        Get the dynamic library aligners for each library

        Returns
        -------
        dict[TaggedLibrary, BarcodeAligner]
            the dynamic library aligners for each library
        """
        aligners: dict[TaggedLibrary, BarcodeAligner] = {}
        for library in self.all_libraries:
            _found_existing: bool = False
            for other_library, aligner in aligners.items():
                if library.barcode_schema == other_library.barcode_schema:
                    aligners[library] = aligner
                    _found_existing = True
                    break
            if not _found_existing:
                aligners[library] = BarcodeAligner(library.barcode_schema)

        return aligners

    def _get_section_seq_aligners(self) -> dict[TaggedLibrary, SectionSequenceAligner]:
        """
        Get the static section sequence aligners for each library

        Since some libraries may have compatible barcode schemas,
        they can share section sequence aligners. This function will
        save memory by reusing section sequence aligners where possible.

        Returns
        -------
        dict[DELibrary, DELSectionSequenceAligner]
            the static section sequence aligners for each library
        """
        aligners: dict[TaggedLibrary, SectionSequenceAligner] = {}
        for query in self._queries:
            for library in query.libraries:
                _found_existing: bool = False
                for other_library, aligner in aligners.items():
                    if library.barcode_schema.is_schema_align_compatible(other_library.barcode_schema):
                        aligners[library] = aligner
                        _found_existing = True
                        break
                if not _found_existing:
                    aligners[library] = SectionSequenceAligner(library.barcode_schema, query.section_names)
        return aligners

    def _get_library_locators(self) -> dict[Q, list[LibraryLocator]]:
        """
        Get the library locators for each query

        For any given query, this will look at the libraries it covers
        and then determine the minimum set of library locators needed
        to cover all libraries to avoid duplicate work during
        decoding.

        Returns
        -------
        dict[_Query, list[LibraryLocators]
            map of queries to the library locators that need to be
            used during decoding
        """
        locators: dict[Q, list[LibraryLocator]] = {}
        for query in self._queries:
            query_locators: list[LibraryLocator] = []
            lib_locator_map: dict[TaggedLibrary, LibraryLocator] = {}
            for library in query.libraries:
                _covered = False
                for other_library, _ in lib_locator_map.items():
                    if library.barcode_schema.is_static_library_locate_compatible(
                        other_library.barcode_schema, list(query.section_names)
                    ):
                        _covered = True
                        break
                if not _covered:
                    new_locator = LibraryLocator(library.barcode_schema, query.section_names)
                    query_locators.append(new_locator)
                    lib_locator_map[library] = new_locator
        return locators

    def _get_best_match(self, sequence: SequenceRecord) -> _Match | None:
        """
        Get the best match for a sequence among all possible queries

        Will exit early if a perfect match is found, otherwise will
        search all queries to find the best match. If no match is found,
        None is returned.

        If revcomp mode is used, will first search the read as if its the positive
        sequence. If no perfect match is found, will then search the reverse complement
        for a possible better (or perfect) match

        Parameters
        ----------
        sequence: SequenceRecord
            the sequence to search

        Returns
        -------
        _Match | None
            the best match found, or None if no match
        """
        best_match: _Match | None = None
        best_score: float = float("inf")

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
            best_revcomp_score: float = float("inf")

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

    def _make_alignments(
        self,
        best_match: _Match,
        library_call: ValidLibraryCall,
        sequence: SequenceRecord,
        _static_alignment_spans: Optional[dict[str, tuple[int, int]]] = None,
    ) -> Iterator[tuple[AlignedSeq, float]]:
        """
        Given the best match and the library call it generated, produce alignments

        Will either use dynamic realignment if specified,
        otherwise it will reuse the alignment(s) generated to
        do the library call.

        Parameters
        ----------
        best_match: _Match
            the best match found for the sequence
        library_call: ValidLibraryCall
            the library call generated from the best match
        sequence: SequenceRecord
            the sequence to align
        _static_alignment_spans: Optional[dict[str, tuple[int, int]]]
            precomputed static alignment spans to use for static alignment
            if not provided will generate them from the best match.
            This is mainly here to avoid recomputing them if already done
            during the library call process.
        """
        if self.realign:
            # this will trim the sequence to the region we care about with a 10bp buffer on each side
            trimmed_seq = best_match.trim_seq(library_call, sequence)
            realigner: BarcodeAligner = self._library_aligners[library_call.obj]
            return realigner.align_sequence(trimmed_seq)
        else:
            _static_alignment_spans = (
                best_match.get_section_spans() if _static_alignment_spans is None else _static_alignment_spans
            )
            return self._section_seq_aligners[library_call.obj].align_sequence(sequence, _static_alignment_spans)


class LibraryTagLibraryDemultiplexer(QueryBasesLibraryDemultiplexer[Q], abc.ABC):
    """
    Base class for library tag based library demultiplexers

    Library tag demultiplexers use *only* the library tag section
    to determine which library a sequence belongs to. This is useful
    in cases where there are no (known) static sections in the barcode schema.

    It is generally rare that this is the case, and library tag demultiplexers
    are both more expensive (more queries needed) and less accurate (fewer nucleotides
    to match) than other query based demultiplexers.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, build_section_seq_aligners_=True, **kwargs)

    def demultiplex(
        self, sequence: SequenceRecord
    ) -> tuple[ValidLibraryCall, Iterator[tuple[AlignedSeq, float]]] | FailedDecodeAttempt:
        """See :meth:`~deli.decode.library_demultiplexer.LibraryDemultiplexer.demultiplex`."""
        best_match = self._get_best_match(sequence)

        # failed to match any static region query
        if best_match is None:
            return FailedStaticAlignment(sequence)
        else:
            library_call = ValidCall(best_match.query.libraries[0], best_match.score)
            return library_call, self._make_alignments(best_match, library_call, sequence)


class LibraryTagRegexLibraryDemultiplexer(LibraryTagLibraryDemultiplexer[_RegexQuery]):
    """
    A library tag based library demultiplexer using regex queries

    Will build a separate regex query for each library using that
    library's tag.

    Parameters
    ----------
    libraries: DELibraryCollection
        the libraries to demultiplex
    tool_compounds: Optional[list[TaggedToolCompoundLibrary]], default = None
        the tool compounds to look for during decoding
    error_tolerance: int, default = 1
        the number of errors to allow in the library tag section
        these can be INDELs or SNPs
    realign: bool, default = False
        whether to realign the barcode region using a semi global
        alignment after the library is called
    revcomp: bool, default = True
        whether to also try the reverse complement of the read
        when searching for matches
    """

    def __init__(
        self,
        libraries: DELibraryCollection,
        tool_compounds: Optional[list[TaggedToolCompoundLibrary]] = None,
        error_tolerance: int = 1,
        realign: bool = False,
        revcomp: bool = True,
    ):
        self.error_tolerance = error_tolerance
        super().__init__(libraries=libraries, tool_compounds=tool_compounds, realign=realign, revcomp=revcomp)

    def _get_queries(self) -> list[_RegexQuery]:
        """
        Get the regex queries for library tag demultiplexing

        Will generate a regex query for each library using that library's tag section
        as the only static section in the query.
        """
        queries: list[_RegexQuery] = list()
        for library in self.all_libraries:
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
                    error_tolerance=self.error_tolerance,
                )
            )
        return queries


class LibraryTagCutadaptLibraryDemultiplexer(LibraryTagLibraryDemultiplexer[_CutadaptQuery]):
    """
    A library tag based library demultiplexer using cutadapt queries

    WARNING: Unless you library tag is very close to the front or back of the barcode,
    cutadapt based library tag demultiplexing may perform poorly. This is because
    cutadapt optimized to find terminal adapters, and will choose matches closer to
    the ends *even if* there are better matches in the middle of the read.

    Will build a separate cutadapt query for each library using that
    library's tag.

    Parameters
    ----------
    libraries: DELibraryCollection
        the libraries to demultiplex
    tool_compounds: Optional[list[TaggedToolCompoundLibrary]], default = None
        the tool compounds to look for during decoding
    error_tolerance: int, default = 1
        the number of errors to allow in the library tag section
        these can be INDELs or SNPs
    realign: bool, default = False
        whether to realign the barcode region using a semi global
        alignment after the library is called
    min_overlap: int, default = 8
        the minimum overlap required for a match
        see the cutadapt docs for more details
    revcomp: bool, default = True
        whether to also try the reverse complement of the read
    """

    def __init__(
        self,
        libraries: DELibraryCollection,
        tool_compounds: Optional[list[TaggedToolCompoundLibrary]] = None,
        error_tolerance: int = 1,
        realign: bool = False,
        min_overlap: int = 8,
        revcomp: bool = True,
    ):
        self.error_tolerance = error_tolerance
        self.min_overlap = min_overlap
        super().__init__(libraries=libraries, tool_compounds=tool_compounds, realign=realign, revcomp=revcomp)

    def _get_queries(self) -> list[_CutadaptQuery]:
        """
        Get the cutadapt queries for library tag demultiplexing

        Will generate a cutadapt query for each library using that library's tag section
        as the only static section in the query. These will be SingleAdapters
        """
        queries: list[_CutadaptQuery] = list()
        for library in self.all_libraries:
            queries.append(
                _CutadaptQuery(libraries=[library], section_names=("library",), error_tolerance=self.error_tolerance)
            )
        return queries


class NoLibraryTagLibraryDemultiplexer(QueryBasesLibraryDemultiplexer[Q], abc.ABC):
    """
    Base class for no library tag based library demultiplexers

    As the name suggests, these demultiplexers use static sections that
    are not the library tag to demultiplex reads. These can be very
    efficient and accurate when there are multiple static sections
    shared by many libraries in the run.

    Since a given query in this set up might cover more than one
    library, a `BarcodeCaller` is used to determine which library.
    This works by using the matched static regions to align the
    sequence and extract the possible library tag sequence.
    This library caller supports error correction.
    """

    def __init__(self, *args, error_correction_mode_str: str = "levenshtein_dist:2,asymmetrical", **kwargs):
        super().__init__(*args, build_section_seq_aligners_=True, build_library_locaters_=True, **kwargs)

        # maps queries to barcode callers for the libraries in that query
        self._library_callers: dict[Q, BarcodeCaller[DELibrary | TaggedToolCompoundLibrary]] = {
            query: get_barcode_caller(
                {library.library_tag: library for library in query.libraries},
                error_correction_mode_str=error_correction_mode_str,
            )
            for query in self._queries
        }

    @abc.abstractmethod
    def _get_query_object(self, libraries: list[TaggedLibrary], section_names: tuple[str, ...]) -> Q:
        """
        Initialize the correct query object for the demultiplexer
        """
        raise NotImplementedError()

    def demultiplex(
        self, sequence: SequenceRecord
    ) -> tuple[ValidLibraryCall, Iterator[tuple[AlignedSeq, float]]] | FailedDecodeAttempt:
        """See :meth:`~deli.decode.library_demultiplexer.LibraryDemultiplexer.demultiplex`."""
        best_match = self._get_best_match(sequence)

        # failed to match any static region query
        if best_match is None:
            return FailedStaticAlignment(sequence)

        # these are the spans of the static regions from the match
        _static_alignment_spans = best_match.get_section_spans()

        # Now find the best matching library within the best regex query
        library_call: ValidLibraryCall | FailedBarcodeLookup = FailedBarcodeLookup("placeholder")
        if len(best_match.query.libraries) == 1:
            library_call = ValidCall(best_match.query.libraries[0], 0.0)  # only one possible library
        else:
            _best_call_score = 100.0  # lower scores are better

            # This could probably be optimized by avoiding re-extraction of the same sequences
            # if we know that the static sections are the same distance from the library tags
            # across several aligners. Instead of doing that now, we just track which sequences
            # we've already tried to avoid some redundant work.
            observed_possible_lib_seqs: set[str] = set()

            for init_aligner in self._library_locaters[best_match.query]:
                _best_internal_call: ValidLibraryCall | FailedBarcodeLookup = FailedBarcodeLookup("placeholder")
                _best_internal_score = 100.0
                # allow some flexability in the library tag location
                for possible_library_seq in init_aligner.iter_possible_library_seqs(sequence, _static_alignment_spans):
                    # track which sequences we've already tried
                    if possible_library_seq in observed_possible_lib_seqs:
                        continue  # skip if we've already tried this sequence
                    else:
                        observed_possible_lib_seqs.add(possible_library_seq)

                    # try and call this one
                    _library_call = self._library_callers[best_match.query].decode_barcode(possible_library_seq)
                    if not isinstance(_library_call, FailedBarcodeLookup):
                        if _library_call.score < _best_internal_score:
                            _best_internal_call = _library_call
                            _best_internal_score = _library_call.score
                    if _best_internal_score == 0:
                        library_call = _best_internal_call
                        break  # perfect match found

                if _best_internal_call is not None:  # library locator made a call
                    if _best_internal_score == 0:
                        library_call = _best_internal_call
                        break  # perfect match found
                    elif _best_internal_score < _best_call_score:
                        library_call = _best_internal_call
                        _best_call_score = _best_internal_score

        # at this point we did everything we can to try and call the library
        if isinstance(library_call, FailedBarcodeLookup):
            return FailedLibraryBarcodeLookup(sequence)
        elif self.realign:
            return library_call, self._library_aligners[library_call.obj].align_sequence(sequence)
        else:
            return library_call, self._make_alignments(
                best_match, library_call, sequence, _static_alignment_spans=_static_alignment_spans
            )


class SinglePrimerLibraryDemultiplexer(NoLibraryTagLibraryDemultiplexer[Q], abc.ABC):
    """
    Demultiplexer for single primer style demultiplexing

    This will use a single static section (generally a primer) to:
    1. determine which libraries the sequence could possibly match to
    2. determine all possible sequences the library tag could have
       given the location of the static sections and barcode schemas
    3. call the library
    4. determine a final alignment of the read to the given library barcode schema
    """

    def __init__(self, *args, error_tolerance: int = 1, **kwargs):
        self.error_tolerance = error_tolerance
        super().__init__(*args, **kwargs)

    def _get_queries(self) -> list[Q]:
        """
        Generate the regex query sets for the single primer demultiplexer

        There are alot of ways to pick the possible static sections
        Right now, DELi does this on it own with no user input.
        It uses a greedy algorithm to pick the minimal set of static sections
        that cover all libraries in the run.
        This may not be "optimal" in the sense that the barcode or the most
        dissimilar or large to help avoid error.
        But should be good enough for most use cases.
        """
        # map all possible static barcode to the barcode schemas they cover
        static_seq_groups: defaultdict[tuple[str, str], set[int]] = defaultdict(set)
        for idx, library in enumerate(self.all_libraries):
            static_section_tags = [
                (sec.get_dna_sequence(), sec.section_name) for sec in library.barcode_schema.static_sections
            ]
            for sec_seq, sec_name in static_section_tags:
                static_seq_groups[(sec_seq, sec_name)].add(idx)

        # greedily pick sequences that cover the most remaining barcode schemas
        candidates = static_seq_groups.copy()
        uncovered = set(range(len(self.all_libraries)))
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
                raise RuntimeError("this is impossible, something went wrong. Please raise an issue")
            else:
                new_idxes = set(candidates[best_seq]) & uncovered
                libraries: list[TaggedLibrary] = list([self.all_libraries[i] for i in new_idxes])

                chosen.append(self._get_query_object(libraries=libraries, section_names=(best_seq[1],)))

                uncovered -= static_seq_groups[best_seq]
                del candidates[best_seq]
        return chosen


class FlankingPrimersLibraryDemultiplexer(NoLibraryTagLibraryDemultiplexer[Q], abc.ABC):
    """
    Demultiplexer for flanking primer style demultiplexing

    This is very similar to single primer demultiplexing, but uses a second
    static region that must be on the opposite side of the library tag.
    This can help improve accuracy of the decoding at a small cost to
    speed and memory.

    This will use the flanking static section (generally a primer) to:
    1. determine which libraries the sequence could possibly match to
    2. determine all possible sequences the library tag could have
       given the location of the static sections and barcode schemas.
       It will use both static sections to do this to help correct potential errors
    3. call the library
    4. determine a final alignment of the read to the given library barcode schema
    """

    def __init__(self, *args, error_tolerance: int = 1, **kwargs):
        self.error_tolerance = error_tolerance
        super().__init__(*args, **kwargs)

    def _get_queries(self) -> list[Q]:
        """
        Generate the regex query sets for the flanking primer demultiplexer

        This is a greedy algorithm to pick the minimal set of static section pairs
        that cover all libraries in the run.
        This may not be "optimal" in the sense that the barcode or the most
        dissimilar or large to help avoid error.

        A static section pair are two static sections, one before
        and one after the library barcode section.
        """
        # map all possible static barcode to the barcode schemas they cover
        static_seq_groups: defaultdict[tuple[tuple[str, str], tuple[str, str]], set[int]] = defaultdict(set)
        for idx, library in enumerate(self.all_libraries):
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

        # greedily pick sequences that cover the most remaining barcode schemas
        candidates = static_seq_groups.copy()
        uncovered = set(range(len(self.all_libraries)))
        chosen: list[Q] = []
        while uncovered:
            best_seq: tuple[tuple[str, str], tuple[str, str]] | None = None
            best_cover: int = 0
            for _sec_seq, grp in candidates.items():
                cover = len(grp & uncovered)
                if cover > best_cover:
                    best_cover = cover
                    best_seq = _sec_seq
            if best_seq is None:
                raise RuntimeError("this is impossible, something went wrong. Please raise an issue")
            else:
                new_idxes = set(candidates[best_seq]) & uncovered
                libraries: list[TaggedLibrary] = [self.all_libraries[i] for i in new_idxes]

                chosen.append(
                    self._get_query_object(
                        libraries=libraries,
                        section_names=tuple(list(sec_name for sec_name, _ in best_seq)),
                    )
                )

                uncovered -= static_seq_groups[best_seq]
                del candidates[best_seq]
        return chosen


class SinglePrimerCutadaptLibraryDemultiplexer(SinglePrimerLibraryDemultiplexer[_CutadaptQuery]):
    """
    Demultiplexer for single primer cutadapt style demultiplexing

    WARNING: Unless your flanking primers are very close to the front and back of the barcode,
    cutadapt based demultiplexing may perform poorly.

    See the docs for more details on how this algorithm works.

    Parameters
    ----------
    libraries: DELibraryCollection
        the libraries to use for demultiplexing
    tool_compounds: Optional[list[TaggedToolCompoundLibrary]], default = None
        the tool compounds to look for during decoding
    revcomp: bool, default = True
        whether to revcomp the sequences when demultiplexing
    realign: bool, default = False
        whether to realign the sequences after demultiplexing
    error_tolerance: int, default = 1
        the number of errors to allow in the static section
        these can be INDELs or SNPs
    error_correction_mode_str: str, default = levenshtein_dist:2,asymmetrical
        the error correction mode to use for library tag decoding
    min_overlap: int, default = 8
        the minimum overlap required for a match
        see the cutadapt docs for more details
    """

    def __init__(
        self,
        libraries: DELibraryCollection,
        tool_compounds: Optional[list[TaggedToolCompoundLibrary]] = None,
        revcomp: bool = True,
        realign: bool = False,
        error_tolerance: int = 1,
        error_correction_mode_str: str = "levenshtein_dist:2,asymmetrical",
        min_overlap: int = 8,
    ):
        self.min_overlap = min_overlap
        super().__init__(
            libraries=libraries,
            tool_compounds=tool_compounds,
            realign=realign,
            error_correction_mode_str=error_correction_mode_str,
            error_tolerance=error_tolerance,
            revcomp=revcomp,
        )

    def _get_query_object(self, libraries: list[TaggedLibrary], section_names: tuple[str, ...]) -> _CutadaptQuery:
        """Get the correct cutadapt query object for single primer demultiplexing"""
        return _CutadaptQuery(
            libraries=libraries,
            section_names=section_names,
            error_tolerance=self.error_tolerance,
            min_overlap=self.min_overlap,
        )


class FlankingPrimersCutadaptLibraryDemultiplexer(FlankingPrimersLibraryDemultiplexer[_CutadaptQuery]):
    """
    Demultiplexer for flanking primer cutadapt style demultiplexing

    WARNING: Unless your flanking primers are very close to the front and back of the barcode,
    cutadapt based demultiplexing may perform poorly.

    See the docs for more details on how this algorithm works.

    Parameters
    ----------
    libraries: DELibraryCollection
        the libraries to use for demultiplexing
    tool_compounds: Optional[list[TaggedToolCompoundLibrary]], default = None
        the tool compounds to look for during decoding
    revcomp: bool, default = True
        whether to revcomp the sequences when demultiplexing
    realign: bool, default = False
        whether to realign the sequences after demultiplexing
    error_tolerance: int, default = 1
        the number of errors to allow in the static section
        these can be INDELs or SNPs
    error_correction_mode_str: str, default = levenshtein_dist:2,asymmetrical
        the error correction mode to use for library tag decoding
    min_overlap: int, default = 8
        the minimum overlap required for a match
        see the cutadapt docs for more details
    """

    def __init__(
        self,
        libraries: DELibraryCollection,
        tool_compounds: Optional[list[TaggedToolCompoundLibrary]] = None,
        revcomp: bool = True,
        realign: bool = False,
        error_tolerance: int = 1,
        error_correction_mode_str: str = "levenshtein_dist:2,asymmetrical",
        min_overlap: int = 8,
    ):
        self.min_overlap = min_overlap
        super().__init__(
            libraries=libraries,
            tool_compounds=tool_compounds,
            realign=realign,
            error_correction_mode_str=error_correction_mode_str,
            error_tolerance=error_tolerance,
            revcomp=revcomp,
        )

    def _get_query_object(self, libraries: list[TaggedLibrary], section_names: tuple[str, ...]) -> _CutadaptQuery:
        """Get the correct cutadapt query object for flanking primer demultiplexing"""
        return _CutadaptQuery(
            libraries=libraries,
            section_names=section_names,
            error_tolerance=self.error_tolerance,
            min_overlap=self.min_overlap,
        )


class SinglePrimerRegexLibraryDemultiplexer(SinglePrimerLibraryDemultiplexer[_RegexQuery]):
    """
    Demultiplexer for single primer regex style demultiplexing

    See the docs for more details on how this algorithm works.

    Parameters
    ----------
    libraries: DELibraryCollection
        the libraries to use for demultiplexing
    tool_compounds: Optional[list[TaggedToolCompoundLibrary]], default = None
        the tool compounds to look for during decoding
    revcomp: bool, default = True
        whether to revcomp the sequences when demultiplexing
    realign: bool, default = False
        whether to realign the sequences after demultiplexing
    error_tolerance: int, default = 1
        the number of errors to allow in the static section
        these can be INDELs or SNPs
    error_correction_mode_str: str, default = levenshtein_dist:2,asymmetrical
        the error correction mode to use for library tag decoding
    """

    def __init__(
        self,
        libraries: DELibraryCollection,
        tool_compounds: Optional[list[TaggedToolCompoundLibrary]] = None,
        realign: bool = False,
        revcomp: bool = True,
        error_tolerance: int = 1,
        error_correction_mode_str: str = "levenshtein_dist:1,asymmetrical",
    ):
        super().__init__(
            libraries=libraries,
            tool_compounds=tool_compounds,
            realign=realign,
            error_correction_mode_str=error_correction_mode_str,
            error_tolerance=error_tolerance,
            revcomp=revcomp,
        )

    def _get_query_object(self, libraries: list[TaggedLibrary], section_names: tuple[str, ...]) -> _RegexQuery:
        """Get the correct regex query object for single primer demultiplexing"""
        return _RegexQuery(
            libraries=libraries,
            section_names=section_names,
            error_tolerance=self.error_tolerance,
        )


class FlankingPrimersRegexLibraryDemultiplexer(FlankingPrimersLibraryDemultiplexer[_RegexQuery]):
    """
    Demultiplexer for flanking primer regex style demultiplexing

    See the docs for more details on how this algorithm works.

    Parameters
    ----------
    libraries: DELibraryCollection
        the libraries to use for demultiplexing
    tool_compounds: Optional[list[TaggedToolCompoundLibrary]], default = None
        the tool compounds to look for during decoding
    revcomp: bool, default = True
        whether to revcomp the sequences when demultiplexing
    realign: bool, default = False
        whether to realign the sequences after demultiplexing
    error_tolerance: int, default = 1
        the number of errors to allow in the static section
        these can be INDELs or SNPs
    error_correction_mode_str: str, default = levenshtein_dist:2,asymmetrical
        the error correction mode to use for library tag decoding
    """

    def __init__(
        self,
        libraries: DELibraryCollection,
        tool_compounds: Optional[list[TaggedToolCompoundLibrary]] = None,
        revcomp: bool = True,
        realign: bool = False,
        error_tolerance: int = 1,
        error_correction_mode_str: str = "levenshtein_dist:1,asymmetrical",
    ):
        super().__init__(
            libraries=libraries,
            tool_compounds=tool_compounds,
            realign=realign,
            error_correction_mode_str=error_correction_mode_str,
            error_tolerance=error_tolerance,
            revcomp=revcomp,
        )

    def _get_query_object(self, libraries: list[TaggedLibrary], section_names: tuple[str, ...]) -> _RegexQuery:
        """Get the correct regex query object for flanking primer demultiplexing"""
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

    This Demultiplexer will align every sequence to every library's
    full barcode schema and call the library that matches best.

    WARNING: This is very slow. It should only be used in scenarios where
    compute is nearly infinite and free and you want to squeeze every
    ounce of data out of your selection *or* your reads are so riddle with
    errors that other demultiplexing methods fail.

    Parameters
    ----------
    libraries: DELibraryCollection
        the libraries to use for demultiplexing
    tool_compounds: Optional[list[TaggedToolCompoundLibrary]], default = None
        the tool compounds to look for during decoding
    revcomp: bool, default = False
        whether to also consider the reverse complement of the sequence
        Note: This is guaranteed to double the amount of work done, since it will
        need to align every read twice. Using an external tool to standardize
        the orientation of all reads before demultiplexing is recommended.

    Attributes
    ----------
    aligners: list[BarcodeAligner]
        the aligners for each library used to align the full sequence
    """

    def __init__(
        self,
        libraries: DELibraryCollection,
        tool_compounds: Optional[list[TaggedToolCompoundLibrary]],
        revcomp: bool = False,
    ):
        super().__init__(libraries=libraries, tool_compounds=tool_compounds, revcomp=revcomp)

        self.aligners = [
            BarcodeAligner(library.barcode_schema, include_library_section=False) for library in self.libraries
        ]

    def demultiplex(
        self, sequence: SequenceRecord
    ) -> tuple[ValidLibraryCall, Iterator[tuple[AlignedSeq, float]]] | FailedDecodeAttempt:
        """See :meth:`~deli.decode.library_demultiplexer.LibraryDemultiplexer.demultiplex`."""
        best_alignments: tuple[ValidLibraryCall, Iterator[tuple[AlignedSeq, float]]] | None = None
        best_score = 0.0
        for library, aligner in zip(self.libraries, self.aligners, strict=False):
            alignments = aligner.align_sequence(sequence)
            alignment, score = alignments.__next__()
            if score > best_score:
                best_alignments = (
                    ValidCall(library, score=0),
                    itertools.chain(
                        iter(
                            [
                                (alignment, score),
                            ]
                        ),
                        alignments,
                    ),
                )
                best_score = score

        # biopython alignments do not track errors; it won't be possible to know
        # if the match that exists already is perfect with processing the alignment
        # but that processing part is so expensive that it is not worth doing here
        if self.revcomp:
            revcomp_seq = sequence.reverse_complement()
            best_revcomp_alignments: tuple[ValidLibraryCall, Iterator[tuple[AlignedSeq, float]]] | None = None
            best_revcomp_score = 0.0
            for library, aligner in zip(self.libraries, self.aligners, strict=False):
                alignments = aligner.align_sequence(revcomp_seq)
                alignment, score = alignments.__next__()
                if score > best_revcomp_score:
                    best_revcomp_alignments = (
                        ValidCall(library, score=0),
                        itertools.chain(
                            iter(
                                [
                                    (alignment, score),
                                ]
                            ),
                            alignments,
                        ),
                    )
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


######################################
# Numba accelerated helper functions #
######################################


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


NoLibraryTagLibraryDemultiplexer.demultiplex.__doc__ = LibraryDemultiplexer.demultiplex.__doc__
LibraryTagLibraryDemultiplexer.demultiplex.__doc__ = LibraryDemultiplexer.demultiplex.__doc__
FullSeqAlignmentLibraryDemultiplexer.demultiplex.__doc__ = LibraryDemultiplexer.demultiplex.__doc__
