"""matching sequences to DEL barcode schemas"""

import enum
import re
from itertools import count
from typing import Iterable, List, Tuple

import regex

from deli.dels import BarcodeSchema
from deli.sequence import FastqSequence
from deli.sequence.fasta import FastaSequence
from deli.constants import FAST_MATCH_PADDING


# global match_id value
MATCH_ID = count()


class MatchOutcome:
    """Abstract class that defines the outcome of a match"""

    def __init__(self):
        self.match_id = next(MATCH_ID)


class NoSequenceMatch(MatchOutcome):
    """
    Match result if there were no matches in the sequence

    Parameters
    ----------
    sequence: FastaSequence
        the sequence searched
    """

    def __init__(self, sequence: FastaSequence):
        super().__init__()
        self.sequence = sequence


class SequenceTooBig(MatchOutcome):
    """
    Match result if the sequence was too big to match

    Parameters
    ----------
    sequence: FastaSequence
        the sequence searched

    Attributes
    ----------
    sequence_length: int
        the length of the sequence
    """

    def __init__(self, sequence: FastaSequence):
        super().__init__()
        self.sequence = sequence
        self.sequence_length = len(sequence)


class SequenceTooSmall(MatchOutcome):
    """
    Match result if the sequence was too small to match

    Parameters
    ----------
    sequence: FastaSequence
        the sequence searched

    Attributes
    ----------
    sequence_length: int
        the length of the sequence
    """

    def __init__(self, sequence: FastaSequence):
        super().__init__()
        self.sequence = sequence
        self.sequence_length = len(sequence)


class FastBarcodeMatch(MatchOutcome):
    """
    Stores minimal information about a match between a sequence and barcode schema

    Attributes
    ----------
    match_id: int
        a unique id for the match
    """

    def __init__(
        self,
        sequence: FastaSequence,
        match_span: Tuple[int, int],
    ):
        """
        Initialize a FastBarcodeMatch object

        Parameters
        ----------
        sequence: FastaSequence
            the full sequence of the read this sequence match is from
        match_span: tuple[int, int]
            the start and end index of the match in the
            `sequence` in format [start, end)
        """
        super().__init__()
        self.sequence = sequence
        self.match_span = match_span

    @property
    def match_sequence(self) -> FastaSequence:
        """The sequence of the match"""
        return self.sequence.subset(self.match_span[0], self.match_span[1])


class FullBarcodeMatch(FastBarcodeMatch):
    """
    Stores all information about a match between a sequence and barcode schema

    Attributes
    ----------
    match_id: int
        a unique id for the match
    substitute_error_count: int
        the number of substitute errors found during matching
    insert_error_count: int
        the number of insert errors found during matching
    delete_error_count: int
        the number of delete errors found during matching
    """

    def __init__(
        self,
        sequence: FastaSequence,
        match_span: Tuple[int, int],
        substitute_error_idx: List[int],
        insert_error_idx: List[int],
        delete_error_idx: List[int],
    ):
        """
        Initialize a FullBarcodeMatch object

        Parameters
        ----------
        sequence: FastaSequence
            the full sequence of the read this sequence match is from
        match_span: tuple[int, int]
            the start and end index of the match in the
            `sequence` in format [start, end)
        substitute_error_idx: list[int]
            the index of each substitute error in the `sequence`
        insert_error_idx: list[int]
            the index of each insert error in the `sequence`
        delete_error_idx: list[int]
            the index of each delete error in the `sequence`
        """
        super().__init__(sequence=sequence, match_span=match_span)
        self.sequence = sequence
        self.match_span = match_span
        self.substitute_error_count = len(substitute_error_idx)
        self.substitute_error_idx = substitute_error_idx
        self.insert_error_count = len(insert_error_idx)
        self.insert_error_idx = insert_error_idx
        self.delete_error_count = len(delete_error_idx)
        self.delete_error_idx = delete_error_idx

    def total_errors(self) -> int:
        """
        Get the total number of errors made during matching

        Returns
        -------
        int
        """
        return sum([self.substitute_error_count, self.insert_error_count, self.delete_error_count])

    def get_match_sequence_substitute_error_idx(self) -> List[int]:
        """
        Get the index positions of each substitute error in the `match_sequence`

        Notes
        -----
        `self.substitute_error_idx` contains the indexes with respect to the `sequence`

        Returns
        -------
        List[int]
            the index position of each substitute error in the `match_sequence`
        """
        return [idx - self.match_span[0] for idx in self.substitute_error_idx]

    def get_match_sequence_insert_error_idx(self) -> List[int]:
        """
        Get the index positions of each insert error in the `match_sequence`

        Notes
        -----
        `self.insert_error_idx` contains the indexes with respect to the `sequence`

        Returns
        -------
        List[int]
            the index position of each insert error in the `match_sequence`
        """
        return [idx - self.match_span[0] for idx in self.insert_error_idx]

    def get_match_sequence_delete_error_idx(self) -> List[int]:
        """
        Get the index positions of each delete error in the `match_sequence`

        Notes
        -----
        `self.substitute_error_idx` contains the indexes with respect to
        the `sequence`

        Returns
        -------
        List[int]
            the index position of each delete error in the `match_sequence`
        """
        return [idx - self.match_span[0] for idx in self.delete_error_idx]

    def get_sorted_match_sequence_error_idx(
        self,
        exclude_substitute: bool = False,
        exclude_insert: bool = False,
        exclude_delete: bool = False,
        normalize: bool = False,
    ) -> Tuple[List[int], List[str]]:
        """
        Get error positions in the `match_sequence` sorted by location.

        Notes
        -----
        Will also return a second list containing the type of error at
         each position ('substitute', 'insert', 'delete')

        Parameters
        ----------
        exclude_substitute: bool, default=False
            if True, exclude the substitute error positions in the return
        exclude_insert: bool, default=False
            if True, exclude the insert error positions in the return
        exclude_delete: bool, default=False
            if True, exclude the delete error positions in the return
        normalize: bool, default=False
            if True, normalize the error positions such that `0` is the
            first element of the match, rather than the full sequence

        Returns
        -------
        error_positions, error_types: Tuple[List[int], List[str]]
            a tuple with element[0] being the positions list
            element[1] is a list of error_type mapped to the position list

        Raises
        ------
        RuntimeError:
            raised when all `exclude_` parameters are True
        """
        _errors = []
        _error_types = []

        if not exclude_substitute:  # add substitute errors to
            _error_pos = self.get_match_sequence_substitute_error_idx()
            _errors.extend(_error_pos)
            _error_types.extend(["substitute"] * len(_error_pos))
        if not exclude_insert:  # add substitute errors to
            _error_pos = self.get_match_sequence_insert_error_idx()
            _errors.extend(_error_pos)
            _error_types.extend(["insert"] * len(_error_pos))
        if not exclude_delete:  # add substitute errors to
            _error_pos = self.get_match_sequence_delete_error_idx()
            _errors.extend(_error_pos)
            _error_types.extend(["delete"] * len(_error_pos))

        _sort_mask = sorted(range(len(_errors)), key=_errors.__getitem__)

        if normalize:
            return (
                [_errors[i]-self.match_span[0] for i in _sort_mask],
                [_error_types[i] for i in _sort_mask]
            )
        else:
            return (
                [_errors[i] for i in _sort_mask],
                [_error_types[i] for i in _sort_mask]
            )

class MatchModes(enum.Enum):
    """
    Enum for matching mode settings

    Notes
    -----
    The string values attached to this MUST match the
    string values used to define that section in the
    schema format. See the schema docs for more info
    """

    preindex = "pre-index"
    primer = "primer"
    overhang = "_overhang"
    preumi = "pre-umi"
    closing_primer = "closing_primer"


class BarcodeMatchPatternError(Exception):
    """Exception for failure to create a matching pattern"""

    pass


class BarcodeMatcher:
    """
    Matches sequences to a DEL barcode schema

    Attributes
    ----------
    pattern: str
        the regex pattern used to find matches

    Notes
    -----
    When building the matching pattern based on the
    barcode schema, there are two sub-patterns:
    - variable patterns
    - static patterns

    static patterns are based on static regions of
    the barcode. Variable patterns are based on size of
    the variable region and match any character

    So a barcode with a primer region of "AGCTGCT" and
    3 building block tag regions of size 5 each, you might
    have a pattern of:

    "AGCTGCT.{5}.{5}.{5}"

    Any static region from the barcode schema not used in
    pattern will be treated as a variable region

    There are several matching modes that can be used.
    Multiple modes can be used, and they will just be combined
    - pre-index:
        this mode will look for the static pre-index region in
        the barcode schema
    - primer:
        this mode will look for the static primer region in
        the barcode schema
    - overhang:
        this mode will look for the overhang regions of all
        DEL tag regions in the barcode schema
    - pre-umi:
        this mode will look for the static pre-umi region in
        the barcode schema
    - closing_primer:
        this mode will look for the static closing_primer region in
        the barcode schema
    - all:
        use all the modes at once. This will return FullMatch results
        and can be used to do quick but unforgiving calling

    It should be noted that the more modes you use, the more static
    regions are looked for. This could be good or bad depending on
    how you are doing your calling.
    Generally, if you are using the align ment caller, you want only
    1 or 2 modes (fewer is better).
    If you are using the match caller you want as many modes as
    possible.

    Matching also has the options for exact or error tolerate matching
    If using exact matching, no base errors in the static regions
    are tolerated, exact matches are needed. If not using exact matching
    you can specify the amount of error you will tolerate in the static
    regions. Errors include insertion, deletion and substitution.

    There is no way to specify insert or deletion errors in the variable
    regions of the pattern; they are covered with error tolerance in
    the static regions OR the padding

    Padding is used to add "junk" DNA to the front and end of the
    sequence. This can prevent the patter from failing to match
    because a variable region is looking for more DNA than exists
    at the head or tail. If not using exact matching, its recommend
    to use padding.
    """

    def __init__(
        self,
        barcode_scheme: BarcodeSchema,
        exact_match: bool = False,
        full_match: bool = False,
        match_index: bool = False,
        match_modes: Iterable[MatchModes] = (MatchModes.primer,),
        error_tolerance: int = 3,
        padding: int = 5,
        rev_comp: bool = True,
    ):
        """
        Initialize BarcodeMatcher object

        Parameters
        ----------
        barcode_scheme: BarcodeSchema
            barcode schema to match too
        exact_match: bool, default=False
            match the static regions of the barcode schema
            exactly
        full_match: bool, default=False
            use all barcode section to make a match
            will override any passed match modes
        match_index: bool, default=False
            search for/match the index regions
        match_modes: Iterable[MatchModes], default=MatchModes.primer
            the matching modes to use when building the pattern
        error_tolerance: int, default=3
            if not exact_match, how many error are tolerated in
            static regions
        padding: int, default=5
            if not exact_match, add padding to the being or end
            to prevent unintentional trimming of the start/end
        rev_comp: bool, default=True
            also search in the reverse complement of the barcode
        """
        self.barcode_scheme = barcode_scheme
        self.exact_match = exact_match
        self.full_match = full_match
        self.match_index = match_index
        self.error_tolerance = error_tolerance
        self.padding = padding
        self.rev_comp = rev_comp

        # if full match, set match modes to all modes:
        self.match_modes = [
            MatchModes.preindex, MatchModes.primer,
            MatchModes.overhang, MatchModes.closing_primer,
            MatchModes.preumi
        ] if self.full_match else match_modes


        if self.match_index and (not self.barcode_scheme.has_index()):
            raise ValueError("`match_index` cannot be `True` when `barcode_scheme` lacks index")

        self.pattern = self._make_pattern()

    def _check_modes(self):
        """check to see if chosen modes are present in barcode"""
        pass

    def _make_pattern(self):
        _pattern = ""

        def _add_static_region(section_name, _pattern):
            """Helper function to add static barcode region to pattern"""
            if self.exact_match:
                _pattern += self.barcode_scheme[section_name]
            else:
                _pattern += (
                    f"(?:{self.barcode_scheme[section_name]})" 
                    f"{{e<={self.error_tolerance}}}"
                )
            return _pattern

        def _add_variable_region(length, _pattern):
            _pattern += f".{{{length}}}"
            return _pattern

        _running_wildcard_count = 0
        for region in self.barcode_scheme.barcode_sections.keys():
            if any([re.search(match_mode.value, region) for match_mode in self.match_modes]):
                if _running_wildcard_count > 0:
                    _pattern = _add_variable_region(_running_wildcard_count, _pattern)
                    _running_wildcard_count = 0
                _pattern = _add_static_region(region, _pattern)
            else:
                _running_wildcard_count += self.barcode_scheme.barcode_lengths[region]

        # add tailing wildcards
        #  not sure if we want this yet....
        if _running_wildcard_count > 0:
            _pattern = _add_variable_region(_running_wildcard_count, _pattern)

        if len({"A", "G", "C", "T"}.intersection(set(_pattern))) == 0:
            raise BarcodeMatchPatternError(f"pattern contains no static regions:\n{_pattern}")

        return _pattern

    def match(
        self,
        sequences: List[FastqSequence],
        min_sequence_size: int = 75,
        max_sequence_size: int = 1000,
    ) -> List[MatchOutcome]:
        """
        Given a set of sequence and a pattern, find all matches of the pattern in the sequence

        Notes
        -----
        Will always return a list of matches even if only 1 or match is found.

        Will loop through sequence to find all non-overlapping matches.
        The search is iterative and will move from start to end of the sequence
        Thus the first match will be found at the front of the sequence.
        This match will then be removed and the sequence searched again for
        another mathc (if the sequence is long enough).
        This is useful for sequencers like the NanoPore that can mistake two
        sequences as a single sequence.

        Matching can fail for 3 reasons:
        1. Sequence was larger than the `max_sequence_size`
        2. Sequence was smaller than the `min_sequence_size`
        3. regex failed to find a match to the barcode schema pattern

        All passed sequence reads will always return a match object,
        even if they fail to find a match for any reason
        This is to assist with logging so that we can trace back which reads
        failed to match and check to see why, and manually review sequences
        if need be

        Parameters
        ----------
        sequences: List[FastqSequence]
            sequence to search for matches in
            can be FastqSequence object or the full sequence as a string
        min_sequence_size: int, default=70
            the minimum length of a sequence to be searched for a match
            cutoff is exclusive
        max_sequence_size: int, default=1000
            the maximum length of a sequence to be searched for a match
            if above this size will not be searched at all
            cutoff is exclusive

        Returns
        -------
        List[BarcodeMatch]
            the list of all matches
            if not matches found will be empty list
        """
        all_matches = []

        if self.rev_comp:
            sequences.extend([sequence.revcomp() for sequence in sequences])

        for sequence in sequences:
            # add padding to both ends of sequence
            if self.padding > 0:
                _sequence = "X" * self.padding + sequence.sequence + "X" * self.padding
            else:
                _sequence = sequence.sequence

            _matches: List[MatchOutcome] = list()

            # check sequence for correct size
            if len(_sequence) > max_sequence_size:
                _matches = [SequenceTooBig(sequence=sequence)]
            if len(_sequence) < min_sequence_size:
                _matches = [SequenceTooSmall(sequence=sequence)]

            # loop until sequence is too small or no matches found
            while len(_sequence) >= min_sequence_size:
                # skip best match flag
                if self.exact_match or self.full_match:
                    hit = regex.search(self.pattern, _sequence, flags=regex.BESTMATCH)
                else:
                    hit = regex.search(self.pattern, _sequence)

                if hit:
                    _sequence = _sequence[
                        (hit.regs[0][1] + 1) :
                    ]  # remove current hit from sequence
                    # add full match if exact mode is used
                    if self.full_match:
                        _matches.append(
                            FullBarcodeMatch(
                                sequence=sequence,
                                match_span=(
                                    hit.regs[0][0] - self.padding,
                                    hit.regs[0][1] - self.padding + (
                                            len(hit.fuzzy_changes[1]) -
                                            len(hit.fuzzy_changes[2])
                                    )
                                ),
                                substitute_error_idx=[
                                    x - self.padding for x in hit.fuzzy_changes[0]
                                ],
                                insert_error_idx=[
                                    x - self.padding for x in hit.fuzzy_changes[1]
                                ],
                                delete_error_idx=[
                                    x - self.padding for x in hit.fuzzy_changes[2]
                                ],
                            )
                        )
                    # fast match is exact matching is not used
                    else:
                        _matches.append(
                            FastBarcodeMatch(
                                sequence=sequence,
                                match_span=(
                                    hit.regs[0][0] - self.padding - FAST_MATCH_PADDING,
                                    hit.regs[0][1] - self.padding + FAST_MATCH_PADDING
                                )
                            )
                        )
                else:
                    break  # exit loop if no matches

            # if no matches found in sequence add a no match object
            if len(_matches) == 0:
                _matches = [NoSequenceMatch(sequence=sequence)]
            all_matches.extend(_matches)
        return all_matches
