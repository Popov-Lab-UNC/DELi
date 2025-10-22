"""code for calling the building blocks of DEL reads"""

import abc
import math
from collections import defaultdict
from typing import TypeAlias, Union

from deli._hamming import BaseQuaternaryHamming
from deli.configure import DELiConfigError, get_deli_config
from deli.dels import BuildingBlockBarcodeSection, TaggedBuildingBlock, TaggedBuildingBlockSet

from .calls import FailedCall, ValidCall


class ValidBuildingBlockCall(ValidCall):
    """base class for all valid building block calls"""

    def __init__(self, building_block: TaggedBuildingBlock, score: Union[float, int]):
        self.building_block = building_block
        self.score = score


class FailedBuildingBlockCall(FailedCall):
    """Base class for all failed building block calls"""

    pass


# type alias for building block calling
BuildingBlockCall: TypeAlias = Union[ValidBuildingBlockCall, FailedBuildingBlockCall]


class ErrorCorrectorException(Exception):
    """base class for all error corrector exceptions"""

    pass


class ErrorCorrector(abc.ABC):
    """Base class for all error correctors"""

    @abc.abstractmethod
    def correct_sequence(self, observed_seq: str) -> str | None:
        """
        Given an observed observed_seq correct it to the correct to the expected observed_seq

        Parameters
        ----------
        observed_seq: str
            the observed DNA observed_seq

        Returns
        -------
        str | None
            the corrected DNA observed_seq; None if fails to correct
        """

    pass


class DummyErrorCorrector(ErrorCorrector):
    """An 'error corrector' that always just return the given sequence"""

    def correct_sequence(self, observed_seq: str) -> str | None:
        """
        Given an observed observed_seq correct it to the correct to the expected observed_seq

        Parameters
        ----------
        observed_seq: str
            the observed DNA observed_seq

        Returns
        -------
        str | None
            the corrected DNA observed_seq; None if fails to correct
        """
        return observed_seq


class HammingDecodeError(Exception):
    """raised when the hamming decoder fails to decode a barcode"""

    pass


class QuaternaryHammingDecoder(BaseQuaternaryHamming, ErrorCorrector):
    """Generates and decodes quaternary hamming codes"""

    def __init__(self, parity_map: list[int], has_extra_parity: bool):
        """
        Initialize a QuaternaryHammingDecoder

        Parameters
        ----------
        parity_map: list[int]
            a list defining where the correctly ordered hamming bits are
        has_extra_parity: bool
            True if the sequences being decoded has extra parity
        """
        super().__init__()

        self.parity_map = parity_map
        self.has_extra_parity = has_extra_parity

        self.real_size = len(self.parity_map)
        self.hamming_size = (2 ** math.ceil(math.log2(self.real_size))) - self.has_extra_parity

        self._parity = math.ceil(math.log2(self.hamming_size))

        self._require_sort = sorted(self.parity_map) != self.parity_map

    def __eq__(self, other) -> bool:
        """Return true if the encoder has the same total length"""
        if isinstance(other, self.__class__):
            return self.parity_map == other.parity_map
        else:
            return False

    @classmethod
    def load(cls, name) -> "QuaternaryHammingDecoder":
        """
        Load a hamming code into a decoder

        NOTE: to load, need a deli.hamming.<name> section
        in the ".deli" config file with both a hamming and custom order field.
        See docs for more details

        Parameters
        ----------
        name: str
            name of the hamming code to load, e.g. "8_4"

        Returns
        -------
        QuaternaryHammingDecoder
        """
        _config = get_deli_config()
        _true_order, _custom_order = _config.get_hamming_code(name)

        if 0 in _true_order:
            _has_extra_parity = True
        else:
            _has_extra_parity = False
            _custom_order = [_ - 1 for _ in _true_order]

        return cls(parity_map=_custom_order, has_extra_parity=_has_extra_parity)

    def correct_sequence(self, observed_seq: str) -> str | None:
        """
        Given an observed observed_seq correct it to the correct to the expected observed_seq

        Parameters
        ----------
        observed_seq: str
            the observed DNA observed_seq

        Returns
        -------
        str | None
            the corrected DNA observed_seq; None if fails to correct
        """
        _tag = [self.nuc_2_int_mapping[char] for char in observed_seq]
        try:
            _decoded_tag = self._hamming_decode(_tag)
            return "".join([self.int_2_nuc_mapping[_] for _ in _decoded_tag])
        except HammingDecodeError:
            return None

    def _hamming_decode(self, bases: list[int]) -> list[int]:
        """
        Decode a tag and correct any correctable errors

        Parameters
        ----------
        bases: Union[List[Base], Tag]
            the tag to decode/correct

        Raises
        ------
        HammingDecodeError
            when more than 2 errors are detected in the barcode
            and it is not correctable

        Returns
        -------
        Tag
        """
        # sort hamming vector if custom order used
        if self._require_sort:
            bases = [bases[i] for i in self.parity_map]

        # need to pad with 0 to next power of 2 length
        _actual_length = len(bases)
        _padded_bases = bases + [0] * (self.hamming_size - _actual_length)

        if not self.has_extra_parity:
            _padded_bases = [0] + _padded_bases

        sub_parity = [
            self._get_sub_parity(_padded_bases, order=p) for p in range(self._parity, 0, -1)
        ][::-1]

        if self.has_extra_parity:
            _parity_set_length = len(set(sub_parity))
            if _parity_set_length > 2:
                raise HammingDecodeError("2 or more errors detected")
            if (_parity_set_length == 2) and (0 not in sub_parity):
                raise HammingDecodeError("2 or more errors detected")

            overall_parity = self._get_global_parity(_padded_bases)
            if (overall_parity == 0) and (len(set(sub_parity)) == 2):
                raise HammingDecodeError("2 or more errors detected")

            if max(sub_parity) != 0:
                self._correct_error(_padded_bases, sub_parity)

            elif overall_parity != 0:
                # need to correct the parity bit in the case that this is the one that is misread
                _padded_bases[0] = (4 - sum(_padded_bases[1:])) % 4

            # final check for overall parity
            overall_parity = self._get_global_parity(_padded_bases)
            if overall_parity != 0:
                raise HammingDecodeError("2 or more errors detected")

        # if not extra parity
        else:
            if max(sub_parity) != 0:
                self._correct_error(_padded_bases, sub_parity)

        return _padded_bases[not self.has_extra_parity : _actual_length]


class HashMapCollisionError(Exception):
    """raised when a hash map error corrector fails to build due to collisions"""

    pass


class HashMapErrorCorrector(ErrorCorrector, abc.ABC):
    """Base class for all hash map based error correctors"""

    def __init__(
        self, bb_set: TaggedBuildingBlockSet, distance_cutoff: int = 1, asymmetrical: bool = False
    ):
        """
        Initialize a HashMapErrorCorrector

        This works by taking each building block tag, creating all its neighbors within a given
        distance,
        and mapping them to the building block itself. This way, a query can be looked up and
        "corrected"
        to its original building block tag in constant time.

        There are some assumptions made here though. First, it assumes that there are no
        collisions. That
        is, no two building blocks have the same tag or the same neighbors within the given
        distance.
        If this is the case, we cannot determine with 100% certainty which building block is
        the correct
        original tag. Building the hash map will fail if a collision is detected by default.

        However, you can set `asymmetrical=True` to allow the hash map to return the best,
        non-ambiguous match instead, thus when building collisions are ignored in favor of
         saving
        the closest tag. Ties are ambiguous, so the hash map will return None.

        Hamming3 or Levenshtein3 distance schemes guarantee all tags
        have a min distance between them (in this case 3), thus the distance cutoff is
        also known:
        (min_distance - 1) / 2. In this case it is easier to construct the hash map by
        using this info.
        Only use asymmetrical when you are sure that the tags follow a known distance
        encoding scheme.

        Parameters
        ----------
        bb_set: TaggedBuildingBlockSet
            the building block set to build the hash map from
        """
        self.distance_cutoff = distance_cutoff
        self.asymmetrical = asymmetrical
        self._hash_map: dict[str, tuple[TaggedBuildingBlock | None, int]] = {}

        for bb in bb_set:
            possible_tags = self._get_neighbors(bb.tag)
            for tag, dist in possible_tags.items():
                if tag in self._hash_map:
                    if not self.asymmetrical:
                        _colliding_bb = self._hash_map[tag][0]
                        if _colliding_bb is None:
                            # this is for type hinting and mypy
                            raise RuntimeError(
                                "this should never happen, please report to the developers"
                            )
                        else:
                            raise HashMapCollisionError(
                                f"collision detected when building {self.__class__.__name__} "
                                f"with distance {self.distance_cutoff} for tag {tag}; "
                                f"collision between building blocks {bb.bb_id} and "
                                f"{_colliding_bb.bb_id} with tags {bb.tag} "
                                f"and {_colliding_bb.tag} respectively"
                            )
                    else:
                        if self._hash_map[tag][1] == dist:
                            self._hash_map[tag] = (None, dist)
                        elif self._hash_map[tag][1] > dist:
                            self._hash_map[tag] = (bb, dist)
                else:
                    self._hash_map[tag] = (bb, dist)

    @abc.abstractmethod
    def _get_neighbors(self, seq: str) -> dict[str, int]:
        """
        Get all neighbors of a sequence within a distance from it

        Parameters
        ----------
        seq: str
            sequence to get neighbors for

        Returns
        -------
        dict[str, int]
            a dictionary mapping each neighbor to its distance from the given sequence
            includes the given sequence with distance 0
        """

    def correct_sequence(self, observed_seq: str) -> str | None:
        """
        Given an observed observed_seq correct it to the correct to the expected observed_seq

        Parameters
        ----------
        observed_seq: str
            the observed DNA observed_seq

        Returns
        -------
        str | None
            the corrected DNA observed_seq; None if fails to correct
        """
        if observed_seq in self._hash_map:
            _existing_bb = self._hash_map[observed_seq][0]
            if _existing_bb is not None:
                return _existing_bb.tag
        return None


class LevenshteinDistHashMap(HashMapErrorCorrector):
    """
    A HashMapErrorCorrector that uses Levenshtein distance to build the hash map

    It is important to note when using this that a distance cutoff of 2 or greater
    can be very expensive in both time and memory. For example, a DNA tag with 11
    bp has 873,900 levenshtein neighbors with a distance of 3 or lower,
    and 8,853 with 2 or lower.

    Also, note that if the DNA tags in your building block set do not have a guaranteed
    levenshein distance of 2*distance_cutoff + 1 or more, then you can have collisions.
    In this case the hash map will fail to build, as there is no way to know which
    Building block should be mapped to that DNA sequence. You can disable this behavior
    with `asymmetrical=True`. In this case, the hash map with return either the best match
    or None if there are more than one possible match with the same distance
    """

    def _get_neighbors(self, seq: str) -> dict[str, int]:
        """
        Get all neighbors of a sequence within a given Levenshtein distance from it

        Recall the Levenshtein distance is the minimum number of single-character edits
        (inserts, deletions or substations) required to change one sequence into the other

        Parameters
        ----------
        seq: str
            sequence to get neighbors for

        Returns
        -------
        dict[str, int]
            a dictionary mapping each neighbor to its distance from the given sequence
            includes the given sequence with distance 0
        """
        neighbors = defaultdict(int, {seq: 0})

        current_neighbors = [seq]

        for d in range(1, self.distance_cutoff + 1):
            next_neighbors = []
            for current_seq in current_neighbors:
                n = len(current_seq)

                # Substitutions
                for i in range(n):
                    for nuc in "AGCT":
                        if current_seq[i] != nuc:
                            neighbor = current_seq[:i] + nuc + current_seq[i + 1 :]
                            if neighbor not in neighbors:
                                neighbors[neighbor] = d
                                next_neighbors.append(neighbor)

                # Insertions
                for i in range(n + 1):
                    for nuc in "AGCT":
                        neighbor = current_seq[:i] + nuc + current_seq[i:]
                        if neighbor not in neighbors:
                            neighbors[neighbor] = d
                            next_neighbors.append(neighbor)

                # Deletions
                for i in range(n):
                    neighbor = current_seq[:i] + current_seq[i + 1 :]
                    if neighbor not in neighbors:
                        neighbors[neighbor] = d
                        next_neighbors.append(neighbor)

            current_neighbors = next_neighbors

        return neighbors


class HammingDistHashMap(HashMapErrorCorrector):
    """
    A HashMapErrorCorrector that uses Hamming distance to build the hash map

    Note that if the DNA tags in your building block set do not have a guaranteed
    levenshein distance of 2*distance_cutoff + 1 or more, then you can have collisions.
    In this case, the hash map will fail to build, as there is no way to know which
    Building block should be mapped to that DNA sequence. You can disable this behavior
    with `asymmetrical=True`. In this case, the hash map with return either the best match
    or None if there is more than one possible match with the same distance
    """

    def _get_neighbors(self, seq: str) -> dict[str, int]:
        """
        Get all neighbors of a sequence within a given Hamming distance from it

        Recall the Hamming distance is the minimum number of single-character
        substations required to change one sequence into the other

        Parameters
        ----------
        seq: str
            sequence to get neighbors for

        Returns
        -------
        dict[str, int]
            a dictionary mapping each neighbor to its distance from the given sequence
            includes the given sequence with distance 0
        """
        neighbors = defaultdict(int, {seq: 0})

        current_neighbors = [seq]

        for d in range(1, self.distance_cutoff + 1):
            next_neighbors = []
            for current_seq in current_neighbors:
                n = len(current_seq)

                # Substitutions
                for i in range(n):
                    for nuc in "AGCT":
                        if current_seq[i] != nuc:
                            neighbor = current_seq[:i] + nuc + current_seq[i + 1 :]
                            if neighbor not in neighbors:
                                neighbors[neighbor] = d
                                next_neighbors.append(neighbor)

            current_neighbors = next_neighbors

        return neighbors


class BuildingBlockSetTagCaller:
    """
    Calls building blocks given a tag query and building block set

    BuildingBlockCallers can attempt to correct errors in the tag query
    based on knowledge they have about the building block set and how
    it was generated (if that info is given). See the docs on
    "Error Correction" for more details on how this works.
    """

    def __init__(
        self,
        building_block_tag_section: BuildingBlockBarcodeSection,
        building_block_set: TaggedBuildingBlockSet,
        disable_error_correction: bool = False,
    ):
        """
        Initialize a BuildingBlockCaller

        Parameters
        ----------
        building_block_tag_section: BuildingBlockBarcodeSection
            the building block tag section from barcode to decode from
        building_block_set: BuildingBlockSet
            the building block set to look for matches in
        disable_error_correction: bool, default=False
            even if an error correction method is provided, do
            not attempt error correction of DNA tags
        """
        self.building_block_tag_section = building_block_tag_section
        self.building_block_set = building_block_set
        self._validate_error_corrector(disable_error_correction)

    def _validate_error_corrector(self, disable_error_correction: bool):
        """
        Validate that the error corrector is set and has a valid method

        There are currently 3 error correction strings that can be parsed:
        - "hamming_dist:<distance>" for Hamming distance based error correction
        - "levenshtein_dist:<distance>" for Levenshtein distance based error correction
        - "hamming_code:<name>" for a QuaternaryHammingDecoder that has been loaded
        from the config file

        "hamming_dist_<distance>" and "levenshtein_dist_<distance>" will create hashmaps
        for rapid error correction, while "<name>" will load a pre-defined
         QuaternaryHammingDecoder.

        Note distance for the hamming_dist and levenshtein_dist modes refers to
        the maximum distance
        for neighbors. So a distance of 1 means that the error corrector will
        look for matches among
        all hamming/levenshtein neighbors that are 1 away from the given tag.

        "hamming_dist_<distance>" and "levenshtein_dist_<distance>" also have an
         asymmetric mode that
        can be triggered by adding "asymmetric" after the distance separated by a ',':
        e.g. "hamming_dist:1,asymmetric".
        """
        self.error_corrector: ErrorCorrector | None = None

        # if disable_error_correction is True, do not set an error corrector
        if disable_error_correction:
            return
        # if no error correction mode is set, do not set an error corrector
        elif self.building_block_tag_section.error_correction_mode is None:
            return
        # else parse the error corrector
        else:
            # parse out the type and info from the error correction mode
            try:
                _correction_type, _correction_info = (
                    self.building_block_tag_section.error_correction_mode.split(":")
                )
            except ValueError as e:
                if ":" not in self.building_block_tag_section.error_correction_mode:
                    raise ErrorCorrectorException(
                        f"error correction mode "
                        f"{self.building_block_tag_section.error_correction_mode} "
                        f"for barcode section {self.building_block_tag_section.section_name} "
                        f"does not have a valid format; missing ':', should be '<type>:<info>'"
                    ) from e
                else:
                    raise ErrorCorrectorException(
                        f"error correction mode "
                        f"{self.building_block_tag_section.error_correction_mode} "
                        f"for barcode section {self.building_block_tag_section.section_name} "
                        f"has too many reserved ':' characters,  should be '<type>:<info>'"
                    ) from e
            if _correction_type in ["hamming_dist", "levenshtein_dist"]:
                # collect info
                _splits = _correction_info.split(",")
                _dist_str = _splits[0]
                if len(_splits) == 1:
                    _asymmetrical = False
                elif len(_splits) == 2:
                    _asymmetrical = True
                else:
                    raise ErrorCorrectorException(
                        f"error correction mode "
                        f"{self.building_block_tag_section.error_correction_mode} "
                        f"for barcode section {self.building_block_tag_section.section_name} "
                        f"has unrecognized information: {_splits[2:]}"
                    )
                # check for valid integer
                try:
                    _dist = int(_splits[0])
                except ValueError as e:
                    raise ErrorCorrectorException(
                        f"error correction mode "
                        f"{self.building_block_tag_section.error_correction_mode} "
                        f"for barcode section {self.building_block_tag_section.section_name} "
                        f"does not have a valid <distance> for {_correction_info}; "
                        f"should be an integer"
                    ) from e

                if _correction_type == "hamming_dist":
                    self.error_corrector = HammingDistHashMap(
                        distance_cutoff=_dist,
                        bb_set=self.building_block_set,
                        asymmetrical=_asymmetrical,
                    )
                elif _correction_type == "levenshtein_dist":
                    self.error_corrector = LevenshteinDistHashMap(
                        distance_cutoff=_dist,
                        bb_set=self.building_block_set,
                        asymmetrical=_asymmetrical,
                    )
                else:
                    raise RuntimeError("UNREACHABLE; Raise issue if observed at runtime")
            elif _correction_type == "hamming_code":
                try:
                    self.error_corrector = QuaternaryHammingDecoder.load(name=_correction_info)
                except DELiConfigError as e:
                    raise ErrorCorrectorException(
                        f"cannot find an error corrector that matched the mode "
                        f"{self.building_block_tag_section.error_correction_mode} "
                        f"for barcode section "
                        f"{self.building_block_tag_section.section_name}"
                    ) from e
            else:
                raise ErrorCorrectorException(
                    f"error correction mode "
                    f"{self.building_block_tag_section.error_correction_mode} "
                    f"for barcode section {self.building_block_tag_section.section_name} "
                    f"has unrecognized error correction type '{_correction_type}'; "
                    f"see docs for valid error correction types"
                )

    def call_building_block(self, tag: str) -> BuildingBlockCall:
        """
        Call a building block given a nucleotide tag query

        Parameters
        ----------
        tag: str
            the nucleotide tag query to search

        Returns
        -------
        BuildingBlockCall
            the building block call
            will be a ValidBuildingBlockCall if match is found
            will be a FailedBuildingBlockCall if match is not found
        """
        if self.error_corrector is not None:
            _tag = self.error_corrector.correct_sequence(tag)
            if _tag is None:
                return FailedBuildingBlockCall()
            else:
                tag = _tag

        _call = self.building_block_set.search_tags(tag, fail_on_missing=False)

        if _call is None:
            return FailedBuildingBlockCall()
        else:
            # score is always 0 for building block matches
            return ValidBuildingBlockCall(_call, 0)
