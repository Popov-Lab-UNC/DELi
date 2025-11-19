"""Classes and functions for barcode calling with/without error correction"""

import abc
import math
from collections import defaultdict
from typing import Generic, TypeVar

from typing_extensions import deprecated

from deli._hamming import BaseQuaternaryHamming
from deli.configure import DELiConfigError, get_deli_config


T = TypeVar("T")


class BarcodeCallerError(Exception):
    """base class for all barcode caller exceptions"""

    pass


class ValidCall(Generic[T]):
    """
    Represents a valid call from a BarcodeCaller

    Notes
    -----
    Can be parameterized with the type of object being called

    Attributes
    ----------
    obj: T
        the object that was called
    score: float
        the score of the call; lower is better
    """

    def __init__(self, obj: T, score: float):
        self.obj = obj
        self.score = score


class FailedBarcodeLookup:
    """Returned when a barcode lookup fails in a BarcodeCaller"""

    def __init__(self, barcode: str):
        self.barcode = barcode


class AmbiguousBarcodeCall(FailedBarcodeLookup):
    """Returned when a barcode call is ambiguous in a BarcodeCaller"""

    def __init__(self, barcode: str):
        super().__init__(barcode=barcode)


class BarcodeCaller(abc.ABC, Generic[T]):
    """
    Base class for all BarcodeCallers

    Notes
    -----
    Can be parameterized with the type of object being called

    Parameters
    ----------
    tag_map: dict[str, T]
        a dictionary mapping each valid tag to its associated object
    """

    def __init__(self, tag_map: dict[str, T]):
        self._tag_map = tag_map

    @abc.abstractmethod
    def decode_barcode(self, observed_barcode: str) -> ValidCall[T] | FailedBarcodeLookup:
        """
        Given an observed barcode sequence, decode it to the associated object

        Notes
        -----
        This method returns an object of the class bound based on the
        tag_map given at initialization

        Parameters
        ----------
        observed_barcode: str
            the observed DNA barcode

        Returns
        -------
        ValidCall | FailedBarcodeLookup
            The called object or FailedBarcodeLookup if no valid call could be made
        """
        raise NotImplementedError()


class SingleItemBarcodeCaller(BarcodeCaller[T]):
    """
    A BarcodeCaller that always just return the single item it contains

    Notes
    -----
    This is for compatability when there is only one possible item to call
    but the decoder needs a BarcodeCaller object

    Is parameterized with the type of object being called

    Parameters
    ----------
    item: T
        the single item this BarcodeCaller will always return
    """

    def __init__(self, item: T):
        super().__init__({"X": item})  # empty tag map

    def decode_barcode(self, observed_barcode: str) -> ValidCall[T]:
        """
        Ignores the observed barcode and just return the single item it contains

        The score of a call from this caller is always 0 (perfect)

        Parameters
        ----------
        observed_barcode: str
            ignored; here for compatability

        Returns
        -------
        ValidCall
            the item used to initialize the SingleItemBarcodeCaller as a valid call
        """
        return ValidCall(self._tag_map["X"], 0)


class GenericBarcodeCaller(BarcodeCaller[T]):
    """
    A BarcodeCaller that has no error correction

    This simply looks up the observed barcode in the tag map and returns
    the associated object if it exists

    Parameters
    ----------
    tag_map: dict[str, T]
        a dictionary mapping each valid tag to its associated object
    """

    def __init__(self, tag_map: dict[str, T]):
        super().__init__(tag_map=tag_map)

    def decode_barcode(self, observed_barcode: str) -> ValidCall[T] | FailedBarcodeLookup:
        """
        Given an observed barcode return the associated object if it exists

        Parameters
        ----------
        observed_barcode: str
            the observed DNA barcode

        Returns
        -------
        ValidCall | FailedBarcodeLookup
            The object associated with the observed barcode;
            FailedBarcodeLookup if no exact match found
        """
        obj = self._tag_map.get(observed_barcode, None)
        if obj is not None:
            return ValidCall(obj, 0)
        else:
            return FailedBarcodeLookup(observed_barcode)


class HammingDecodeError(BarcodeCallerError):
    """Raised when the hamming decoder fails to decode an observed barcode"""

    pass


@deprecated(
    "QuaternaryHammingBarcodeCaller is deprecated and will be removed in a future release; "
    "use a HammingDistBarcodeCaller instead"
)
class QuaternaryHammingBarcodeCaller(BaseQuaternaryHamming, BarcodeCaller[T]):
    """
    Generates and decodes quaternary hamming codes

    Notes
    -----
    It is externally unlikely that this is the class you want to use, it
    requires that you built your DNA tags using a valid hamming code (not
    just creating a set of tags with a certain minimum hamming distance).

    Instead, consider using a HammingDistBarcodeCaller or
    LevenshteinDistBarcodeCaller for more general error correction.
    While more memory intensive, they perform the same function and are
    just as fast for decoding.

    Parameters
    ----------
    tag_map: dict[str, T]
        a dictionary mapping each valid tag to its associated object
    parity_map: list[int]
        a list defining where the correctly ordered hamming bits are.
        See the documentation for more details
    has_extra_parity: bool
        True if the sequences being decoded has extra parity.
        See the documentation for more details

    Attributes
    ----------
    real_size: int
        the real size of the hamming code (number of data bits)
    hamming_size: int
        the size of the hamming code (number of data + parity bits)
    """

    def __init__(self, tag_map: dict[str, T], parity_map: list[int], has_extra_parity: bool):
        super().__init__()

        self._tag_map = tag_map  # since MRO would call BarcodeCaller init
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
    def load(cls, name, tag_map: dict[str, T]) -> "QuaternaryHammingBarcodeCaller":
        """
        Load a hamming code into a decoder

        NOTE: to load, need a deli.hamming.<name> section
        in the ".deli" config file with both a hamming and custom order field.
        See docs for more details

        Parameters
        ----------
        name: str
            name of the hamming code to load, e.g. "8_4"
        tag_map: dict[str, Object]
            a dictionary mapping each valid tag to its associated object

        Returns
        -------
        QuaternaryHammingBarcodeCaller
        """
        _config = get_deli_config()
        _true_order, _custom_order = _config.get_hamming_code(name)

        if 0 in _true_order:
            _has_extra_parity = True
        else:
            _has_extra_parity = False
            _custom_order = [_ - 1 for _ in _true_order]

        return cls(tag_map=tag_map, parity_map=_custom_order, has_extra_parity=_has_extra_parity)

    def decode_barcode(self, observed_seq: str) -> ValidCall[T] | FailedBarcodeLookup:
        """
        Given an observed barcode correct it to the correct to the expected barcode

        Parameters
        ----------
        observed_seq: str
            the observed DNA barcode

        Returns
        -------
        ValidCall | FailedBarcodeLookup
            The object associated with the corrected seq;
            FailedBarcodeLookup if fails to correct
        """
        _tag = [self.nuc_2_int_mapping[char] for char in observed_seq]
        try:
            _decoded_tag = self._hamming_decode(_tag)
            return ValidCall(self._tag_map["".join([self.int_2_nuc_mapping[_] for _ in _decoded_tag])], 0)
        except HammingDecodeError:
            return FailedBarcodeLookup(observed_seq)

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

        sub_parity = [self._get_sub_parity(_padded_bases, order=p) for p in range(self._parity, 0, -1)][::-1]

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


class HashMapCollisionError(BarcodeCallerError):
    """Raised when a hash map error corrector fails to build due to collisions"""

    pass


class HashMapBarcodeCaller(BarcodeCaller[T], abc.ABC):
    """
    Base class for all hash map based error correctors

    Hash map decoders work by taking each DNA barcode for each object,
    creating all its neighbors within a given distance function
    and mapping them to the building block itself.
    This way, a query can be looked up and "corrected"
    to its original building block tag in constant time at
    the cost of a some memory and longer initialization time.

    The distance of the match is also returned as the score of the call.

    By default, collisions are not allowed when building the hash map.
    So if more than one object has the same "neighboring" tag within the
    given distance cutoff, the hash map will fail to build.
    If necessary (and often it is) an asymmetrical mode where collisions
    are ignored is supported. This results in the colliding tag being mapped
    to object with the shortest distance, *or* if both are the same distance
    it will map to a None

    Will also keep track of the distance of the query to the object, asigning
    that as the call score

    Notes
    -----
    Distance cutoff is the distance from the origin. If you had a Hamming3
    set (meaning every element has a hamming distance of 3 or greater from
    all other members), then the distance cutoff should be 1, since this
    code can correct up to 1 error. For a Hamming5 set it would be 2, 7
    would be 3 and so on.

    Parameters
    ----------
    tag_map: dict[str, T]
        a dictionary mapping each valid tag to its associated object
    distance_cutoff: int, default=1
        the maximum distance for neighbors to be included in the hash map
    asymmetrical: bool, default=False
        whether to allow asymmetrical collisions when building the hash map
    """

    def __init__(self, tag_map: dict[str, T], distance_cutoff: int = 1, asymmetrical: bool = False):
        super().__init__(tag_map=tag_map)
        self.distance_cutoff = distance_cutoff
        self.asymmetrical = asymmetrical
        self._hash_map: dict[str, tuple[T | None, int]] = {}

        for tag_, obj in self._tag_map.items():
            # loop through all tags to get neighbors, saving smallest dist if conflict since
            #  these two map to the same compound
            possible_tags = self._get_neighbors(tag_)
            for tag, dist in possible_tags.items():
                if tag in self._hash_map:
                    if not self.asymmetrical:
                        colliding_obj = self._hash_map[tag][0]
                        if colliding_obj is not None:
                            # if the colliding bb is the same as the current one not an issue
                            if colliding_obj == obj:
                                continue
                            else:
                                raise HashMapCollisionError(
                                    f"collision detected when building {self.__class__.__name__} "
                                    f"with distance {self.distance_cutoff} for tag {tag}"
                                )
                        else:
                            # this is for type hinting and mypy
                            raise RuntimeError("this should never happen, please report to the developers")
                    else:
                        if self._hash_map[tag][1] == dist:
                            if self._hash_map[tag][0] == obj:
                                continue
                            else:
                                self._hash_map[tag] = (None, dist)
                        elif self._hash_map[tag][1] > dist:
                            self._hash_map[tag] = (obj, dist)
                else:
                    self._hash_map[tag] = (obj, dist)

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

    def decode_barcode(self, observed_barcode: str) -> ValidCall[T] | FailedBarcodeLookup:
        """
        Given an observed barcode correct it to the correct to the expected barcode

        Parameters
        ----------
        observed_barcode: str
            the observed DNA barcode

        Returns
        -------
        ValidCall | FailedBarcodeLookup
            The object associated with observed barcode;
            FailedDecodeAttempt if fails to correct.
            Will be an AmbiguousBarcodeCall if a collusion is
            detected and asymmetrical mode is on, otherwise a
            FailedBarcodeLookup
        """
        obj = self._hash_map.get(observed_barcode, None)
        if obj is not None:
            if obj[0] is None:
                return AmbiguousBarcodeCall(observed_barcode)
            return ValidCall(obj[0], obj[1])
        else:
            return FailedBarcodeLookup(observed_barcode)


class LevenshteinDistBarcodeCaller(HashMapBarcodeCaller[T]):
    """
    A HashMapBarcodeCaller that uses Levenshtein distance to build the hash map

    Notes
    -----
    It is important to note when using this that a distance cutoff of 2 or greater
    can be very expensive in both time and memory. For example, a DNA tag with 11
    bp has 873,900 levenshtein neighbors with a distance of 3 or lower,
    and 8,853 with 2 or lower.
    """

    # TODO could this have some jit with numba?
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


class HammingDistBarcodeCaller(HashMapBarcodeCaller[T]):
    """
    A HashMapBarcodeCaller that uses Hamming distance to build the hash map

    Notes
    -----
    The Hamming distance is only defined for sequences of equal length.
    Thus, only SNPs can be corrected for, not INDELs
    """

    # TODO could this have some jit with numba?
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


def get_barcode_caller(tag_map: dict[str, T], error_correction_mode_str: str) -> BarcodeCaller[T]:
    """
    Parse the error correction mode string and return a BarcodeCaller with those settings

    Notes
    -----
    There are currently 4 error correction strings that can be parsed:
    - "hamming_dist:<distance>" for Hamming distance based error correction
    - "levenshtein_dist:<distance>" for Levenshtein distance based error correction
    - "hamming_code:<name>" for a QuaternaryHammingBarcodeCaller that has been loaded.
      Note: this is deprecated and will be removed in a future release
    - "disable" to disable error correction

    '<distance>' for the hamming_dist and levenshtein_dist modes refers to
    the maximum distance cutoff for neighbors.

    "hamming_dist_<distance>" and "levenshtein_dist_<distance>" also have an
    asymmetric mode that can be triggered by adding "asymmetric" after the
    distance separated by a ','. For example "hamming_dist:1,asymmetric".

    Parameters
    ----------
    tag_map: dict[str, Object]
        a dictionary mapping each valid tag to its associated object
    error_correction_mode_str: str
        the error correction mode string to parse

    Returns
    -------
    BarcodeCaller[Object]
        the appropriate BarcodeCaller for the given error correction mode string
    """
    # parse out the type and info from the error correction mode
    if len(tag_map) == 1:
        # if only one item, no need for error correction
        return SingleItemBarcodeCaller(list(tag_map.values())[0])

    try:
        _correction_type, _correction_info = error_correction_mode_str.split(":")
    except ValueError as e:
        if ":" not in error_correction_mode_str:
            raise BarcodeCallerError(
                f"error correction mode "
                f"{error_correction_mode_str} "
                f"does not have a valid format; "
                f"missing ':', should be '<type>:<info>'"
            ) from e
        else:
            raise BarcodeCallerError(
                f"error correction mode "
                f"{error_correction_mode_str} "
                f"has too many reserved ':' characters; "
                f" should be '<type>:<info>'"
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
            raise BarcodeCallerError(
                f"error correction mode {error_correction_mode_str} has unrecognized information: {_splits[2:]}"
            )
        # check for valid integer
        try:
            _dist = int(_splits[0])
        except ValueError as e:
            raise BarcodeCallerError(
                f"error correction mode "
                f"{error_correction_mode_str} "
                f"does not have a valid <distance> for {_correction_info}; "
                f"should be an integer"
            ) from e

        if _correction_type == "hamming_dist":
            return HammingDistBarcodeCaller(
                tag_map=tag_map,
                distance_cutoff=_dist,
                asymmetrical=_asymmetrical,
            )
        elif _correction_type == "levenshtein_dist":
            return LevenshteinDistBarcodeCaller(
                tag_map=tag_map,
                distance_cutoff=_dist,
                asymmetrical=_asymmetrical,
            )
        else:
            # for the type checker
            raise RuntimeError("UNREACHABLE; Raise issue if observed at runtime")

    elif _correction_type == "hamming_code":
        try:
            return QuaternaryHammingBarcodeCaller.load(name=_correction_info, tag_map=tag_map)
        except DELiConfigError as e:
            raise BarcodeCallerError(
                f"cannot find an QuaternaryHamming error corrector that matched the mode {error_correction_mode_str}"
            ) from e
    else:
        raise BarcodeCallerError(
            f"error correction mode "
            f"{error_correction_mode_str} "
            f"has unrecognized error correction type '{_correction_type}'; "
            f"see docs for valid error correction types"
        )
