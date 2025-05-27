"""Handles hamming correction of DEL barcode tags"""

import math
from typing import List, Optional, Self

from deli.configure import get_deli_config, DeliDataLoadable, accept_deli_data_name


class DecodeError(Exception):
    """raised when the decoder fails to decode a barcode"""

    pass


class BaseQuaternaryHamming:
    """Gives a object the functionality to do hamming parity checks"""

    def __init__(self):
        """Initialize a BaseQuaternaryHamming"""
        self.nuc_2_int_mapping = get_deli_config().nuc_2_int
        self.int_2_nuc_mapping = {val: key for key, val in self.nuc_2_int_mapping.items()}

    @staticmethod
    def _get_sub_parity(bases: list[int], order: int = 0) -> int:
        """Get the sub parity of the bases base on the parity bit order"""
        _parity = 0
        s = math.ceil(len(bases) / (2**order))
        for i in range((2**order)):
            if i % 2 == 1:
                _parity += sum(bases[s * i : s * (i + 1)])
        return (4 - _parity) % 4

    @staticmethod
    def _get_global_parity(bases: list[int]) -> int:
        """Get the global parity of the base list"""
        return (4 - sum(bases)) % 4

    @staticmethod
    def _correct_error(bases: list[int], parity_list: List[int]):
        """Correct a bit given the current parity status"""
        pos_corr = 0
        for i in range(len(parity_list)):
            pos_corr += (parity_list[-(i + 1)] > 0) * (2**i)
        bases[pos_corr] = (bases[pos_corr] - (4 - max(parity_list))) % 4


class QuaternaryHammingDecoder(BaseQuaternaryHamming, DeliDataLoadable):
    """Generates and decodes quaternary hamming codes"""

    def __init__(self, parity_map: List[int], has_extra_parity: bool):
        """
        Initialize a QuaternaryHammingDecoder

        Parameters
        ----------
        parity_map: List[int]
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

        self.nuc_2_int_mapping = get_deli_config().nuc_2_int
        self.int_2_nuc_mapping = {val: key for key, val in self.nuc_2_int_mapping.items()}

        self._require_sort = sorted(self.parity_map) != self.parity_map

    def __eq__(self, other) -> bool:
        """Return true if the encoder has the same total length"""
        if isinstance(other, self.__class__):
            return self.parity_map == other.parity_map
        else:
            return False

    @classmethod
    @accept_deli_data_name("hamming", "txt")
    def load(cls, path) -> Self:
        """
        Load a hamming code into a decoder

        Notes
        -----
        Can be the name of a hamming config in the
        DELI data dir subfolder 'hamming'

        Parameters
        ----------
        path: str
            path or name of file (if in DELI data dir)

        Returns
        -------
        Self
        """
        _true_order: Optional[str] = None
        _real_order: Optional[str] = None
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                elif line.startswith("hamming_order:"):
                    _true_order = line.split(":")[-1].strip()
                elif line.startswith("custom_order:"):
                    _real_order = line.split(":")[-1].strip()

        if _true_order is None:
            raise RuntimeError(
                f"hamming file '{path}' missing the 'hamming_order' field; "
                f"see DELi hamming docs for details"
            )

        if _real_order is None:
            raise RuntimeError(
                f"hamming file '{path}' missing the 'custom_order' field; "
                f"see DELi hamming docs for details"
            )

        _true_order_nums = [int(_[1:]) for _ in _true_order.split(",")]
        if _true_order_nums != sorted(_true_order_nums):
            raise RuntimeError(
                f"hamming file '{path}' has an improper "
                f"'true parity order'; see DELi hamming docs for details"
            )

        _real_order_nums = [int(_[1:]) for _ in _real_order.split(",")]
        if len(_real_order_nums) != len(_true_order_nums):
            raise RuntimeError(
                f"hamming file '{path}' has a parity "
                f"length mismatch; see DELi hamming docs for details"
            )

        if 0 in _true_order_nums:
            _has_extra_parity = True
        else:
            _has_extra_parity = False
            _real_order_nums = [_ - 1 for _ in _true_order_nums]

        return cls(parity_map=_real_order_nums, has_extra_parity=_has_extra_parity)

    def decode_sequence(self, sequence: str) -> Optional[str]:
        """
        Given a string sequence of nucleotides, correct the sequence and return it

        Notes
        -----
        If decoding fails, will just return the original sequence
        and a False bool. If it is decoded without issue will return
        the decoded/correct seq but also a True bool

        Parameters
        ----------
        sequence: str
            the sequence to decode/correct

        Returns
        -------
        Optional[str]
            decoded/correct sequence
            will be None if hamming decoding fails
        """
        _tag = [self.nuc_2_int_mapping[char] for char in sequence]
        try:
            _decoded_tag = self._hamming_decode(_tag)
            return "".join([self.int_2_nuc_mapping[_] for _ in _decoded_tag])
        except DecodeError:
            return None

    def _hamming_decode(self, bases: List[int]) -> List[int]:
        """
        Decode a tag and correct any correctable errors

        Parameters
        ----------
        bases: Union[List[Base], Tag]
            the tag to decode/correct

        Raises
        ------
        DecodeError
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
                raise DecodeError("2 or more errors detected")
            if (_parity_set_length == 2) and (0 not in sub_parity):
                raise DecodeError("2 or more errors detected")

            overall_parity = self._get_global_parity(_padded_bases)
            if (overall_parity == 0) and (len(set(sub_parity)) == 2):
                raise DecodeError("2 or more errors detected")

            if max(sub_parity) != 0:
                self._correct_error(_padded_bases, sub_parity)

            elif overall_parity != 0:
                # need to correct the parity bit in the case that this is the one that is misread
                _padded_bases[0] = (4 - sum(_padded_bases[1:])) % 4

            # final check for overall parity
            overall_parity = self._get_global_parity(_padded_bases)
            if overall_parity != 0:
                raise DecodeError("2 or more errors detected")

        # if not extra parity
        else:
            if max(sub_parity) != 0:
                self._correct_error(_padded_bases, sub_parity)

        return _padded_bases[not self.has_extra_parity : _actual_length]
