"""Base classes for handing hamming codes"""

import math

from deli.configure import get_deli_config


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
    def _correct_error(bases: list[int], parity_list: list[int]):
        """Correct a bit given the current parity status"""
        pos_corr = 0
        for i in range(len(parity_list)):
            pos_corr += (parity_list[-(i + 1)] > 0) * (2**i)
        bases[pos_corr] = (bases[pos_corr] - (4 - max(parity_list))) % 4
