"""handles loading a parsing of DELi run settings"""

from dataclasses import dataclass
from typing import Any


class DELiSettingError(Exception):
    """raised when DELi has an issue with reading/loading settings"""

    pass


def _check_bool(var: Any, var_name: str) -> bool:
    """If a variable is a bool, check that it can be loaded as such"""
    if isinstance(var, bool):
        return var
    if isinstance(var, str):
        if var.upper() == "TRUE":
            return True
        elif var.upper() == "FALSE":
            return False
        elif var == "":
            raise DELiSettingError(f"missing value for setting '{var_name}'")
        else:
            raise DELiSettingError(f"invalid boolean value '{var}' for setting '{var_name}'")
    if isinstance(var, int):
        return bool(var)
    raise DELiSettingError(f"unrecognized value '{var}' for setting '{var_name}'")


def _check_int(var: Any, var_name: str) -> int:
    """If a variable is a integer, check that it can be loaded as such"""
    if isinstance(var, bool):
        return int(var)
    if isinstance(var, str):
        try:
            return int(var)
        except ValueError as e:
            raise DELiSettingError(
                f"invalid integer value '{var}' " f"for setting '{var_name}'"
            ) from e
    if isinstance(var, int):
        return var
    raise DELiSettingError(f"unrecognized value '{var}' for setting '{var_name}'")


@dataclass
class MatcherSettings:
    """Holds setting for the `BarcodeMatcher` object"""

    error_tolerance: int = 3
    primer_match_length: int = -1
    rev_comp: bool = True

    def __post_init__(self):
        """Check types post init"""
        self.error_tolerance = _check_int(self.error_tolerance, "error_tolerance")
        self.primer_match_length = _check_int(self.primer_match_length, "primer_match_length")
        self.rev_comp = _check_bool(self.rev_comp, "rev_comp")


@dataclass
class CallerSettings:
    """Holds setting for the `BarcodeCaller` object"""

    hamming: bool = True

    def __post_init__(self):
        """Check types post init"""
        self.hamming = _check_bool(self.hamming, "hamming")


@dataclass
class CubeGenSettings:
    """Holds setting for the `CubeGenerator` object"""

    disynthon: bool = True
    monosynthon: bool = False
    split_by_index: bool = False
    split_by_lib: bool = True
    umi_cluster: bool = False
    normalize: bool = False

    def __post_init__(self):
        """Check types post init"""
        self.disynthon = _check_bool(self.disynthon, "disynthon")
        self.monosynthon = _check_bool(self.monosynthon, "monosynthon")
        self.split_by_index = _check_bool(self.split_by_index, "split_by_index")
        self.split_by_lib = _check_bool(self.split_by_lib, "split_by_lib")
        self.umi_cluster = _check_bool(self.umi_cluster, "umi_cluster")
        self.normalize = _check_bool(self.normalize, "normalize")
