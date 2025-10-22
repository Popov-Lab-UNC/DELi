"""Code for handling loading/validating DELi configs for decoding"""

import abc
import configparser
import functools
import inspect
import os
import shutil
import warnings
from pathlib import Path
from typing import Any, Callable, Literal, Optional, ParamSpec, TypeVar, Union
from typing_extensions import Self


P = ParamSpec("P")
R = TypeVar("R")

# DEFAULT SETTINGS
_BB_MASK_TOKEN_DEFAULT = "###"
_NUC_2_INT_DEFAULT = "A:0,T:1,C:2,G:3"

# DEFAULT HAMMING CODE ORDERS
_hamming_order_8_4 = "p0,p1,p2,d3,p4,d5,d6,d7"
_custom_order_8_4 = "p0,p1,p2,d3,p4,d5,d6,d7"

_hamming_order_16_5 = "p0,p1,p2,d3,p4,d5,d6,d7,p8,d9,d10,d11,d12,d13,d14,d15"
_custom_order_16_5 = "p0,p1,p2,d3,p4,d5,d6,d7,p8,d9,d10,d11,d12,d13,d14,d15"

_hamming_order_7_3 = "p1,p2,d3,p4,d5,d6,d7"
_custom_order_7_3 = "p1,p2,d3,p4,d5,d6,d7"

_hamming_order_15_4 = "p1,p2,d3,p4,d5,d6,d7,p8,d9,d10,d11,d12,d13,d14,d15"
_custom_order_15_4 = "p1,p2,d3,p4,d5,d6,d7,p8,d9,d10,d11,d12,d13,d14,d15"

DELI_DATA_SUB_DIRS = ["libraries", "building_blocks"]

DELI_CONFIG = None


class DELiConfigError(Exception):
    """raised when a DELi config is invalid or missing"""

    pass


def set_deli_data_dir(data_dir: Union[str, Path, None]) -> None:
    """Sets the deli data directory path"""
    get_deli_config().deli_data_dir = Path(data_dir) if isinstance(data_dir, str) else data_dir


def get_deli_config():
    """Get the DELi config, loading it lazily if not already loaded"""
    global DELI_CONFIG
    if DELI_CONFIG is None:
        _deli_config_dir = os.environ.get("DELI_CONFIG", None)
        if (_deli_config_dir is not None) and (_deli_config_dir != ""):
            DELI_CONFIG = _DeliConfig.load_config(_deli_config_dir)
        elif (Path.home() / ".deli").exists():
            DELI_CONFIG = _DeliConfig.load_config(Path.home() / ".deli")
        else:
            warnings.warn(
                f"no DELi config file in home directory; "
                f"creating default DELi config: {Path.home() / '.deli'}",
                stacklevel=1,
            )
            init_deli_config(Path.home() / ".deli", fail_on_exist=False)
            DELI_CONFIG = _DeliConfig.load_config(Path.home() / ".deli")

        # take the data_dir from the environment variable if it exists
        _deli_data_dir = os.environ.get("DELI_DATA_DIR", None)
        if (isinstance(_deli_data_dir, str)) and (_deli_data_dir != ""):
            _deli_data_dir_path = Path(_deli_data_dir).resolve()
            DELI_CONFIG.deli_data_dir = _deli_data_dir_path

    return DELI_CONFIG


def load_deli_config(path: Union[str, Path]) -> None:
    """
    Load a new deli config from a given path

    Overrides any existing DELi config settings
    """
    global DELI_CONFIG
    DELI_CONFIG = _DeliConfig.load_config(path)


class DeliDataNotFound(Exception):
    """raised when a file cannot be found in DELi data directory"""

    pass


class DeliDataDirError(Exception):
    """raised when a DELi data directory cannot be found or is invalid (missing sub-dirs)"""

    pass


class _DeliConfig:
    """
    Struct to hold info on DELi settings

    Parameters
    ----------
    deli_data_dir: Path | None
        where the DELi data directory is located
    bb_mask: str
        the 3-char long token to replace any masked building blocks
        masked building blocks result from synthon-based analysis
    nuc_2_int: dict[str, int]
        the nucleotide-to-integer mapping
    """

    def __init__(
        self,
        bb_mask: str,
        nuc_2_int: dict[str, int],
        deli_data_dir: Optional[Path] = None,
        hamming_codes: Optional[dict[str, tuple[list[int], list[int]]]] = None,
    ):
        # process deli data director
        self._deli_data_dir: Path | None = None
        if deli_data_dir is not None:
            self.deli_data_dir = deli_data_dir
        self._bb_mask: str = bb_mask
        self._nuc_2_int: dict[str, int] = nuc_2_int

        self.hamming_codes: dict[str, tuple[list[int], list[int]]]
        if hamming_codes is not None:
            self.hamming_codes = hamming_codes
        else:
            self.hamming_codes = {}

    @property
    def deli_data_dir(self) -> Path:
        if self._deli_data_dir is None:
            raise DeliDataDirError("DELi data directory is not set")
        else:
            return self._deli_data_dir

    @deli_data_dir.setter
    def deli_data_dir(self, value) -> None:
        if value is None:
            self._deli_data_dir = None
        else:
            if not isinstance(value, Path):
                try:
                    _path = Path(value).resolve()
                except TypeError as e:
                    raise TypeError(
                        f"'deli_data_directory' must be a str, bytes or "
                        f"os.PathLike object, not {type(value)}"
                    ) from e
            else:
                _path = value.resolve()

            # validate the deli data directory
            validate_deli_data_dir(_path)
            self._deli_data_dir = _path

    @deli_data_dir.deleter
    def deli_data_dir(self) -> None:
        self._deli_data_dir = None

    @property
    def bb_mask(self) -> str:
        return self._bb_mask

    @bb_mask.setter
    def bb_mask(self, value) -> None:
        if len(value) == 3 and isinstance(value, str):
            self._bb_mask = value
        else:
            if len(value) != 3:
                raise DELiConfigError(
                    f"'bb_mask' must be a string of 3 characters, not {len(value)}: '{value}'"
                )
            if not isinstance(value, str):
                raise DELiConfigError(f"'bb_mask' must be a 'string', not '{type(value)}'")

    @bb_mask.deleter
    def bb_mask(self) -> None:
        raise RuntimeError("Cannot delete 'bb_mask'")

    @property
    def nuc_2_int(self) -> dict[str, int]:
        return self._nuc_2_int

    @nuc_2_int.setter
    def nuc_2_int(self, value):
        if not isinstance(value, dict):
            try:
                _dict = dict(value)
            except TypeError as e:
                raise TypeError(
                    f"'nuc_2_int' must be a type 'dict' or castable "
                    f"to a dict, found type '{type(value)}'"
                ) from e
        else:
            _dict = value

        # check that the keys are the nucleotides A, T, G, C
        if {"A", "T", "G", "C"} != set(_dict.keys()):
            raise DELiConfigError(
                f"'nuc_2_int' must contain the only the nucleotides 'A', 'T', 'G', 'C'; "
                f"found nucleotides '{set(_dict.keys())}'"
            )

        # check that all values are cast/castable to int
        for key, val in _dict.items():
            try:
                _dict[key] = int(val)
            except ValueError as e:
                raise ValueError(
                    f"nucleotide '{key}' in 'nuc_2_int' must be mapped to an integer, "
                    f"found value '{val}' of type '{type(val)}'"
                ) from e

        # check that only 0, 1, 2, 3 are used as values
        if {0, 1, 2, 3} != set(_dict.values()):
            raise DELiConfigError(
                f"'nuc_2_int' must map nucleotides to values 0, 1, 2, 3; "
                f"found values '{set(_dict.values())}'"
            )
        self._nuc_2_int = _dict

    @nuc_2_int.deleter
    def nuc_2_int(self):
        raise RuntimeError("Cannot delete 'nuc_2_int'")

    @classmethod
    def load_config(cls, path: Union[str, Path]) -> Self:
        """Helper func to load in config data"""
        config = configparser.RawConfigParser()
        try:
            if config.read(os.path.normpath(path)) is None:
                raise FileNotFoundError(
                    f"cannot find deli config file at '{os.path.normpath(path)}'"
                )
        except configparser.Error as e:
            raise DELiConfigError(
                f"error reading config file at '{os.path.normpath(path)}';\n"
                f"you can use `deli config init --overwrite` to generate a new valid "
                f"config file with default settings"
            ) from e

        try:
            _bb_mask = config.get("deli.buildingblocks", "BB_MASK")
        except configparser.NoOptionError as e:
            raise DELiConfigError(
                f"missing 'BB_MASK' option in 'deli.buildingblocks' "
                f"config file at '{os.path.normpath(path)}';\n"
                f"you can use `deli config init --overwrite` to generate a new valid "
                f"config file with default settings"
            ) from e
        except configparser.NoSectionError as e:
            raise DELiConfigError(
                f"missing 'deli.buildingblocks' section in "
                f"config file at '{os.path.normpath(path)}';\n"
                f"you can use `deli config init --overwrite` to generate a new valid "
                f"config file with default settings"
            ) from e

        try:
            _nuc_2_int = {
                pair.split(":")[0].strip(): int(pair.split(":")[1].strip())
                for pair in config.get("deli.hamming", "nuc_2_int").strip().split(",")
            }
        except configparser.NoOptionError as e:
            raise DELiConfigError(
                f"missing 'nuc_2_int' option in 'deli.hamming' "
                f"config file at '{os.path.normpath(path)}';\n"
                f"you can use `deli config init --overwrite` to generate a new valid "
                f"config file with default settings"
            ) from e
        except configparser.NoSectionError as e:
            raise DELiConfigError(
                f"missing 'deli.hamming' section in config file at '{os.path.normpath(path)}';\n"
                f"you can use `deli config init --overwrite` to generate a new valid "
                f"config file with default settings"
            ) from e

        try:
            _deli_data_dir = config.get("deli.data", "deli_data_dir")
            if _deli_data_dir == "":
                _deli_data_dir_path = None
            else:
                _deli_data_dir_path = Path(config.get("deli.data", "deli_data_dir"))
        except (configparser.NoOptionError, configparser.NoSectionError):
            _deli_data_dir_path = None

        # extract all hamming codes from the config
        _codes: dict[str, tuple[list[int], list[int]]] = {}
        for section in config.sections():
            if section.startswith("deli.hamming."):
                name = section.split(".")[-1]
                try:
                    hamming_order = config.get(section, "hamming_order")
                except configparser.NoOptionError as e:
                    raise DELiConfigError(
                        f"missing 'hamming_order' option in hamming section '{section}'"
                    ) from e
                try:
                    custom_order = config.get(section, "custom_order")
                except configparser.NoOptionError as e:
                    raise DELiConfigError(
                        f"missing 'custom_order' option in hamming section '{section}'"
                    ) from e

                try:
                    _true_order_nums = [int(_[1:]) for _ in hamming_order.split(",")]
                except ValueError as e:
                    raise DELiConfigError(
                        f"invalid hamming order tokens in '{hamming_order}' of section "
                        f"'{section}'; see DELi hamming docs for details about how to "
                        f"define hamming code order"
                    ) from e
                if (_true_order_nums != sorted(_true_order_nums)) or (
                    _true_order_nums[0] not in {0, 1}
                ):
                    raise DELiConfigError(
                        f"invalid hamming order '{hamming_order}' in section '{section}'; "
                        f"see DELi hamming docs for details about how to define "
                        f"hamming code order"
                    )

                try:
                    _real_order_nums = [int(_[1:]) for _ in custom_order.split(",")]
                except ValueError as e:
                    raise DELiConfigError(
                        f"invalid custom order tokens in '{custom_order}' of section "
                        f"'{section}'; see DELi hamming docs for details about how to "
                        f"define custom order"
                    ) from e
                if len(_real_order_nums) != len(_true_order_nums):
                    raise DELiConfigError(
                        f"hamming section '{section}' has a parity "
                        f"length mismatch; see DELi hamming docs for details"
                    )
                _codes[name] = (_true_order_nums, _real_order_nums)

        return cls(
            deli_data_dir=_deli_data_dir_path,
            bb_mask=_bb_mask,
            nuc_2_int=_nuc_2_int,
            hamming_codes=_codes,
        )

    def get_hamming_code(self, name: str) -> tuple[list[int], list[int]]:
        """
        Get the hamming code order and custom order for a given name

        Parameters
        ----------
        name: str
            name of the hamming code to load, e.g. "8_4"

        Returns
        -------
        tuple[list[int], list[int]]
            (hamming_order ints, custom_order ints)
        """
        try:
            return self.hamming_codes[name]
        except KeyError as e:
            raise DELiConfigError(
                f"unrecognized hamming code '{name}'; "
                f"available codes are {list(self.hamming_codes.keys())}; "
                f"new codes can be added to the deli config file; see 'Hamming' docs for details"
            ) from e


def validate_deli_data_dir(deli_data_dir_: Path) -> bool:
    """
    Validate that the given DELi data directory exists and has all the required sub-directories

    Parameters
    ----------
    deli_data_dir_: Path
        path to the DELi data directory to validate

    Returns
    -------
    bool
        True if the directory is valid, raises DeliDataDirError otherwise
    """
    if not deli_data_dir_.exists():
        raise FileNotFoundError(
            f"directory at '{deli_data_dir_}' does not exist; "
            f"you can create a new DELi data directory using 'deli data init {deli_data_dir_}'"
        )

    if not deli_data_dir_.is_dir():
        raise NotADirectoryError(f"'{deli_data_dir_}' is not a directory")

    sub_dirs = {p.stem for p in deli_data_dir_.iterdir()}
    missing_sub_dirs = set(DELI_DATA_SUB_DIRS) - sub_dirs
    if len(missing_sub_dirs) > 0:
        raise DeliDataDirError(
            f"DELi data directory '{deli_data_dir_}' is invalid; "
            f"missing sub-directories {list(missing_sub_dirs)}"
            f"use 'deli data init --fix-missing {deli_data_dir_}' to add missing sub-directories"
        )
    return True


def init_deli_data_directory(path: Path, fail_on_exist: bool = True, overwrite: bool = False):
    """
    Create a new DELi Data Directory with all the necessary sub folders

    Notes
    -----
    If `fail_on_exist` is True, will raise an exception if the directory already exists
    even if `overwrite` is set to True.

    Parameters
    ----------
    path: Path
        path to create deli data dir at
    fail_on_exist: bool, default = True
        if True, will raise an exception if the directory path already exists
        if False, will try and create the sub-dirs in the existing directory
    overwrite: bool, default = False
        if a (sub) directory already exists, delete it and create a new directory
    """
    path = path.resolve()
    if path.exists():
        if fail_on_exist:
            raise FileExistsError(
                f"'{path}' already exists; set 'fail_on_exist' to False to overwrite"
            )
        elif not path.is_dir():
            raise NotADirectoryError(f"'{path}' is not a directory")
    else:
        os.makedirs(path)

    for sub_dir in DELI_DATA_SUB_DIRS:
        _sub_path = path / sub_dir
        if _sub_path.exists() and overwrite:
            shutil.rmtree(_sub_path)
        os.makedirs(path / sub_dir, exist_ok=True)


def init_deli_config(
    path: Optional[Path] = None,
    fail_on_exist: bool = True,
):
    """
    Create a default Deli config directory

    Parameters
    ----------
    path: Optional[Union[str, os.PathLike]]
        path to create deli config dir at
        will be at $USER/.deli if left as None
    fail_on_exist: bool, default = True
        if True, will raise an exception if the directory path already exists
        if False, will try and create the sub-dirs in the existing directory
    """
    if path is None:
        _path = Path.home() / ".deli"
    else:
        _path = path

    if _path.exists() and fail_on_exist:
        raise FileExistsError(
            f"'{_path}' already exists; "
            f"set 'fail_on_exist' to False to overwrite OR use `deli config init --overwrite`"
        )

    _config = (
        f"[DEFAULT]\n\n"
        f"[deli.data]\n"
        f"deli_data_dir = \n\n"
        f"[deli.hamming]\n"
        f"nuc_2_int = {_NUC_2_INT_DEFAULT}\n\n"
        f"[deli.hamming.8_4]\n"
        f"hamming_order = {_hamming_order_8_4}\n"
        f"custom_order = {_custom_order_8_4}\n\n"
        f"[deli.hamming.16_5]\n"
        f"hamming_order = {_hamming_order_16_5}\n"
        f"custom_order = {_custom_order_16_5}\n\n"
        f"[deli.hamming.7_3]\n"
        f"hamming_order = {_hamming_order_7_3}\n"
        f"custom_order = {_custom_order_7_3}\n\n"
        f"[deli.hamming.15_4]\n"
        f"hamming_order = {_hamming_order_15_4}\n"
        f"custom_order = {_custom_order_15_4}\n\n"
        f"[deli.buildingblocks]\n"
        f"BB_MASK = {_BB_MASK_TOKEN_DEFAULT}\n"
    )
    with open(_path, "w") as f:
        f.write(_config)


def _build_default_hamming_code_strings(
    include_extra_parity: bool = False,
) -> dict[str, tuple[str, str]]:
    """
    Build the hamming code and custom order strings for various length tags

    Only builds codes of length 7 to 16

    Parameters
    ----------
    include_extra_parity: bool, default = False
        use an additional parity bit in the hamming code

    Returns
    -------
    dict[str, tuple[str, str]]
        a dictionary mapping hamming code names to tuples of (hamming_order, custom_order)
    """
    from math import log2

    _codes = {"7_3": ("p1,p2,d3,p4,d5,d6,d7", "p1,p2,d3,p4,d5,d6,d7")}
    if include_extra_parity:
        _codes["8_4"] = ("p0,p1,p2,d3,p4,d5,d6,d7", "p0,p1,p2,d3,p4,d5,d6,d7")

    for i in range(9, 16):
        _hamming_code_str = ",".join(
            [f"p{j}" if log2(j).is_integer() else f"d{j}" for j in range(1, i + 1)]
        )
        _codes[f"{i}_4"] = (_hamming_code_str, _hamming_code_str)
        if include_extra_parity:
            _codes[f"{i + 1}_5"] = ("p0," + _hamming_code_str, "p0," + _hamming_code_str)
    return _codes


def accept_deli_data_name(
    sub_dir: Literal["building_blocks", "libraries", "indexes", "barcodes", "hamming"],
    extension: str,
) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """
    Decorator to allow load functions to take name of DELi data object

    Notes
    -----
    This should be used decorate the `load` function of any
    class that has the `DeliDataLoadable` mixin

    It will result in the function taking in, rather than a full path,
    the name of the object (e.g. "MyMegaLibrary") and the
    DELi config object that defines where deli_data_dir is

    Set the sub-dir (libraries, indexes, building_blocks, barcodes)
    when assigning the decorator based on where the data should be
    See `Storing DEL info` in docs for more details

    Also set the file `extension` so that the write file loaders
    can be used and file type checked

    Not all DEL objects use the same file format.
    DocStrings of a given object will define which
    file formats can be used

    Parameters
    ----------
    sub_dir: str
        the name of the DELi data sub directory to look in
    extension: str
        the file extension for matching files

    Returns
    -------
    decorated_function: Callable[P, R]
    """

    def _build_deli_data_path(path: Union[str, os.PathLike]) -> Path:
        if os.path.exists(path):
            return Path(path)

        _sub_dir_path = get_deli_config().deli_data_dir / sub_dir
        if not os.path.exists(_sub_dir_path):
            raise DeliDataNotFound(
                f"cannot find DELi data subdirectory at `{_sub_dir_path}`; "
                f"did you check that deli_data_dir is set correctly?"
            )

        _file_name = os.path.basename(path).split(".")[0] + "." + extension
        file_path = get_deli_config().deli_data_dir / sub_dir / _file_name

        if not os.path.exists(file_path):
            raise DeliDataNotFound(f"cannot find file '{path}.{extension}' in {_sub_dir_path}")

        return file_path

    try:
        decorator = _build_argument_validation_decorator(_build_deli_data_path, "path")
    except ValueError as err:  # will throw this error if the function lacks a "path" argument
        raise RuntimeError(
            "cannot decorate function without 'path' parameter "
            "with the `accept_deli_data_name` decorator"
        ) from err
    return decorator


def _build_argument_validation_decorator(
    validator_func: Callable[[Any], Any], target_arg_name: str
):
    """
    General decorator generator for validating a given argument with a validator function.
    """

    def outer(func: Callable[P, R]) -> Callable[P, R]:
        @functools.wraps(func)
        def inner(*args: P.args, **kwargs: P.kwargs) -> R:
            pos = None
            if target_arg_name in kwargs:
                target_arg_value = kwargs[target_arg_name]
            else:
                signature = inspect.signature(func)
                pos = list(signature.parameters.keys()).index(target_arg_name)
                target_arg_value = args[pos]
            new_arg = validator_func(target_arg_value)
            if pos:
                _new_args = (*args[:pos], new_arg, *args[pos + 1 :])
            else:
                kwargs[target_arg_name] = new_arg
                _new_args = (*args,)
            return func(*_new_args, **kwargs)

        return inner

    return outer


class DeliDataLoadable(abc.ABC):
    """Mixin for objects that can be loaded from DeliDataDir"""

    loaded_from: Union[os.PathLike, str] = ""

    @classmethod
    @abc.abstractmethod
    @accept_deli_data_name(sub_dir="barcodes", extension="json")
    def load(cls, name_or_path: str):
        """Load the file into the object"""
        raise NotImplementedError()
