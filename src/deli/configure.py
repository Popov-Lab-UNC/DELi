"""Code for handling loading/validating DELi configs for decoding"""

import abc
import configparser
import functools
import inspect
import math
import os
from pathlib import Path
from typing import Any, Callable, Literal, Optional, ParamSpec, Self, TypeVar, Union


P = ParamSpec("P")
R = TypeVar("R")

BB_MASK_TOKEN_DEFAULT: str = "###"
MAX_INDEX_RISK_DIST_THRESHOLD_DEFAULT: int = 3
MAX_LIBRARY_RISK_DIST_THRESHOLD_DEFAULT: int = 4
NUC_2_INT_DEFAULT: dict[str, int] = {"A": 0, "T": 1, "C": 2, "G": 3}


DELI_DATA_SUB_DIRS = ["hamming", "barcodes", "libraries", "indexes", "building_blocks"]


def set_deli_data_dir(data_dir: Union[str, Path]) -> None:
    """Sets the deli data directory path"""
    DELI_CONFIG.deli_data_dir = Path(data_dir) if isinstance(data_dir, str) else data_dir


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

    Attributes
    ----------
    deli_data_dir : str
        where the DELi data directory is located
        default location is `~/.deli/data`
        should follow format specific in 'Storing DEL info'
    bb_mask: str
        the 3 char long token to replace any masked building blocks
        masked building blocks result from synthon-based analysis
    max_index_risk_dist_threshold: int
        the maximum risk willing to make in an index call
        if risk is above this number will fail the call
    max_library_risk_dist_threshold: int
        the maximum risk willing to make in a library call
        if risk is above this number will fail the call
    nuc_2_int: dict[str, int]
        the nucleotide to integer mapping
    """

    def __init__(self, **kwargs):
        # process deli data director
        self.deli_data_dir: Path
        if kwargs.get("DELI_DATA_DIR", None) is None:
            self.deli_data_dir = Path(os.path.join(os.path.expanduser("~"), ".deli", "deli_data"))
        else:
            self.deli_data_dir = Path(os.fspath(os.path.expanduser(kwargs.get("DELI_DATA_DIR"))))
        _validate_deli_data_dir(self.deli_data_dir)

        self.bb_mask: str = (
            str(kwargs["BB_MASK"]) if kwargs.get("BB_MASK") is not None else BB_MASK_TOKEN_DEFAULT
        )
        self.max_index_risk_dist_threshold: int = (
            int(kwargs["MAX_INDEX_RISK_DIST_THRESHOLD"])
            if (kwargs.get("MAX_INDEX_RISK_DIST_THRESHOLD") is not None)
            else (MAX_INDEX_RISK_DIST_THRESHOLD_DEFAULT)
        )
        self.max_library_risk_dist_threshold: int = (
            int(kwargs["MAX_LIBRARY_RISK_DIST_THRESHOLD"])
            if (kwargs.get("MAX_LIBRARY_RISK_DIST_THRESHOLD") is not None)
            else (MAX_LIBRARY_RISK_DIST_THRESHOLD_DEFAULT)
        )
        self.nuc_2_int: dict[str, int] = (
            kwargs["NUC_2_INT"]
            if (kwargs.get("NUC_2_INT", None) is not None)
            else (NUC_2_INT_DEFAULT)
        )

    @classmethod
    def load_config(cls, path: Union[str, Path], use_env: bool = False) -> Self:
        """Helper func to load in config data"""
        config = configparser.RawConfigParser()
        config.read(os.path.normpath(path))

        settings = dict(config.items("DEFAULT"))
        if len(settings) == 0:
            raise FileNotFoundError(f"cannot find config file at '{os.path.normpath(path)}'")

        nuc_2_int = (
            {
                pair.split(":")[0].strip(): int(pair.split(":")[1].strip())
                for pair in settings["NUC_2_INT"].strip().split(",")
            }
            if settings.get("NUC_2_INT", None) is not None
            else NUC_2_INT_DEFAULT
        )

        __deli_data_dir: Optional[str] = settings.get("DELI_DATA_DIR")
        if use_env:
            deli_data_dir_env = os.environ.get("DELI_DATA_DIR", None)
            if isinstance(deli_data_dir_env, str):
                __deli_data_dir = deli_data_dir_env

        return cls(
            deli_data_dir=__deli_data_dir,
            bb_mask=settings.get("BB_MASK"),
            max_index_risk_dist_threshold=settings.get("MAX_INDEX_RISK_DIST_THRESHOLD"),
            max_library_risk_dist_threshold=settings.get("MAX_LIBRARY_RISK_DIST_THRESHOLD"),
            nuc_2_int=nuc_2_int,
        )


def _validate_deli_data_dir(deli_data_dir: Union[str, Path]) -> bool:
    if not os.path.exists(deli_data_dir):
        raise DeliDataDirError(
            f"DELi data directory '{deli_data_dir}' does not exist; "
            f"create it using 'deli create_data_dir {deli_data_dir}'"
        )

    sub_dirs = os.listdir(deli_data_dir)
    missing_sub_dirs = set(DELI_DATA_SUB_DIRS) - set(sub_dirs)
    if len(missing_sub_dirs) > 0:
        raise DeliDataDirError(
            f"DELi data directory '{deli_data_dir}' is invalid; "
            f"missing sub-directories {list(missing_sub_dirs)}"
            f"use 'deli create_data_dir {deli_data_dir} --fix'"
        )
    return True


def init_deli_data_directory(
    path: Union[str, os.PathLike],
    fail_on_exist: bool = True,
    create_default_hamming_files: bool = True,
    use_extra_parity: bool = True,
):
    """
    Create a new Deli Data Directory with all the necessary sub folders

    Parameters
    ----------
    path: Union[str, os.PathLike]
        path to create deli data dir at
    fail_on_exist: bool, default = True
        if True, will raise an exception if the directory path already exists
        if False, will try and create the sub-dirs in the existing directory
    create_default_hamming_files: bool, default = True
        if True, will create the default hamming files in the new directory
        if False, will not create default hamming files
    use_extra_parity: bool, default = True
        if True, will create extra parity bit hamming files
        if False, will not create parity bit hamming files
        only relevant if create_default_hamming_files is True
    """
    _path = Path(path)
    os.makedirs(_path, exist_ok=not fail_on_exist)
    os.makedirs(_path / "libraries", exist_ok=not fail_on_exist)
    os.makedirs(_path / "indexes", exist_ok=not fail_on_exist)
    os.makedirs(_path / "barcodes", exist_ok=not fail_on_exist)
    os.makedirs(_path / "hamming", exist_ok=not fail_on_exist)
    os.makedirs(_path / "building_blocks", exist_ok=not fail_on_exist)

    # create hamming files
    if create_default_hamming_files:
        for i in range(5, 16):
            file_path = _path / "hamming" / f"hamming3_{i}.txt"

            hamming_order = []
            for j in range(1, i + 1):
                hamming_order.append(f"p{j}" if math.log2(j).is_integer() else f"d{j}")
            extra_hamming_order = ["p0"] + hamming_order[:-1]

            with open(file_path, "w") as f:
                f.write(f"hamming_order: {','.join(hamming_order)}\n")
                f.write(f"custom_order: {','.join(hamming_order)}\n")

            if use_extra_parity:
                extra_parity_file_path = _path / "hamming" / f"hamming4_{i}.txt"
                with open(extra_parity_file_path, "w") as f:
                    f.write(f"hamming_order: {','.join(extra_hamming_order)}\n")
                    f.write(f"custom_order: {','.join(extra_hamming_order)}\n")


def init_deli_config_dir(
    path: Optional[Union[str, os.PathLike]] = None,
    fail_on_exist: bool = True,
    include_deli_data_dir: bool = True,
    create_default_hamming_files: bool = True,
    use_extra_parity: bool = True,
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
    include_deli_data_dir: bool, default = True
        if True create a deli data directory within the deli config dir
        will be named 'deli_data'
    create_default_hamming_files: bool, default = True
        if True, will create the default hamming files in the new directory
        if False, will not create default hamming files
    use_extra_parity: bool, default = True
        if True, will create extra parity bit hamming files
        if False, will not create parity bit hamming files
        only relevant if create_default_hamming_files is True
    """
    if path is not None:
        _path = Path(path)
    else:
        _path = Path(os.path.join(os.path.expanduser("~"), ".deli"))

    os.makedirs(_path, exist_ok=not fail_on_exist)

    _config = (
        f"[DEFAULT]\n"
        f"BB_MASK = {BB_MASK_TOKEN_DEFAULT}\n"
        f"MAX_INDEX_RISK_DIST_THRESHOLD = {MAX_INDEX_RISK_DIST_THRESHOLD_DEFAULT}\n"
        f"MAX_LIBRARY_RISK_DIST_THRESHOLD = {MAX_LIBRARY_RISK_DIST_THRESHOLD_DEFAULT}\n"
        f"NUC_2_INT = "
        + ",".join([f"{key}:{val}" for key, val in NUC_2_INT_DEFAULT.items()])
        + "\n"
    )
    with open(_path / ".deli", "w") as f:
        f.write(_config)

    if include_deli_data_dir:
        init_deli_data_directory(
            path=_path / "deli_data",
            create_default_hamming_files=create_default_hamming_files,
            use_extra_parity=use_extra_parity,
        )


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
    DELi config object that defines where DELI_DATA_DIR is

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

        _sub_dir_path = DELI_CONFIG.deli_data_dir / sub_dir
        if not os.path.exists(_sub_dir_path):
            raise DeliDataNotFound(
                f"cannot find DELi data subdirectory at `{_sub_dir_path}`; "
                f"did you check that DELI_DATA_DIR is set correctly?"
            )

        _file_name = os.path.basename(path).split(".")[0] + "." + extension
        file_path = DELI_CONFIG.deli_data_dir / sub_dir / _file_name

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


# load the DELi config
_deli_config_dir = os.environ.get("DELI_CONFIG", None)
if _deli_config_dir is not None and _deli_config_dir != "":
    DELI_CONFIG = _DeliConfig.load_config(_deli_config_dir, use_env=True)
elif os.path.exists(os.path.join(os.path.expanduser("~"), ".deli")):
    DELI_CONFIG = _DeliConfig.load_config(
        Path(os.path.join(os.path.expanduser("~"), ".deli")), use_env=True
    )
else:
    _deli_data_dir = os.environ.get("DELI_DATA_DIR", None)
    deli_data_dir: Optional[str]
    if isinstance(_deli_data_dir, str) and _deli_data_dir != "":
        deli_data_dir = os.fspath(os.path.expanduser(_deli_data_dir))
    else:
        deli_data_dir = None
    DELI_CONFIG = _DeliConfig(DELI_DATA_DIR=deli_data_dir)
