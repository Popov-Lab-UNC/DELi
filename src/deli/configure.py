"""Code for handling loading/validating DELi configs for decoding"""

import abc
import configparser
import dataclasses
import functools
import getpass
import inspect
import math
import os
import warnings
from pathlib import Path
from typing import Any, Callable, Dict, Literal, ParamSpec, Self, TypeVar, Union


P = ParamSpec("P")
R = TypeVar("R")

CONFIG_DIR_PATH = Path(os.path.join(os.path.expanduser("~"), ".deli"))

DELI_DATA_SUB_DIRS = ["hamming", "barcodes", "libraries", "indexes", "building_blocks"]


class DeliDataNotFound(Exception):
    """raised when a file cannot be found in DELi data directory"""

    pass


class DeliDataDirError(Exception):
    """raised when a DELi data directory cannot be found or is invalid (missing sub-dirs)"""

    pass


@dataclasses.dataclass
class DeliConfig:
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
    bb_null: str
        the 3 char long token to replace any null building blocks
    max_index_risk_dist_threshold: int
        the maximum risk willing to make in an index call
        if risk is above this number will fail the call
    max_library_risk_dist_threshold: int
        the maximum risk willing to make in a library call
        if risk is above this number will fail the call
    """

    deli_data_dir: Path

    bb_mask: str
    bb_null: str

    max_index_risk_dist_threshold: int
    max_library_risk_dist_threshold: int

    nuc_2_int: Dict[str, int]

    def __post_init__(self):
        """Validate the passed DELi config parameters"""
        self.deli_data_dir = Path(os.fspath(os.path.expanduser(self.deli_data_dir)))
        self.bb_mask = str(self.bb_mask)
        self.bb_null = str(self.bb_null)
        self.max_index_risk_dist_threshold = int(self.max_index_risk_dist_threshold)
        self.max_library_risk_dist_threshold = int(self.max_library_risk_dist_threshold)
        self.nuc_2_int = {
            pair.split(":")[0].strip(): int(pair.split(":")[1].strip())
            for pair in self.nuc_2_int.strip().split(",")
        }

        # check for a valid Deli Data Dir
        if not os.path.exists(self.deli_data_dir):
            raise DeliDataDirError(
                f"DELi data directory '{self.deli_data_dir}' does not exist; "
                f"create it using 'deli create_data_dir {self.deli_data_dir}' OR "
                f"change your '.deli' config file to point to a valid DELi data directory"
            )

        sub_dirs = os.listdir(self.deli_data_dir)
        missing_sub_dirs = set(DELI_DATA_SUB_DIRS) - set(sub_dirs)
        if len(missing_sub_dirs) > 0:
            raise DeliDataDirError(
                f"DELi data directory '{self.deli_data_dir}' is invalid; "
                f"missing sub-directories {list(missing_sub_dirs)}"
                f"create missing sub-directories manually OR "
                f"use 'deli create_data_dir {self.deli_data_dir} --fix' OR "
                f"change your '.deli' config file to point to a valid DELi data directory"
            )

        os.environ["DELI_DATA_DIR"] = str(self.deli_data_dir)

    @classmethod
    def _load_config(cls, path: Union[str, Path]) -> dict[str, Any]:
        """Helper func to load in config data"""
        config = configparser.RawConfigParser()
        config.read(os.path.normpath(path))

        settings = dict(config.items("DEFAULT"))
        if len(settings) == 0:
            raise FileNotFoundError(f"cannot find config file at '{os.path.normpath(path)}'")

        for field in dataclasses.fields(cls):
            if settings.get(field.name) is None:
                _system_var = os.getenv(field.name.upper())
                if _system_var is None:
                    raise ValueError(f"Missing DELi configuration setting: {field.name}")
                else:
                    settings[field.name] = _system_var

        return settings

    @classmethod
    def load_defaults(cls) -> Self:
        """
        Load default DELi configuration settings

        Notes
        -----
        Default DELi configuration files are saved in
        `~/.deli/.deli` where `~` if your USER home directory

        If a given DELI parameter is unfilled in the config
        file then DELI will try and load it from the system
        arguments (if they are set). System args should be
        all caps versions of the variables

        Returns
        -------
        DeliConfig
        """
        return cls(**cls._load_config(CONFIG_DIR_PATH / ".deli"))


def init_deli_data_dir(path: Union[str, os.PathLike], fail_on_exist: bool = True):
    """
    Create a new Deli Data Directory with all the necessary sub folders

    Parameters
    ----------
    path: Union[str, os.PathLike]
        path to create deli data dir at
    fail_on_exist: bool, default = True
        if True, will raise an exception if the directory path already exists
        if False, will try and create the sub-dirs in the existing directory
    """
    _path = Path(path)
    os.makedirs(_path, exist_ok=not fail_on_exist)
    os.makedirs(_path / "libraries", exist_ok=not fail_on_exist)
    os.makedirs(_path / "indexes", exist_ok=not fail_on_exist)
    os.makedirs(_path / "barcodes", exist_ok=not fail_on_exist)
    os.makedirs(_path / "hamming", exist_ok=not fail_on_exist)
    os.makedirs(_path / "building_blocks", exist_ok=not fail_on_exist)
    _create_default_hamming_files(_path / "hamming")


def _create_default_hamming_files(path: Path):
    """Will generate generic hamming files for all hamming coded tags length 5-16"""
    for i in range(5, 16):
        file_path = path / f"hamming3_{i}.txt"
        extra_parity_file_path = path / f"hamming4_{i}.txt"

        hamming_order = []
        for j in range(1, i + 1):
            hamming_order.append(f"p{j}" if math.log2(j).is_integer() else f"d{j}")
        extra_hamming_order = ["p0"] + hamming_order[:-1]

        with open(file_path, "w") as f:
            f.write(f"hamming_order: {','.join(hamming_order)}\n")
            f.write(f"custom_order: {','.join(hamming_order)}\n")

        with open(extra_parity_file_path, "w") as f:
            f.write(f"hamming_order: {','.join(extra_hamming_order)}\n")
            f.write(f"custom_order: {','.join(extra_hamming_order)}\n")


def _create_deli_config_default(path: Path):
    """Create a default .deli config file"""
    _config = (
        f"[DEFAULT]\n"
        f"DELI_DATA_DIR = {CONFIG_DIR_PATH / 'deli_data'}\n"
        f"BB_MASK = ###\n"
        f"BB_NULL = NUL\n"
        f"MAX_INDEX_RISK_DIST_THRESHOLD = 3\n"
        f"MAX_LIBRARY_RISK_DIST_THRESHOLD = 4\n"
        f"NUC_2_INT = A:0,T:1,C:2,G:3\n"
    )
    with open(path, "w") as f:
        f.write(_config)


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


# load in the deli config file
if not os.path.exists(CONFIG_DIR_PATH):
    warnings.warn(
        f"user {getpass.getuser()} has no '.deli' config directory in their home directory; "
        f"creating one",
        stacklevel=0,
    )
    os.makedirs(CONFIG_DIR_PATH, exist_ok=True)
    _create_deli_config_default(CONFIG_DIR_PATH / ".deli")
    init_deli_data_dir(CONFIG_DIR_PATH / "deli_data")

try:
    DELI_CONFIG = DeliConfig.load_defaults()
except FileNotFoundError:
    warnings.warn(
        f"user {getpass.getuser()} has no '.deli' config file in their user directory; "
        f"creating default",
        stacklevel=0,
    )
    _create_deli_config_default(CONFIG_DIR_PATH / ".deli")
