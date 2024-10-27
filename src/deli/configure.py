"""Code for handling loading/validating DELi configs for decoding"""

import configparser
import dataclasses
import functools
import inspect
import os
from pathlib import Path
from typing import Any, Callable, Literal, ParamSpec, Self, TypeVar, Union


P = ParamSpec("P")
R = TypeVar("R")


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

    def __post_init__(self):
        """Validate the passed DELi config parameters"""
        # make sure fields are not `None` (missing config data)
        for field in dataclasses.fields(self.__class__):
            if self.__getattribute__(field.name) is None:
                _system_var = os.getenv(field.name.upper())
                if _system_var is None:
                    raise ValueError(f"Missing DELi configuration setting: {field.name}")
                else:
                    self.__setattr__(field.name, _system_var)

    @classmethod
    def _load_config(cls, path: str) -> dict[str, Any]:
        """Helper func to load in config data"""
        config = configparser.RawConfigParser()
        config.read(os.path.normpath(path))

        settings = dict(config.items("DEFAULT"))

        # clean up the params
        _clean_settings: dict[str, Any] = {
            "deli_data_dir": os.fspath(os.path.expanduser(settings["deli_data_dir"])),
            "max_index_risk_dist_threshold": int(settings["max_index_risk_dist_threshold"]),
            "max_library_risk_dist_threshold": int(settings["max_library_risk_dist_threshold"]),
        }
        for key, val in settings.items():
            if key not in _clean_settings.keys():
                _clean_settings[key] = val

        return _clean_settings

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
        return cls(**cls._load_config(os.path.join(os.path.expanduser("~"), ".deli", ".deli")))

    @classmethod
    def load_from_file(cls, file_path: str) -> Self:
        """
        Load DELi configuration settings from a given file path

        Notes
        -----
        DELi configuration settings follow the python configure
        format

        Parameters
        ----------
        file_path: str
            path to DELi configuration file

        Returns
        -------
        DeliConfig
        """
        return cls(**cls._load_config(path=file_path))


DELI_CONFIG = DeliConfig.load_defaults()


class DeliDataNotFound(Exception):
    """raised when a file cannot be found in DELi data directory"""

    pass


def accept_deli_data_name(
    sub_dir: Literal["building_blocks", "libraries", "indexes", "barcodes"], extension: str
) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """
    Decorator to allow load functions to take name of DELi data object

    Notes
    -----
    This should be used decorate the `load` function of any
    class that has the `DeliDataLoadableMixin` mixin

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
