"""Code for handling loading/validating DELi configs for decoding"""

import configparser
import dataclasses
import functools
import os
from typing import Any, Callable, Concatenate, Literal, ParamSpec, Self, TypeVar


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

    deli_data_dir: str

    bb_mask: str
    bb_null: str

    max_index_risk_dist_threshold: int
    max_library_risk_dist_threshold: int

    def __post_init__(self):
        """Validate the passed DELi config parameters"""
        # make sure fields are not `None` (missing config data)
        for field in dataclasses.fields(self.__class__):
            if self.__getattribute__(field.name) is None:
                raise ValueError(f"Missing DELi configuration setting: {field.name}")

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


class DeliDataNotFound(Exception):
    """raised when a file cannot be found in DELi data directory"""

    pass


def accept_deli_data(
    sub_dir: Literal["building_blocks", "libraries", "indexes", "barcodes"], extension: str
) -> Callable[
    [Callable[Concatenate[Any, str, P], R]], Callable[Concatenate[Any, str, DeliConfig, P], R]
]:
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

    def _decorator(
        func: Callable[Concatenate[Any, str, P], R],
    ) -> Callable[Concatenate[Any, str, DeliConfig, P], R]:
        @functools.wraps(func)
        def _inner_func(
            cls: Any, name_or_path: str, deli_config: DeliConfig, *args: P.args, **kwargs: P.kwargs
        ) -> R:
            # check if passed name_or_path is a valid path
            if not (os.path.exists(name_or_path) and name_or_path.endswith(f".{extension}")):
                # check for valid subdir

                name_or_path = os.path.basename(name_or_path)

                _sub_path = os.path.join(os.path.abspath(deli_config.deli_data_dir), sub_dir)
                if not os.path.exists(_sub_path):
                    raise DeliDataNotFound(
                        f"cannot find DELi data subdirectory at `{_sub_path}`; "
                        f"did you check that DELI_DATA_DIR is set correctly?"
                    )
                _path = os.path.join(_sub_path, name_or_path + "." + extension)

                # check that a file with the correct name exists
                if not os.path.exists(_path):
                    raise DeliDataNotFound(
                        f"cannot find file '{name_or_path}.{extension}' in {_sub_path}"
                    )
            # name_or_path was a valid path
            else:
                _path = name_or_path

            return func(cls, _path, *args, **kwargs)

        return _inner_func

    return _decorator
