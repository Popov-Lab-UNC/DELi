"""defines base DEL selection objects"""

import abc
import random
import time
from os import PathLike
from typing import Any, Self

import yaml

from deli.dels import DELibrary, DELibraryPool

from .settings import DecodingSettings


class SelectionSyntaxError(Exception):
    """Exception to raise when a syntax error if found in a selection file"""

    pass


class BaseSelection(abc.ABC):
    """
    Base class for all children that deal with selection data/settings

    This can go beyond just a "selection".
    For example, all decoding and analysis experiment definitions are
    children of selection, as they are a selection with extra setting.
    """

    def __init__(
        self,
        library_pool: DELibraryPool,
        target_id: str | None = None,
        selection_id: str | None = None,
    ):
        """
        Initialize the BaseSelection object

        Notes
        -----
        DELi assumes all data is already demultiplex at the target level
        so every selection should only have one target ID

        DELi assumes uniqueness of the target ID, but it is never
        strictly enforced. It would not impact anything with decoding
        but it could mess with some downstream analysis

        Parameters
        ----------
        library_pool: DELibraryPool
            the libraries used in the selection
        target_id: str
            the id of the target used in the selection
            if `None` will default to "NA"
        selection_id: str, optional
            the id of the selection.
            if `None` will default to a random number based on
            timestamp the object was created
        """
        self.library_pool = library_pool
        self.target_id: str = target_id if target_id is not None else "NA"
        self.selection_id: str = (
            selection_id
            if selection_id
            else (str(time.time()).replace(".", "") + f"{random.randint(0, 1000000):06d}")
        )

    @abc.abstractmethod
    def to_file(self, out_path: str | PathLike):
        """Covert the object to a human readable file"""
        raise NotImplementedError()

    @classmethod
    @abc.abstractmethod
    def from_file(cls, path: str | PathLike) -> Self:
        """Load the object from a human readable file"""
        raise NotImplementedError()


class DecodingExperiment(BaseSelection):
    """
    Defines a decoding experiment for a DEL selection

    Configure using the decoding settings
    """

    def __init__(
        self,
        library_pool: DELibraryPool,
        target_id: str,
        decode_settings: DecodingSettings,
        selection_id: str | None = None,
    ):
        """
        Initialize the experiment with the given settings

        Parameters
        ----------
        library_pool: DELibraryPool
            the library pool used in the selection
        target_id: str
            the id of the target used in the selection
        decode_settings: DecodingSettings
            Settings to use for decoding
        selection_id: str or None, default = None
            the id of the selection
            if `None` will default to a random number based on
            timestamp the object was created
        """
        super().__init__(library_pool, target_id, selection_id=selection_id)
        self.decode_settings = decode_settings

    def to_file(self, out_path: str | PathLike):
        """
        Write experiment to human readable file

        Parameters
        ----------
        out_path: str or PathLike
            path to save experiment to
        """
        data = {
            "target_id": self.target_id,
            "selection_id": self.selection_id,
            "libraries": [lib.loaded_from for lib in self.library_pool],
            "decode_settings": self.decode_settings.__dict__,
        }
        yaml.safe_dump(data, open(out_path, "w"))

    @classmethod
    def from_file(cls, file_path: str | PathLike) -> Self:
        """
        Load the experiment from a human readable file

        Parameters
        ----------
        file_path: str or PathLike
            path to load experiment from

        Returns
        -------
        DecodingExperiment
        """
        data = yaml.safe_load(open(file_path, "r"))

        _selection_id = data.get("selection_id", None)
        _target_id = data.get("target_id", None)

        try:
            _libraries: list[str] = data["libraries"]
        except KeyError as e:
            raise SelectionSyntaxError(
                f"{file_path} decoding file does not contain a 'libraries' section"
            ) from e

        _library_pool = DELibraryPool([DELibrary.load(lib_path) for lib_path in _libraries])

        _decode_settings: dict[str, Any] | None = data.get("decode_settings", None)
        if _decode_settings is None:
            _decode_setting_obj = DecodingSettings()
        else:
            try:
                _decode_setting_obj = DecodingSettings(**_decode_settings)
            except TypeError as e:
                _unknown_arg = e.args[0].split()[-1]
                raise SelectionSyntaxError(
                    f"unrecognized decoding settings: {_unknown_arg}"
                ) from e

        return cls(
            selection_id=_selection_id,
            target_id=_target_id,
            library_pool=_library_pool,
            decode_settings=_decode_setting_obj,
        )
