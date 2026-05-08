"""Parsing interface for loading DEL objects from the DELi data directory"""

import abc
from typing import Self

import fsspec

from deli.config import DELI_DATA_SUB_DIRS, get_deli_config


class Parsable(abc.ABC):
    """
    Interface for objects that can be parsed from a string representation, either from a URI or from an ID

    DEL objects that can be loaded from files in the DELi data directory should implement this interface.
    Classes that implement it only need to implement logic to construct the object from a string representation
    of the file contents, and the interface will handle reading the file contents for the correct locations.

    When defining a new Parsable subclass, you must specify the subdirectory in the DELi data directory where
    objects of that class are stored.
    """

    _sub_dir: str

    def __init_subclass__(cls, sub_dir: str | None = None, **kwargs):
        """Set the subdirectory for a Parsable subclass"""
        super().__init_subclass__(**kwargs)
        if sub_dir is None:
            raise TypeError("Subdirectory must be specified for Parsable subclasses")
        if sub_dir not in DELI_DATA_SUB_DIRS:
            raise ValueError(f"Subdirectory '{sub_dir}' is a known DELi data subdirectory")
        cls._sub_dir = sub_dir

    # when implementing a Parsable subclass, you should override the load method if you are modifying its signature.
    # This makes its easier for IDEs and users to understand how to use the load method for that subclass.
    @classmethod
    def load(
        cls,
        uri_or_id: str | None = None,
        uri: str | None = None,
        id: str | None = None,
        **kwargs
    ) -> Self:
        """
        Load class from the DELi data directory or from a file

        Is capable of deciding whether the input is a URI or an ID and loading accordingly.
        If input types is known at time of calling, can specify directly using `uri` or `id`
        parameters to avoid ambiguity.

        Only one of `uri_or_id`, `uri`, or `id` should be provided.

        Notes
        -----
        When loading with an ID, DELi will search for the objects config file in
        the given subdirectory assigned to the object class.

        Parameters
        ----------
        uri_or_id: str | None
            Either a URI or an ID. DELi will attempt to resolve whether the input is a URI or an ID.
        uri: str | None
            A URI (or local path) to load object from.
        id: str | None
            The ID of the object to load from the DELi data directory.
        validate_chemicals: bool
            Whether to validate the chemicals in the object after loading. This will assert that
            elements within the object that could have chemical structures do have them and they
            are valid.
        """
        _config = get_deli_config()

        if sum(x is not None for x in [uri_or_id, uri, id]) != 1:
            raise ValueError("Exactly one of `uri_or_id`, `uri`, or `id` should be provided.")
        elif uri_or_id is not None:
            # attempt to resolve whether input is a URI or an ID
            try:
                data = _config.data_directory_filesystem.read_subdir_file_by_id(cls._sub_dir, uri_or_id)
            except Exception:
                # if loading by ID fails, attempt to load by URI
                data = fsspec.open(uri_or_id).read_text()
        elif uri is not None:
            data = fsspec.open(uri).read_text()
        elif id is not None:
            data = _config.data_directory_filesystem.read_subdir_file_by_id(cls._sub_dir, id)
        else:
            raise ValueError("Exactly one of `uri_or_id`, `uri`, or `id` should be provided.")

        return cls._loads(data, **kwargs)

    @classmethod
    @abc.abstractmethod
    def _loads(cls, data: str, **kwargs) -> Self:
        """
        Construct object from a string representation

        When defining in a subclass, any additional kwargs that are given in
        load through load will be passed through to this function.
        """
        raise NotImplementedError()
