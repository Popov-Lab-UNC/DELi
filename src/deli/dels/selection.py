"""handles DEL selection information"""

import os
from datetime import datetime
from typing import Any

import yaml

from deli.dna import SequenceGlobReader, SequenceReader

from .library import DELibrary, DELibraryCollection, Library, LibraryCollection


class SectionCondition:
    """Represents a selection condition for a DEL selection"""

    def __init__(
        self,
        target_id: str | None = None,
        selection_condition: str | None = None,
        additional_info: str | None = None,
    ):
        """
        Initialize a SelectionCondition object.

        Parameters
        ----------
        target_id: str | None
            the target id for the selection, defaults to None if not provided
        selection_condition: str | None
            the selection condition, defaults to None if not provided
        additional_info: str | None
            any additional information about the selection, defaults to None if not provided
        """
        self.target_id = target_id or "Unknown"
        self.selection_condition = selection_condition or "Unknown"
        self.additional_info = additional_info or ""

    @classmethod
    def from_dict(cls, data: dict):
        """
        Create a SelectionCondition object from a dictionary.

        Parameters
        ----------
        data: dict
            the dictionary containing selection condition data

        Returns
        -------
        Selection
            the Selection object created from the dictionary
        """
        return cls(
            target_id=data.get("target_id"),
            selection_condition=data.get("selection_condition"),
            additional_info=data.get("additional_info"),
        )

    def to_dict(self) -> dict:
        """
        Convert the Selection object to a dictionary.

        Returns
        -------
        dict
        """
        return {
            "target_id": self.target_id,
            "selection_condition": self.selection_condition,
            "additional_info": self.additional_info,
        }


class Selection:
    """Represents a selection made with a DEL library collection"""

    def __init__(
        self,
        library_collection: LibraryCollection,
        date_ran: datetime | None = None,
        target_id: str | None = None,
        selection_condition: str | None = None,
        selection_id: str | None = None,
        additional_info: str | None = None,
    ):
        """
        Initialize a Selection object.

        Parameters
        ----------
        library_collection: LibraryCollection
            the library collection used in the selection
        date_ran: datetime | None
            the date the selection was run, defaults to now if None
        target_id: str | None
            the target id for the selection, defaults to None if not provided
        selection_condition: str | None
            the selection condition, defaults to None if not provided
        selection_id: str | None
            the unique id for this selection, defaults to None if not provided
        additional_info: str | None
            any additional information about the selection, defaults to None if not provided
        """
        self.library_collection: LibraryCollection = library_collection
        self.selection_id = selection_id if selection_id else "Unknown"
        self.date_ran = date_ran

        self.selection_condition = SectionCondition(
            target_id=target_id,
            selection_condition=selection_condition,
            additional_info=additional_info,
        )

    @classmethod
    def from_dict(cls, data: dict):
        """
        Create a Selection object from a dictionary.

        Parameters
        ----------
        data: dict
            the dictionary containing selection data

        Returns
        -------
        Selection
            the Selection object created from the dictionary
        """
        lib_collection = LibraryCollection([Library.load(lib) for lib in data["libraries"]])

        return cls(
            library_collection=lib_collection,
            date_ran=datetime.fromisoformat(data["date_ran"])
            if data.get("date_ran") is not None
            else None,
            target_id=data.get("target_id"),
            selection_condition=data.get("selection_condition"),
            selection_id=data.get("selection_id"),
            additional_info=data.get("additional_info"),
        )

    @classmethod
    def from_yaml(cls, path: str | os.PathLike) -> "Selection":
        """
        Create a Selection object from a YAML file

        Parameters
        ----------
        path: str | os.PathLike
            the path to the YAML file containing selection data

        Returns
        -------
        Selection
        """
        data = yaml.safe_load(open(path, "r"))
        return cls.from_dict(data)

    def to_json_str(self) -> dict:
        """
        Convert the Selection object to a dictionary of JSON serializable objects.

        Returns
        -------
        dict
        """
        _data: dict[str, Any] = {
            "date_ran": self.get_run_date_as_str(),
            "target_id": self.selection_condition.target_id,
            "selection_condition": self.selection_condition,
            "selection_id": self.selection_id,
            "additional_info": self.selection_condition.additional_info,
            "libraries": [lib.library_id for lib in self.library_collection.libraries],
        }
        return _data

    def get_run_date_as_str(self) -> str:
        """
        Get the date the selection was run as a string

        Returns
        -------
        str
            the date the selection was run in ISO format, or "Unknown" if not set
        """
        return self.date_ran.isoformat() if self.date_ran is not None else "Unknown"


class SequencedSelection(Selection):
    """Represents a selection made with a DEL collection that has been sequenced"""

    def __init__(
        self,
        library_collection: DELibraryCollection,
        sequence_files: list[str | os.PathLike],
        date_ran: datetime | None = None,
        target_id: str | None = None,
        selection_condition: str | None = None,
        selection_id: str | None = None,
        additional_info: str | None = None,
    ):
        """
        Initialize a Selection object.

        Parameters
        ----------
        library_collection: DELibraryCollection
            the library collection used in the selection
        date_ran: datetime | None
            the date the selection was run, defaults to now if None
        target_id: str | None
            the target id for the selection, defaults to None if not provided
        selection_condition: str | None
            the selection condition, defaults to None if not provided
        selection_id: str | None
            the unique id for this selection, defaults to None if not provided
        additional_info: str | None
            any additional information about the selection, defaults to None if not provided
        """
        super().__init__(
            library_collection=library_collection,
            date_ran=date_ran,
            target_id=target_id,
            selection_condition=selection_condition,
            selection_id=selection_id,
            additional_info=additional_info,
        )
        self.library_collection: DELibraryCollection = library_collection  # for type checking
        self.sequence_files = sequence_files

        if len(self.sequence_files) < 1:
            raise ValueError(
                "A SequencedSelection must have at least 1 observed_seq file; found 0"
            )

        # validate observed_seq files
        SequenceGlobReader(self.sequence_files)

    @classmethod
    def from_dict(cls, data: dict) -> "SequencedSelection":
        """
        Create a Selection object from a dictionary.

        Parameters
        ----------
        data: dict
            the dictionary containing selection data

        Returns
        -------
        Selection
            the Selection object created from the dictionary
        """
        lib_collection = DELibraryCollection([DELibrary.load(lib) for lib in data["libraries"]])
        return cls(
            library_collection=lib_collection,
            sequence_files=data.get("sequence_files", []),
            date_ran=datetime.fromisoformat(data["date_ran"]),
            target_id=data.get("target_id"),
            selection_condition=data.get("selection_condition"),
            selection_id=data.get("selection_id"),
            additional_info=data.get("additional_info"),
        )

    def to_json_str(self) -> dict:
        """
        Convert the Selection object to a dictionary.

        Returns
        -------
        dict
        """
        _data = super().to_json_str()
        _data["sequence_files"] = self.sequence_files
        return _data

    def get_sequence_reader(self) -> SequenceReader:
        """
        Get a observed_seq reader for the selection's observed_seq files.

        Returns
        -------
        SequenceReader
            a SequenceReader object for the selection's observed_seq files
        """
        return SequenceGlobReader(self.sequence_files)
