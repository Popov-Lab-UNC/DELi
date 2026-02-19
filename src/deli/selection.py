"""handles DEL selection information"""

import os
from pathlib import Path
from typing import Optional, Sequence

import yaml

from deli.dels.combinatorial import DELibrary, DELibraryCollection
from deli.dels.tool_compound import TaggedToolCompound
from deli.dna.io import SequenceReader, get_reader


class Selection:
    """Represents a selection made with a DEL library collection"""

    def __init__(
        self,
        selection_id: str,
        library_collection: DELibraryCollection,
        tool_compounds: Optional[Sequence[TaggedToolCompound]] = None,
    ):
        """
        Initialize a Selection object.

        Parameters
        ----------
        selection_id: str
            the unique id for this selection; unqiueness is not programmatically enforced
            but recommended for tracking selections
        library_collection: DELibraryCollection
            the DEL collection used in the selection
        tool_compounds: Sequence[TaggedToolCompound], optional
            the tagged tool compound libraries used in the selection
        """
        self.library_collection: DELibraryCollection = library_collection
        self.tool_compounds: Sequence[TaggedToolCompound] = tool_compounds if tool_compounds else []
        self.selection_id = selection_id


class SequencedSelection(Selection):
    """
    Represents a selection made with a DEL collection that has been sequenced

    Parameters
    ----------
    selection_id: str
        the unique id for this selection; unqiueness is not programmatically enforced
        but recommended for tracking selections
    library_collection: DELibraryCollection
        the library collection used in the selection
    sequence_reader: SequenceReader
        the list of barcode files for the selection
    tool_compounds: Sequence[TaggedToolCompoundLibrary], optional
        the tagged tool compound libraries used in the selection
    """

    def __init__(
        self,
        selection_id: str,
        library_collection: DELibraryCollection,
        sequence_reader: SequenceReader,
        tool_compounds: Optional[Sequence[TaggedToolCompound]] = None,
    ):
        """Initialize a Selection object."""
        super().__init__(
            selection_id=selection_id,
            library_collection=library_collection,
            tool_compounds=tool_compounds,
        )
        self.sequence_reader: SequenceReader = sequence_reader

    @property
    def sequence_files(self) -> tuple[Path, ...]:
        """The sequence files associated with this selection"""
        return self.sequence_reader.get_sequence_files()


def load_selection(path: os.PathLike, load_chemical_info: bool = True) -> Selection:
    """
    Load a Selection object from a YAML file.

    See the Selection file docs for info on the expected YAML format.
    Will load as a `SequencedSelection` if "sequence_files" key is present in the YAML,
    otherwise loads as a `Selection`.

    Notes
    -----
    If the selection file is missing an ID, the ID will be generated as the file hash.
    This will not be saved anywhere, but may be used in naming output files related to
    the selection. It is recommended to provide an ID in the selection file to avoid
    confusion.

    Parameters
    ----------
    path: os.PathLike
        the path to the YAML file containing the selection data
    load_chemical_info: bool, default = True
        Whether to load any present chemical information for libraries and tool compounds in selection.

    Returns
    -------
    Selection
        the loaded Selection object
    """
    data = yaml.safe_load(open(path, "r"))

    lib_collection = DELibraryCollection(
        [DELibrary.load(lib, load_chemical_info=load_chemical_info) for lib in data["libraries"]]
    )
    tool_compounds = [
        TaggedToolCompound.load(tc, load_smiles=load_chemical_info) for tc in data.get("tool_compounds", [])
    ]

    if "selection_id" not in data:
        data["selection_id"] = str(hash(str(data)))

    if "sequence_files" in data:
        sequence_reader = get_reader(data["sequence_files"])
        return SequencedSelection(
            selection_id=data["selection_id"],
            library_collection=lib_collection,
            sequence_reader=sequence_reader,
            tool_compounds=tool_compounds,
        )
    else:
        return Selection(
            selection_id=data["selection_id"],
            library_collection=lib_collection,
            tool_compounds=tool_compounds,
        )
