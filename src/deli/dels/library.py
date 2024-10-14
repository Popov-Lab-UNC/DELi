"""defines DEL library functions and classes"""

import json
from typing import Iterator, List, Optional, Self, Union

from Levenshtein import distance

from .barcode import BarcodeSchema
from .building_block import BuildingBlockSet
from .configure import check_file_path


class LibraryBuildError(Exception):
    """error raised when a library build fails"""

    pass


class Reaction:
    """struct to contain info on DEL reactions"""

    def __init__(self, cycle_id_1: str, cycle_id_2: str, reaction: str):
        self.cycle_id_1 = cycle_id_1
        self.bb_set_id_2 = cycle_id_2
        self.reaction = reaction


class DELibrary:
    """
    Contains information about a DNA-encoded library

    Attributes
    ----------
    num_cycles : int
        Number of building block cycles in the library
    library_size : int
        Size of the enumerated library
    """

    def __init__(
        self,
        library_id: str,
        library_dna_tag: str,
        barcode_schema: BarcodeSchema,
        bb_sets: List[BuildingBlockSet],
        reactions: List[Reaction],
        dna_barcode_on: str,
        scaffold: Optional[str] = None,
    ):
        """
        Initialize a DELibrary object

        Parameters
        ----------
        library_id : str
            name/id of the library
        library_dna_tag : str
            the DNA sequence associated with this library
        barcode_schema : BarcodeSchema
            The barcode schema defining how the barcodes are designed
        bb_sets : List[BuildingBlockSet]
            the sets of building-block used to build this library
            order in list should be order of synthesis
            must have length >= 2
        reactions : List[str]
            the reaction SMARTS/SMIRKS that connect bb_set cycles
            reaction at index i is the reaction between bb_sets[i] and bb_sets[i+1]
            reactions must have length equal to the the number of `bb_sets` minus 1
        dna_barcode_on: str
            the id of the bb_set that is linked to the DNA tag
            can be 'scaffold' if DNA tag is linked to the scaffold
        scaffold: optional, str
            SMILES of the scaffold
            if no scaffold in library should be `None`

        Raises
        ------
        LibraryBuildError
            the parameters passed to build the library are not compatible with each other
            error message will contain specific details of build issue
        """
        self.library_id = library_id
        self.library_tag = library_dna_tag
        self.barcode_schema = barcode_schema
        self.bb_sets = bb_sets
        self.reactions = reactions
        self.dna_barcode_on = dna_barcode_on
        self.scaffold = scaffold

        if len(self.bb_sets) < 2:
            raise LibraryBuildError(
                f"Library requires at least 2 cycle (bb_sets);" f" found {len(self.bb_sets)}"
            )
        self.num_cycles = len(self.bb_sets)

        if self.num_cycles - (1 if self.scaffold is None else 0) < len(self.reactions):
            raise LibraryBuildError(
                f"Library requires at least N-1 reactions for N cycles+scaffolds; "
                f"found {len(self.bb_sets)} cycles and {len(self.reactions)} reactions"
            )

        if self.num_cycles != barcode_schema.num_cycles:
            raise LibraryBuildError(
                f"Number of library cycles does not match barcode schema cycles; "
                f"got {self.num_cycles} and {barcode_schema.num_cycles}"
            )

        self.library_size = sum([len(bb_set) for bb_set in self.bb_sets])

        # handle the dna tag location
        if self.dna_barcode_on not in [bb_set.bb_set_id for bb_set in self.bb_sets]:
            if scaffold is not None and self.dna_barcode_on != "scaffold":
                raise LibraryBuildError(
                    f"cannot find cycle {self.dna_barcode_on} to put DNA barcode on"
                )
            if scaffold is None and self.dna_barcode_on == "scaffold":
                raise LibraryBuildError("no scaffold to attach DNA barcode to")

    @classmethod
    def from_dict(cls, lib_dict: dict) -> Self:
        """
        Load a DEL from a dict

        Parameters
        ----------
        lib_dict: dict
            DEL info as a dictionary

        Returns
        -------
        DELibrary
        """
        return cls(
            library_id=lib_dict["library_id"],
            library_dna_tag=lib_dict["library_tag"],
            dna_barcode_on=lib_dict["dna_barcode_on"],
            barcode_schema=BarcodeSchema.load_from_json(lib_dict["barcode_schema"]),
            bb_sets=[BuildingBlockSet.from_csv(bb) for bb in lib_dict["bb_sets"]],
            reactions=[Reaction(**react) for react in lib_dict["reactions"]],
            scaffold=lib_dict["scaffold"],
        )

    @classmethod
    def read_json(cls, path: str) -> Self:
        """
        Load a DEL from a json file

        Parameters
        ----------
        path: str
            path to file with DEL json

        Returns
        -------
        DELibrary
        """
        path = check_file_path(path, "libraries")
        data = json.load(open(path))

        return cls.from_dict(data)

    def iter_bb_sets(self) -> Iterator[BuildingBlockSet]:
        """
        Iterate through building block sets that make up library

        Yields
        ------
        BuildingBlockSet
        """
        for bb_set in self.bb_sets:
            if bb_set:
                yield bb_set


class MegaDELibrary:
    """A set of many DELibraries"""

    def __init__(self, libraries: List[DELibrary]):
        """
        Initialize a MegaDELibrary object

        Parameters
        ----------
        libraries: List[DELibrary]
            libraries to include in the mega library
        """
        self.libraries = libraries

    def __len__(self) -> int:
        """Return number of libraries in the mega library"""
        return len(self.libraries)

    def __iter__(self) -> Iterator[DELibrary]:
        """Iterate through all libraries in the mega library"""
        return iter(self.libraries)

    @classmethod
    def read_json(cls, path: str) -> Self:
        """
        Read a mega library from a json file

        Parameters
        ----------
        path: str
            path to mega library json

        Returns
        -------
        MegaDELibrary
        """
        path = check_file_path(path, "libraries")
        data = json.load(open(path))

        return cls(libraries=[DELibrary.from_dict(d) for d in data])


def get_min_library_tag_distance(
    included_libraries: Optional[Union[list[DELibrary], MegaDELibrary]],
) -> int:
    """
    Determine the minimum Levenshtein distance between all library tags

    Parameters
    ----------
    included_libraries: list[DELibrary] or MegaDELibrary
        the libraries to use

    Returns
    -------
    min_distance: int
        the minimum Levenshtein distance between the passed library's tag
    """
    # if only one index just give it a constant threshold
    if included_libraries is None or len(included_libraries) == 1:
        return 0
    _lib_dna_sequences = [i.library_tag for i in included_libraries]
    return min(
        [distance(s1, s2) for s1 in _lib_dna_sequences for s2 in _lib_dna_sequences if s1 != s2]
    )
