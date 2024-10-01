"""defines DEL library functions and classes"""
import json
from typing import Iterator, List, Optional, Self

from .barcode import BarcodeSchema
from .building_block import BuildingBlockSet
from .configure import check_file_path


class LibraryBuildError(Exception):
    """error raised when a library build fails"""

    pass


class Reaction:
    def __init__(
            self,
            cycle_id_1: str,
            cycle_id_2: str,
            reaction: str
    ):
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
        reactions: List[Reaction    ],
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
    def read_json(cls, path: str) -> Self:
        path = check_file_path(path, "libraries")
        data = json.load(open(path))

        return cls(
            library_id=data["library_id"],
            library_dna_tag=data["library_tag"],
            dna_barcode_on=data["dna_barcode_on"],
            barcode_schema=BarcodeSchema.load_from_json(data["barcode_schema"]),
            bb_sets=[BuildingBlockSet.from_csv(bb) for bb in data["bb_sets"]],
            reactions=[Reaction(**react) for react in data["reactions"]],
            scaffold=data["scaffold"]
        )

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


def get_min_index_distance(included_index: list[DELibrary]) -> int:
    """
    Given a list of index IDs, determine the minimum Levenshtein distance between all of them

    Parameters
    ----------
    included_index: list[Index]
        the index_ids to use

    Returns
    -------
    min_distance: int
        the minimum Levenshtein distance between the passed indexes
    """
    # if only one index just give it a constant threshold
    if len(included_index) == 1:
        return MAX_INDEX_RISK_DIST_THRESHOLD
    if included_index is None:
        return 0
    _index_sequences = [i.dna_tag for i in included_index]
    return min([distance(s1, s2) for s1 in _index_sequences for s2 in _index_sequences if s1 != s2])
