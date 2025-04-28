"""defines DEL library functions and classes"""

import json
import os
from functools import reduce
from operator import mul
from pathlib import Path
from typing import Iterator, List, Optional, Self, Union

from deli.configure import DeliDataLoadable, accept_deli_data_name

from .barcode import BarcodeSchema, BuildingBlockBarcodeSection
from .building_block import BuildingBlock, BuildingBlockSet
from .enumerator import DELEnumerator
from .reaction import ReactionWorkflow


class LibraryBuildError(Exception):
    """error raised when a library build fails"""

    pass


class DELibrary(DeliDataLoadable):
    """
    Contains information about a DNA-encoded library

    Attributes
    ----------
    num_cycles : int
        Number of building block cycles in the library
    library_size : int
        Size of the enumerated library
    enumerator : DELEnumerator
        the enumerator attached to the library
    """

    def __init__(
        self,
        library_id: str,
        barcode_schema: BarcodeSchema,
        bb_sets: List[BuildingBlockSet],
        library_reaction_workflow: Optional[ReactionWorkflow] = None,
        dna_barcode_on: Optional[str] = None,
        scaffold: Optional[str] = None,
        anchored: bool = False,
    ):
        """
        Initialize a DELibrary object

        Parameters
        ----------
        library_id : str
            name/id of the library
        barcode_schema : BarcodeSchema
            The barcode schema defining how the barcodes are designed
        bb_sets : List[BuildingBlockSet]
            the sets of building-block used to build this library
            order in list should be order of synthesis
            must have length >= 2
        library_reaction_workflow : ReactionWorkflow
            The reaction workflow/schema used to build this library
        dna_barcode_on: str
            the id of the bb_set that is linked to the DNA bases
            can be 'scaffold' if DNA bases is linked to the scaffold
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
        self.barcode_schema = barcode_schema
        self.bb_sets = bb_sets
        self.library_reaction_workflow = library_reaction_workflow
        self.dna_barcode_on = dna_barcode_on
        self.scaffold = scaffold
        self.anchored = anchored

        self.library_tag = self.barcode_schema.library_section.get_dna_sequence()
        self.library_size = reduce(mul, [len(bb_set) for bb_set in self.bb_sets])
        self.num_cycles = len(self.bb_sets)
        self.enumerator: Optional[DELEnumerator] = None
        if isinstance(self.library_reaction_workflow, ReactionWorkflow):
            self.enumerator = DELEnumerator(
                self.library_reaction_workflow, self.bb_sets, self.scaffold
            )

        ### VALIDATION ###

        if len(self.bb_sets) < 2:
            raise LibraryBuildError(
                f"Library requires at least 2 cycle (bb_sets);" f" found {len(self.bb_sets)}"
            )

        if self.num_cycles != barcode_schema.get_num_building_block_sections():
            raise LibraryBuildError(
                f"Number of library cycles does not match barcode schema cycles; "
                f"got {self.num_cycles} and {barcode_schema.get_num_building_block_sections()}"
            )

        if isinstance(self.dna_barcode_on, str):
            if self.dna_barcode_on not in [bb_set.bb_set_id for bb_set in self.bb_sets]:
                if scaffold is not None and self.dna_barcode_on != "scaffold":
                    raise LibraryBuildError(
                        f"cannot find cycle {self.dna_barcode_on} to put DNA barcode on"
                    )
                if scaffold is None and self.dna_barcode_on == "scaffold":
                    raise LibraryBuildError("no scaffold to attach DNA barcode to")

    def __repr__(self):
        """Represent the library as its name"""
        return self.library_id

    def __eq__(self, other):
        """Check if two libraries are equal by their library_id"""
        if isinstance(other, DELibrary):
            return self.library_id == other.library_id
        return False

    @classmethod
    @accept_deli_data_name("libraries", "json")
    def load(cls, path: str) -> Self:
        """
        Load a library from the DELi data directory

        Notes
        -----
        This is decorated by `accept_deli_data`
        which makes this function actually take
          path_or_name: str

        `path_or_name` can be the full path to the file
        or it can be the name of the object to load

        See `Storing DEL info` in docs for more details


        Parameters
        ----------
        path: str
            path of the library to load

        Returns
        -------
        DELibrary
        """
        _cls = cls.read_json(path)
        _cls.loaded_from = path
        return _cls

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
        data = json.load(open(path))

        # load bb sets (needed for reaction setup)
        _observed_sets: list[tuple[int, BuildingBlockSet]] = list()
        for i, bb_data in enumerate(data["bb_sets"]):
            cycle = bb_data.get("cycle", None)
            if cycle is None:
                raise LibraryBuildError(
                    f"build block sets require a cycle number;" f"set at index {i} lacks a cycle"
                )
            bb_set_name = bb_data.get("bb_set_name", None)
            file_path = bb_data.get("bb_set_path", None)
            if file_path is None and bb_set_name is None:
                raise LibraryBuildError(
                    f"either 'bb_set_name' or 'bb_set_path' must be provided "
                    f"for a building block set;"
                    f"set in index {i} lacks a both"
                )
            if file_path is None:
                file_path = bb_set_name
            if bb_set_name is None:
                bb_set_name = os.path.basename(file_path).split(".")[0]
            _observed_sets.append((cycle, BuildingBlockSet.load(file_path, set_id=bb_set_name)))

        # check for right order of sets
        _bb_cycles = [_[0] for _ in _observed_sets]
        if _bb_cycles != list(range(1, len(_observed_sets) + 1)):
            raise LibraryBuildError(
                f"building block sets must be in consecutive ascending "
                f"order starting from 1 (1, 2, 3...); "
                f"observed order: '{_bb_cycles}'"
            )

        bb_sets: list[BuildingBlockSet] = [_[1] for _ in _observed_sets]
        bb_set_ids = set([bb_set.bb_set_id for bb_set in bb_sets] + ["scaffold"])

        if "reactions" in data.keys():
            reaction_workflow = ReactionWorkflow.load_from_json_list(data["reactions"], bb_set_ids)
        else:
            reaction_workflow = None

        return cls(
            library_id=data["id"],
            anchored=data.get("anchored", False),
            dna_barcode_on=data["dna_barcode_on"],
            barcode_schema=BarcodeSchema.from_dict(data["barcode_schema"]),
            bb_sets=bb_sets,
            library_reaction_workflow=reaction_workflow,
            scaffold=data.get("scaffold"),
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

    def iter_bb_barcode_sections_and_sets(
        self,
    ) -> Iterator[tuple[BuildingBlockBarcodeSection, BuildingBlockSet]]:
        """
        Iterate through building block sets and their respective barcode sections sections

        Yields
        ------
        tuple[BuildingBlockBarcodeSection, BuildingBlockSet]
        """
        for bb_section, bb_set in zip(self.barcode_schema.building_block_sections, self.bb_sets):
            yield bb_section, bb_set

    def enumerate_library_to_file(
        self, out_path: Union[str, Path], use_tqdm: bool = False
    ) -> None:
        """
        Enumerate the compound encoded in the DEL to a csv file

        Will auto generate DEL ids as <LIB_ID>-[<BB_ID>]
        for all building block sets in the lib
        e.g. L04-234-567-789 for a 3 cycle library with the id 'L04'

        Parameters
        ----------
        out_path: Union[str, Path]
            path to save csv file
        use_tqdm: bool, default False
            whether to use tqdm progress bar
        """
        bb_set_order = [bb_set.bb_set_id for bb_set in self.bb_sets]

        def del_id_func(bb_id_mapping: dict[str, BuildingBlock]) -> str:
            return f"{self.library_id}-" + "-".join(
                [bb_id_mapping[bb_set_id].bb_id for bb_set_id in bb_set_order]
            )

        if isinstance(self.enumerator, DELEnumerator):
            self.enumerator.enumerate_to_csv_file(out_path, del_id_func, use_tqdm=use_tqdm)
        else:
            raise RuntimeError(
                f"cannot enumerate library {self.library_id} without reaction information"
            )


class DELibraryPool:
    """base class for any class that holds a group of DEL libraries"""

    def __init__(self, libraries: List[DELibrary]):
        """
        Initialize a DELibrarySchemaGroup object

        Parameters
        ----------
        libraries: List[DELibrary]
            libraries to include in the library schema group
        """
        self.libraries = libraries
        self._library_map = {lib.library_id: lib for lib in self.libraries}

        self.pool_size = sum([lib.library_size for lib in self.libraries])

        ### VALIDATE ###
        _ids: list[str] = []
        _tags: list[str] = []
        for _library in self.libraries:
            # check id uniqueness
            if _library.library_id in _ids:
                raise LibraryBuildError(
                    f"multiple libraries share identical `library_id` '{_library.library_id}'"
                )
            else:
                _ids.append(_library.library_id)

            # check the bases uniqueness
            if _library.library_tag in _tags:
                _idx = _tags.index(_library.library_tag)
                raise LibraryBuildError(
                    f"library '{_library.library_id}' and library '{self.libraries[_idx]}' "
                    f"have identical DNA tag: '{_library.library_tag}'"
                )
            else:
                _tags.append(_library.library_tag)

    def __len__(self) -> int:
        """Return the number of libraries in the library pool"""
        return len(self.libraries)

    def __iter__(self) -> Iterator[DELibrary]:
        """Iterate through all libraries in the library pool"""
        return iter(self.libraries)

    def __eq__(self, other):
        """Check if two library pools are equal by their libraries"""
        if isinstance(other, DELibraryPool):
            return set(self.libraries) == set(other.libraries)
        return False

    def get_library(self, library_id: str) -> DELibrary:
        """
        return the library from the pool with the same ID

        Parameters
        ----------
        library_id: str
            id of the library to get

        Returns
        -------
        DELibrary

        Raises
        ------
        KeyError
            if library_id not in the pool
        """
        try:
            return self._library_map[library_id]
        except KeyError as e:
            raise KeyError(KeyError(f"cannot find library with id '{library_id}' in pool")) from e
