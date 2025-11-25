"""defines DEL library functions and classes"""

import json
import os
from collections.abc import Iterator
from functools import reduce
from operator import mul
from pathlib import Path
from typing import Any, Generic, Literal, Optional, Sequence, TypeVar, overload

from deli.configure import DeliDataLoadable, accept_deli_data_name
from deli.enumeration.enumerator import EnumeratedDELCompound, Enumerator
from deli.enumeration.reaction import ReactionTree
from deli.utils import to_smi

from .barcode import BarcodedMixin, BuildingBlockBarcodeSection, DELBarcodeSchema
from .base import Library
from .building_block import BuildingBlock, BuildingBlockSet, TaggedBuildingBlockSet
from .compound import DELCompound
from .tool_compounds import DopedToolCompound, ToolCompound


RESERVED_CONFIG_KEYS = [
    "linker",
    "truncated_linker",
    "scaffold",
    "reactions",
    "barcode_schema",
    "dna_barcode_on",
    "bb_sets",
]


def _parse_library_json(data: dict, load_dna: bool) -> dict[str, Any]:
    """
    Load from a JSON dict

    Parameters
    ----------
    data: dict
        data read in from JSON
    load_dna: bool
        whether to load information about DNA barcodes

    Returns
    -------
    dict[str, Any]
        arguments to construct the library object
    """
    # load bb sets and check for hamming decoding
    _observed_sets: list[tuple[int, BuildingBlockSet | TaggedBuildingBlockSet]] = list()
    for i, bb_data in enumerate(data["bb_sets"]):
        cycle = bb_data.get("cycle", None)
        if cycle is None:
            raise LibraryBuildError(f"build block sets require a cycle number;set at index {i} lacks a cycle")
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

        if bb_set_name in RESERVED_CONFIG_KEYS:
            raise LibraryBuildError(
                f"building block set name '{bb_set_name}' is a reserved keyword; please choose a different name"
            )

        _observed_sets.append(
            (
                cycle,
                BuildingBlockSet.load(file_path) if not load_dna else TaggedBuildingBlockSet.load(file_path),
            )
        )

    # check for right order of sets
    _bb_cycles = [_[0] for _ in _observed_sets]
    if _bb_cycles != list(range(1, len(_observed_sets) + 1)):
        raise LibraryBuildError(
            f"building block sets must be in consecutive ascending "
            f"order starting from 1 (1, 2, 3...); "
            f"observed order: '{_bb_cycles}'"
        )

    bb_sets: list[BuildingBlockSet] = [_[1] for _ in _observed_sets]
    bb_set_ids = set([bb_set.bb_set_id for bb_set in bb_sets])

    # check for doped compounds and validate
    tool_compounds: list[ToolCompound] = list()
    if "doped" in data.keys():
        for doped_bb in data["doped"]:
            doped_compound_id = doped_bb.get("compound_id", None)
            if doped_compound_id is None:
                raise LibraryBuildError("doped building blocks must have a 'doped_compound_id' field")
            smiles = doped_bb.get("smiles", None)
            if load_dna:
                _bb_tags: list[str] = list()
                for _bb_cycle in _bb_cycles:
                    bb_cycle_tag = doped_bb.get(f"bb{_bb_cycle}", None)
                    if bb_cycle_tag is not None:
                        raise LibraryBuildError(
                            f"doped compound {doped_compound_id} missing 'bb{_bb_cycle}' "
                            f"for {len(_bb_cycles)} cycle library"
                        )
                tool_compounds.append(
                    DopedToolCompound(
                        compound_id=doped_compound_id,
                        smiles=smiles,
                        bb_tags=tuple(_bb_tags),
                    )
                )
            else:
                tool_compounds.append(
                    ToolCompound(
                        compound_id=doped_compound_id,
                        smiles=smiles,
                    )
                )

    # check for scaffold and linker
    linker = data.get("linker", None)
    truncated_linker = data.get("truncated_linker", None)
    scaffold = data.get("scaffold", None)

    if "reactions" in data.keys():
        rxn_tree = ReactionTree.load_from_dict(
            data["reactions"],
            frozenset(bb_set_ids),
            static_comp_lookup={
                "linker": linker,
                "scaffold": scaffold,
                "truncated_linker": truncated_linker,
            },
            use_deli_data_dir=True,
        )
        enumerator = Enumerator(reaction_tree=rxn_tree, building_block_sets=bb_sets)
    else:
        enumerator = None

    res = {
        "bb_sets": bb_sets,
        "enumerator": enumerator,
        "scaffold": data.get("scaffold"),
    }

    if len(tool_compounds) > 0:
        res["tool_compounds"] = tool_compounds

    if load_dna:
        res.update(
            {
                "barcode_schema": DELBarcodeSchema.from_dict(data["barcode_schema"]),
                "dna_barcode_on": data.get("dna_barcode_on", None),
            }
        )
    return res


class LibraryBuildError(Exception):
    """error raised when a library build fails"""

    pass


class CombinatorialLibrary(Library, DeliDataLoadable):
    """
    A Library of compounds built combinatorially from building blocks

    This is primarily to provide support for a situation where everything
    about the library is known, but no decoding by DELi is needed
    (and the user might not have info about the DNA barcodes).

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
        bb_sets: Sequence[BuildingBlockSet],
        tool_compounds: Optional[Sequence[ToolCompound]] = None,
        enumerator: Optional[Enumerator] = None,
        scaffold: Optional[str] = None,
    ):
        """
        Initialize a DELibrary object

        Parameters
        ----------
        library_id : str
            name/id of the library
        bb_sets : Sequence[BuildingBlockSet]
            the sets of building-block used to build this library
            order in list should be the same as the order of synthesis
            must have length >= 2
        tool_compounds : optional, List[ToolCompound]
            list of tool compounds (doped compounds) in the library
        enumerator : optional, Enumerator
            The Enumerator used to build this library
        scaffold: optional, str
            SMILES of the scaffold
            if no scaffold in library should be `None`

        Raises
        ------
        LibraryBuildError
            the parameters passed to build the library are not compatible with each other
            error message will contain specific details of build issue
        """
        super().__init__(library_id=library_id)
        self.bb_sets = bb_sets
        self.scaffold = scaffold
        self.tool_compounds: Sequence[ToolCompound] = tool_compounds if tool_compounds else []

        self.library_size = (
            reduce(mul, [len(bb_set) for bb_set in self.bb_sets])
            if enumerator is None
            else enumerator.get_enumeration_size()
        )
        self.num_cycles = len(self.bb_sets)

        self._enumerator = enumerator
        self._can_enumerate = (self._enumerator is not None) and self.building_blocks_have_smi()

        ### VALIDATION ###

        if len(self.bb_sets) < 2:
            raise LibraryBuildError(f"Library requires at least 2 cycle (bb_sets); found {len(self.bb_sets)}")

    def __repr__(self):
        """Represent the library as its name"""
        return f"{self.__class__}({self.library_id})"

    def __hash__(self):
        """Hash the library by its id"""
        return hash(self.library_id)

    def __eq__(self, other):
        """Check equality of libraries by their id"""
        if isinstance(other, CombinatorialLibrary):
            return True
        return False

    @classmethod
    @accept_deli_data_name("libraries", "json")
    def load(cls, path: str) -> "CombinatorialLibrary":
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
        CombinatorialLibrary
        """
        library_id = Path(path).stem.replace(" ", "_")
        _cls = cls(library_id=library_id, **_parse_library_json(json.load(open(path)), load_dna=False))
        _cls.loaded_from = path
        return _cls

    @property
    def enumerator(self) -> Enumerator:
        """
        Return the enumerator for this library

        Returns
        -------
        Enumerator

        Raises
        ------
        RuntimeError
            if the library does not have an enumerator
        """
        if self._enumerator is None:
            raise RuntimeError(f"library {self.library_id} does not have an enumerator")
        if not self._can_enumerate:
            raise RuntimeError(f"library {self.library_id} cannot enumerate; missing building block SMILES")
        return self._enumerator

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

    def building_blocks_have_smi(self) -> bool:
        """
        Check if any building block sets in the library are missing SMILES

        Returns
        -------
        bool
            True if all building block sets have SMILES, False otherwise
        """
        return all([bb_set.has_smiles for bb_set in self.bb_sets])

    def can_enumerate(self) -> bool:
        """
        Return True if the library object can enumerate SMILES

        Returns
        -------
        bool
        """
        return self._can_enumerate

    @overload
    def enumerate(
        self, dropped_failed: Literal[True], fail_on_error: bool, use_tqdm: bool
    ) -> Iterator[EnumeratedDELCompound]: ...

    @overload
    def enumerate(
        self, dropped_failed: bool, fail_on_error: Literal[True], use_tqdm: bool
    ) -> Iterator[EnumeratedDELCompound]: ...

    @overload
    def enumerate(
        self, dropped_failed: Literal[False], fail_on_error: Literal[False], use_tqdm: bool
    ) -> Iterator[EnumeratedDELCompound | DELCompound]: ...

    @overload
    def enumerate(
        self, dropped_failed: bool, fail_on_error: bool, use_tqdm: bool
    ) -> Iterator[EnumeratedDELCompound | DELCompound]: ...

    def enumerate(
        self, dropped_failed: bool = False, fail_on_error: bool = False, use_tqdm: bool = False
    ) -> Iterator[EnumeratedDELCompound | DELCompound]:
        """
        Run reactions for all possible products of the building block set

        Notes
        -----
        The default behavior is to yield all compounds,
        whether they failed enumeration or not.
        This behavior can be controlled with the `dropped_failed` and `fail_on_error` parameters.

        Parameters
        ----------
        dropped_failed: bool, default = False
            skip yielding compounds that failed enumeration
            will only be utilized if `fail_on_error` is False
            otherwise an error will be raised if any compound fails enumeration
        fail_on_error: bool, default = False
            if True, will raise an error if ANY compound fails enumeration
            Unless you are confident all compounds will enumerate successfully,
            better to leave as False and instead yield the failed compounds without SMILES
        use_tqdm : bool
            turn on tqdm progress bar

        Yields
        ------
        EnumeratedDELCompound | DELCompound
            if dropped_failed is True, will yield only EnumeratedDELCompound

        Raises
        ------
        EnumerationRunError
            if `fail_on_error` is `True` and enumeration fails for any compound
            will raise an error with the failed compound's bb_id_map

        Warnings
        --------
        If `dropped_failed` is set to True and `fail_on_error` is also set to True,
        `fail_on_error` will take precedence
        """
        _enumerator = self.enumerator  # check it exists

        for building_blocks, enumerated_mol in _enumerator.enumerate(dropped_failed, fail_on_error, use_tqdm):
            if enumerated_mol is not None:
                enumerated_smile = to_smi(enumerated_mol)
                yield EnumeratedDELCompound(
                    library=self,
                    building_blocks=building_blocks,
                    smiles=enumerated_smile,
                    mol=enumerated_mol,
                )
            else:
                yield DELCompound(
                    library=self,
                    building_blocks=building_blocks,
                )

    def enumerate_by_bbs(self, bbs: list[BuildingBlock]) -> EnumeratedDELCompound:
        """
        Enumerate a single compound by providing the building blocks

        Parameters
        ----------
        bbs: list[BuildingBlock]
            list of building block ids to enumerate

        Returns
        -------
        EnumeratedDELCompound
            the enumerated compound with the provided building blocks

        Raises
        ------
        EnumerationRunError
            if enumeration fails
        """
        enumerated_mol = self.enumerator.enumerate_by_bbs(bbs)
        enumerated_smile = to_smi(enumerated_mol)
        return EnumeratedDELCompound(
            library=self,
            building_blocks=bbs,
            smiles=enumerated_smile,
            mol=enumerated_mol,
        )

    def enumerate_by_bb_ids(self, bb_ids: list[str]) -> EnumeratedDELCompound:
        """
        Enumerate a single compound by providing the building block ids

        Parameters
        ----------
        bb_ids: list[str]
            list of building block ids to enumerate

        Returns
        -------
        EnumeratedDELCompound
            the enumerated compound with the provided building blocks

        Raises
        ------
        EnumerationRunError
            if enumeration fails
        """
        return self.enumerate_by_bbs(
            [
                bb_set.get_bb_by_id(bb_id, fail_on_missing=True)
                for bb_set, bb_id in zip(self.bb_sets, bb_ids, strict=False)
            ]
        )

    def enumerate_to_file(
        self,
        out_path: str | Path,
        separator: str = "\t",
        dropped_failed: bool = False,
        fail_on_error: bool = False,
        use_tqdm: bool = False,
    ) -> None:
        r"""
        Enumerate the compound encoded in the DEL to a file

        Will have the following columns:
        - `DEL_ID`: the DEL id of the compound
        - `SMILES`: the SMILES of the compound
        - `LIB_ID`: the id of the library the compound belongs to
        - * for each building block set in the library * :
            - `BB<cycle#>_ID`: the ids of the building blocks used to build the compound
            - `BB<cycle#>_SMILES`: the smiles of the building block sets used to build the compound

        The default behavior is to yield all compounds,
        whether they failed enumeration or not.
        This behavior can be controlled with the `dropped_failed` and `fail_on_error` parameters.
        If failed compounds are included, the SMILES column will contain the text "NA"
        (NA is not a parsable SMILES)

        Parameters
        ----------
        out_path: Union[str, Path]
            path to save csv file
        separator: str, default = "\t"
            the char to be used as the separator in the output file
        dropped_failed: bool, default = False
            skip yielding compounds that failed enumeration
            will only be utilized if `fail_on_error` is False
            otherwise an error will be raised if any compound fails enumeration
        fail_on_error: bool, default = False
            if True, will raise an error if ANY compound fails enumeration
            Unless you are confident all compounds will enumerate successfully,
            better to leave as False and instead yield the failed compounds without SMILES
        use_tqdm : bool
            turn on tqdm progress bar

        Raises
        ------
        EnumerationRunError
            if `fail_on_error` is `True` and enumeration fails for any compound
            will raise an error with the failed compound's bb_id_map

        Warnings
        --------
        If `dropped_failed` is set to True and `fail_on_error` is also set to True,
        `fail_on_error` will take precedence
        """
        output_file = open(out_path, "w")

        # write the header
        _header = ["DEL_ID", "SMILES", "LIB_ID"]
        for cycle in range(1, len(self.bb_sets) + 1):
            _header.append(f"BB{cycle}_ID")
            _header.append(f"BB{cycle}_SMILES")
        output_file.write(separator.join(_header) + "\n")

        # enumerate and write each compound to the file
        for compound in self.enumerate(dropped_failed=dropped_failed, fail_on_error=fail_on_error, use_tqdm=use_tqdm):
            # get the row
            row: list[str] = [
                compound.compound_id,
                compound.smi if isinstance(compound, EnumeratedDELCompound) else "NA",
                compound.library.library_id,
            ]
            for bb in compound.building_blocks:
                row.append(bb.bb_id)
                row.append(bb.smi)

            output_file.write(separator.join(row) + "\n")


class DELibrary(CombinatorialLibrary, BarcodedMixin[DELBarcodeSchema]):
    """
    Contains information about a DNA-encoded library

    Attributes
    ----------
    num_cycles : int
        Number of building block cycles in the library
    library_size : int
        Size of the enumerated library
    library_tag : str
        the DNA tag mapped to this library
    """

    def __init__(
        self,
        library_id: str,
        barcode_schema: DELBarcodeSchema,
        bb_sets: Sequence[TaggedBuildingBlockSet],
        tool_compounds: Optional[Sequence[DopedToolCompound]] = None,
        enumerator: Optional[Enumerator] = None,
        dna_barcode_on: Optional[str] = None,
        linker: Optional[str] = None,
        truncated_linker: Optional[str] = None,
        scaffold: Optional[str] = None,
    ):
        """
        Initialize a DELibrary object

        Parameters
        ----------
        library_id : str
            name/id of the library
        barcode_schema : DELBarcodeSchema
            The barcode schema defining how the barcodes are designed
        bb_sets : Sequence[TaggedBuildingBlockSet]
            the sets of building-block used to build this library
            order in list should be the same as order of synthesis
            must have length >= 2
        tool_compounds : optional, Sequence[DopedToolCompound] = None
            list of tool compounds (doped compounds) in the library
        enumerator : ReactionTree
            The enumerator used to build this library
        dna_barcode_on: optional, str
            the id of the bb_set that is linked to the DNA bases
            can be 'scaffold' if DNA bases is linked to the scaffold
        linker: optional, str
            SMILES of the full linker
        truncated_linker: optional, str
            SMILES of the truncated linker; usually [C13]
        scaffold: optional, str
            SMILES of the scaffold
            if no scaffold in library should be `None`

        Raises
        ------
        LibraryBuildError
            the parameters passed to build the library are not compatible with each other
            error message will contain specific details of build issue
        """
        super().__init__(
            library_id=library_id,
            bb_sets=bb_sets,
            tool_compounds=tool_compounds,
            enumerator=enumerator,
            scaffold=scaffold,
        )
        # type hint again
        self.tool_compounds: Sequence[DopedToolCompound] = tool_compounds if tool_compounds else []
        self.bb_sets: Sequence[TaggedBuildingBlockSet] = bb_sets  # for type checker

        self.barcode_schema = barcode_schema
        self.dna_barcode_on = dna_barcode_on
        self.linker = linker
        self.truncated_linker = truncated_linker

        self.library_tag = self.barcode_schema.library_section.get_dna_sequence()

        ### VALIDATION ###
        # num bb sets matches barcode schema
        if self.num_cycles != barcode_schema.get_num_building_block_sections():
            raise LibraryBuildError(
                f"Number of library cycles does not match observed_barcode schema cycles; "
                f"got {self.num_cycles} and {barcode_schema.get_num_building_block_sections()}"
            )

        # Validate dna_barcode_on is attachd to an existing cycle or scaffold
        if isinstance(self.dna_barcode_on, str):
            if self.dna_barcode_on not in [bb_set.bb_set_id for bb_set in self.bb_sets]:
                if scaffold is not None and self.dna_barcode_on != "scaffold":
                    raise LibraryBuildError(f"cannot find cycle {self.dna_barcode_on} to put DNA barcode on")
                if scaffold is None and self.dna_barcode_on == "scaffold":
                    raise LibraryBuildError("no scaffold to attach DNA observed_barcode to")

        # Validate that all doped compounds have valid bb tags
        # and do not overlap with existing building block tags
        for tool_compound in self.tool_compounds:
            if len(tool_compound.bb_tags) != self.num_cycles:
                raise LibraryBuildError(
                    f"doped compound {tool_compound.compound_id} has {len(tool_compound.bb_tags)} "
                    f"bb tags; expected {self.num_cycles} for library {self.library_id}"
                )
            for i, (bb_tag, bb_set) in enumerate(zip(tool_compound.bb_tags, self.bb_sets, strict=True)):
                # if not `None` means tag exists already
                if bb_set.search_tags(bb_tag, fail_on_missing=False) is not None:
                    raise LibraryBuildError(
                        f"doped compound {tool_compound.compound_id} in library {self.library_id} "
                        f"has bb cycle {i} tag '{bb_tag}' which is found cycle {i} building block "
                        f"set {bb_set.bb_set_id}"
                    )

    @classmethod
    @accept_deli_data_name("libraries", "json")
    def load(cls, path: str) -> "DELibrary":
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
        library_id = Path(path).stem.replace(" ", "_")
        _cls = cls(library_id=library_id, **_parse_library_json(json.load(open(path)), load_dna=True))
        _cls.loaded_from = path
        return _cls

    def iter_bb_barcode_sections_and_sets(
        self,
    ) -> Iterator[tuple[BuildingBlockBarcodeSection, TaggedBuildingBlockSet]]:
        """
        Iterate through building block sets and their respective barcode sections

        Yields
        ------
        tuple[BuildingBlockBarcodeSection, BuildingBlockSet]
        """
        for bb_section, bb_set in zip(self.barcode_schema.building_block_sections, self.bb_sets, strict=False):
            yield bb_section, bb_set


LibType = TypeVar("LibType", bound=CombinatorialLibrary)


class LibraryCollection(Generic[LibType]):
    """
    base class for any class that holds a group of DEL libraries
    """

    def __init__(self, libraries: Sequence[LibType]):
        """
        Initialize a DELibrarySchemaGroup object

        Parameters
        ----------
        libraries: List[CombinatorialLibrary]
            libraries to include in the library schema group
        """
        self.libraries: Sequence[LibType] = libraries
        self._library_map = {lib.library_id: lib for lib in self.libraries}

        self.collection_size = sum([lib.library_size for lib in self.libraries])

        ### VALIDATE ###
        _ids: list[str] = []
        for _library in self.libraries:
            # check id uniqueness
            if _library.library_id in _ids:
                raise LibraryBuildError(f"multiple libraries share identical `library_id` '{_library.library_id}'")
            else:
                _ids.append(_library.library_id)

    def __len__(self) -> int:
        """Return the number of libraries in the library collection"""
        return len(self.libraries)

    def __iter__(self) -> Iterator[LibType]:
        """Iterate through all libraries in the library collection"""
        return iter(self.libraries)

    def get_library(self, library_id: str) -> LibType:
        """
        Return the library from the collection with the same ID

        Parameters
        ----------
        library_id: str
            id of the library to get

        Returns
        -------
        CombinatorialLibrary

        Raises
        ------
        KeyError
            if `library_id` not in the collection
        """
        try:
            return self._library_map[library_id]
        except KeyError as e:
            raise KeyError(KeyError(f"cannot find library with id '{library_id}' in collection")) from e

    def all_libs_can_enumerate(self) -> bool:
        """
        Check if all libraries in the collection have valid DEL enumerators

        Returns
        -------
        bool
        """
        return all([lib.can_enumerate() for lib in self.libraries])

    def all_libs_have_building_block_smiles(self) -> bool:
        """
        Check if all libraries in the collection have building block SMILES

        Returns
        -------
        bool
        """
        return all([lib.building_blocks_have_smi() for lib in self.libraries])

    def max_cycle_size(self) -> int:
        """
        Get the maximum cycle size of all libraries in the collection

        Returns
        -------
        int
            The maximum cycle size
        """
        return max([lib.num_cycles for lib in self.libraries])


class DELibraryCollection(LibraryCollection[DELibrary]):
    """base class for any class that holds a group of DEL libraries"""

    def __init__(self, libraries: Sequence[DELibrary]):
        """
        Initialize a DELibrarySchemaGroup object

        Parameters
        ----------
        libraries: List[DELibrary]
            libraries to include in the library schema group
        """
        super().__init__(libraries)
        self.libraries: Sequence[DELibrary] = libraries  # for type checker

        ### VALIDATE ###
        _tags: list[str] = []
        for _library in self.libraries:
            # check the bases uniqueness
            if _library.library_tag in _tags:
                _idx = _tags.index(_library.library_tag)
                raise LibraryBuildError(
                    f"library '{_library.library_id}' and library '{self.libraries[_idx]}' "
                    f"have identical DNA tag: '{_library.library_tag}'"
                )
            else:
                _tags.append(_library.library_tag)
