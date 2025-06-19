"""defines DEL library functions and classes"""

import json
import os
import warnings
from collections.abc import Iterator
from copy import deepcopy
from functools import reduce
from itertools import product as iter_product
from operator import mul
from pathlib import Path
from typing import Any, Generic, Literal, Optional, Sequence, TypeVar, Union, overload

from tqdm import tqdm

from deli.configure import DeliDataLoadable, accept_deli_data_name
from deli.utils import to_smi

from .barcode import BarcodeSchema, BuildingBlockBarcodeSection
from .building_block import BuildingBlock, BuildingBlockSet, TaggedBuildingBlockSet
from .compound import DELCompound, EnumeratedDELCompound
from .reaction import ReactionError, ReactionVial, ReactionWorkflow


class LibraryBuildError(Exception):
    """error raised when a library build fails"""

    pass


class EnumerationRunError(Exception):
    """error raised when a enumerator run fails"""

    pass


BBSetType = TypeVar("BBSetType", bound=BuildingBlockSet)


class Library(DeliDataLoadable):
    """
    A Library is a DEL that lacks information about the DNA barcodes

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
        reaction_workflow: Optional[ReactionWorkflow] = None,
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
            order in list should be order of synthesis
            must have length >= 2
        reaction_workflow : ReactionWorkflow
            The reaction workflow/schema used to build this library
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
        self.bb_sets = bb_sets
        self.scaffold = scaffold

        self.library_size = reduce(mul, [len(bb_set) for bb_set in self.bb_sets])
        self.num_cycles = len(self.bb_sets)

        self._reaction_workflow = reaction_workflow
        if isinstance(reaction_workflow, ReactionWorkflow) and self.building_blocks_have_smi():
            self._valid_reaction_workflow: ReactionWorkflow = reaction_workflow

        ### VALIDATION ###

        if len(self.bb_sets) < 2:
            raise LibraryBuildError(
                f"Library requires at least 2 cycle (bb_sets); found {len(self.bb_sets)}"
            )

    def __repr__(self):
        """Represent the library as its name"""
        return self.library_id

    @classmethod
    @accept_deli_data_name("libraries", "json")
    def load(cls, path: str) -> "Library":
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
        Library
        """
        _cls = cls(**cls.read_json(path))
        _cls.loaded_from = path
        return _cls

    @classmethod
    def read_json(cls, path: str) -> dict[str, Any]:
        """
        Load a DEL from a json file

        Parameters
        ----------
        path: str
            path to file with DEL json

        Returns
        -------
        dict[str, Any]
            arguments to construct the library object
        """
        data = json.load(open(path))

        # load bb sets and check for hamming decoding
        _observed_sets: list[tuple[int, BuildingBlockSet]] = list()
        for i, bb_data in enumerate(data["bb_sets"]):
            cycle = bb_data.get("cycle", None)
            if cycle is None:
                raise LibraryBuildError(
                    f"build block sets require a cycle number;set at index {i} lacks a cycle"
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
            _observed_sets.append((cycle, BuildingBlockSet.load(file_path)))

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

        # get library id from path
        library_id = Path(path).stem.replace(" ", "_")

        return {
            "library_id": library_id,
            "bb_sets": bb_sets,
            "reaction_workflow": reaction_workflow,
            "scaffold": data.get("scaffold"),
        }

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
        return any([bb_set.has_smiles for bb_set in self.bb_sets])

    def _get_reaction_workflow(self) -> ReactionWorkflow:
        """
        Check and get the library has a valid reaction workflow

        Will fail if library cannot enumerate using the reaction workflow

        Raises
        ------
        EnumerationRunError
            if the library lacks enough information to run the reaction workflow
        """
        if hasattr(self, "_valid_reaction_workflow"):
            return self._valid_reaction_workflow
        else:
            if self._reaction_workflow is None:
                raise EnumerationRunError(
                    f"cannot enumerate library {self.library_id}; lacking reaction information"
                )
            else:
                for bb_set in self.bb_sets:
                    if not bb_set.has_smiles:
                        raise EnumerationRunError(
                            f"cannot enumerate library {self.library_id}; "
                            f"building block set {bb_set.bb_set_id} lacks SMILES"
                        )
                # this error is unreachable, but here to help with clarity,
                # as the above loop must raise an enumeration error
                raise EnumerationRunError(
                    f"cannot enumerate library {self.library_id}; building blocks missing SMILES"
                )

    def can_enumerate(self) -> bool:
        """
        Return True if the library object can enumerate SMILES

        Returns
        -------
        bool
        """
        return hasattr(self, "_valid_reaction_workflow")

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
            if dropped_failed is False, will yield EnumeratedDELCompound

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
        _reaction_workflow = self._get_reaction_workflow()

        if dropped_failed and fail_on_error:
            warnings.warn(
                "'dropped_failed' is set to True, "
                "but 'fail_on_error' is also set to True; "
                "ignoring 'dropped_failed' and will raise "
                "an error if any compound fails enumeration",
                stacklevel=1,
            )

        bb_set_dict = {bb_set.bb_set_id: bb_set for bb_set in self.bb_sets}

        # get all possible combinations of building blocks
        bb_combos: iter_product[tuple[tuple[str, BuildingBlock], ...]] = iter_product(
            *[[(key, bb) for bb in val.building_blocks] for key, val in bb_set_dict.items()]
        )

        return self._enumerate(
            reaction_workflow=_reaction_workflow,
            bb_combos=bb_combos,
            dropped_failed=dropped_failed,
            fail_on_error=fail_on_error,
            use_tqdm=use_tqdm,
        )

    @overload
    def _enumerate(
        self,
        reaction_workflow: ReactionWorkflow,
        bb_combos,
        dropped_failed: Literal[True],
        fail_on_error: bool,
        use_tqdm: bool,
    ) -> Iterator[EnumeratedDELCompound]: ...

    @overload
    def _enumerate(
        self,
        reaction_workflow: ReactionWorkflow,
        bb_combos,
        dropped_failed: bool,
        fail_on_error: Literal[True],
        use_tqdm: bool,
    ) -> Iterator[EnumeratedDELCompound]: ...

    @overload
    def _enumerate(
        self,
        reaction_workflow: ReactionWorkflow,
        bb_combos,
        dropped_failed: Literal[False],
        fail_on_error: Literal[False],
        use_tqdm: bool,
    ) -> Iterator[EnumeratedDELCompound | DELCompound]: ...

    @overload
    def _enumerate(
        self,
        reaction_workflow: ReactionWorkflow,
        bb_combos,
        dropped_failed: bool,
        fail_on_error: bool,
        use_tqdm: bool,
    ) -> Iterator[EnumeratedDELCompound | DELCompound]: ...

    def _enumerate(
        self,
        reaction_workflow: ReactionWorkflow,
        bb_combos,
        dropped_failed: bool,
        fail_on_error: bool,
        use_tqdm: bool,
    ) -> Iterator[EnumeratedDELCompound | DELCompound]:
        """Run the enumeration after check that enumeration is possible"""
        for bb_combo in tqdm(bb_combos, total=self.library_size, disable=not use_tqdm):
            # map to cycle number (required for the ReactionVial)
            bb_id_map = {bb[0]: bb[1] for bb in bb_combo}
            try:
                enumerated_mol = reaction_workflow.run_workflow(
                    ReactionVial({bb[0]: deepcopy(bb[1].mol) for bb in bb_combo})
                )
                enumerated_smile = to_smi(enumerated_mol)
                yield EnumeratedDELCompound(
                    library=self,
                    building_blocks=list(bb_id_map.values()),
                    smiles=enumerated_smile,
                    mol=enumerated_mol,
                )
            except ReactionError as e:  # crash is fail_on_error is True
                if fail_on_error:
                    raise EnumerationRunError(
                        f"enumeration failed for compound from library {self.library_id} "
                        f"with building block: {bb_id_map}; "
                    ) from e
                elif dropped_failed:  # skip yield if dropped_failed is True
                    continue
                else:  # yield the compound with no SMILES otherwise
                    yield DELCompound(library=self, building_blocks=list(bb_id_map.values()))

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
        _reaction_workflow = self._get_reaction_workflow()

        bb_map = {bb_set.bb_set_id: bb for bb_set, bb in zip(self.bb_sets, bbs)}
        try:
            enumerated_mol = _reaction_workflow.run_workflow(
                ReactionVial({bb_cycle: deepcopy(bb.mol) for bb_cycle, bb in bb_map.items()})
            )
            enumerated_smile = to_smi(enumerated_mol)
            return EnumeratedDELCompound(
                library=self,
                building_blocks=list(bb_map.values()),
                smiles=enumerated_smile,
                mol=enumerated_mol,
            )
        except ReactionError as e:
            raise EnumerationRunError(
                f"enumeration failed for compound from library {self.library_id} "
                f"with building blocks: {[bb.bb_id for bb in bbs]}; "
            ) from e

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
                for bb_set, bb_id in zip(self.bb_sets, bb_ids)
            ]
        )

    def enumerate_to_file(
        self,
        out_path: Union[str, Path],
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
        for compound in self.enumerate(
            dropped_failed=dropped_failed, fail_on_error=fail_on_error, use_tqdm=use_tqdm
        ):
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


class DELibrary(Library):
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
        barcode_schema: BarcodeSchema,
        bb_sets: Sequence[TaggedBuildingBlockSet],
        reaction_workflow: Optional[ReactionWorkflow] = None,
        dna_barcode_on: Optional[str] = None,
        scaffold: Optional[str] = None,
    ):
        """
        Initialize a DELibrary object

        Parameters
        ----------
        library_id : str
            name/id of the library
        barcode_schema : BarcodeSchema
            The barcode schema defining how the barcodes are designed
        bb_sets : Sequence[TaggedBuildingBlockSet]
            the sets of building-block used to build this library
            order in list should be order of synthesis
            must have length >= 2
        reaction_workflow : ReactionWorkflow
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
        super().__init__(
            library_id=library_id,
            bb_sets=bb_sets,
            reaction_workflow=reaction_workflow,
            scaffold=scaffold,
        )
        self.bb_sets: Sequence[TaggedBuildingBlockSet] = bb_sets  # for type checker

        self.barcode_schema = barcode_schema
        self.dna_barcode_on = dna_barcode_on

        self.library_tag = self.barcode_schema.library_section.get_dna_sequence()

        ### VALIDATION ###
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
        _cls = cls(**cls.read_json(path))
        _cls.loaded_from = path
        return _cls

    @classmethod
    def read_json(cls, path: str) -> dict[str, Any]:
        """
        Load a DEL from a json file

        Parameters
        ----------
        path: str
            path to file with DEL json

        Returns
        -------
        dict[str, Any]
            argument to construct the DELibrary object
        """
        data = json.load(open(path))

        # load bb sets and check for hamming decoding
        _observed_sets: list[tuple[int, TaggedBuildingBlockSet]] = list()
        for i, bb_data in enumerate(data["bb_sets"]):
            cycle = bb_data.get("cycle", None)
            if cycle is None:
                raise LibraryBuildError(
                    f"build block sets require a cycle number;set at index {i} lacks a cycle"
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
            _observed_sets.append((cycle, TaggedBuildingBlockSet.load(file_path)))

        # check for right order of sets
        _bb_cycles = [_[0] for _ in _observed_sets]
        if _bb_cycles != list(range(1, len(_observed_sets) + 1)):
            raise LibraryBuildError(
                f"building block sets must be in consecutive ascending "
                f"order starting from 1 (1, 2, 3...); "
                f"observed order: '{_bb_cycles}'"
            )

        bb_sets: list[TaggedBuildingBlockSet] = [_[1] for _ in _observed_sets]
        bb_set_ids = set([bb_set.bb_set_id for bb_set in bb_sets] + ["scaffold"])

        if "reactions" in data.keys():
            reaction_workflow = ReactionWorkflow.load_from_json_list(data["reactions"], bb_set_ids)
        else:
            reaction_workflow = None

        # get library id from path
        library_id = Path(path).stem.replace(" ", "_")

        return {
            "library_id": library_id,
            "bb_sets": bb_sets,
            "reaction_workflow": reaction_workflow,
            "scaffold": data.get("scaffold"),
            "barcode_schema": BarcodeSchema.from_dict(data["barcode_schema"]),
            "dna_barcode_on": data["dna_barcode_on"],
        }

    def iter_bb_barcode_sections_and_sets(
        self,
    ) -> Iterator[tuple[BuildingBlockBarcodeSection, TaggedBuildingBlockSet]]:
        """
        Iterate through building block sets and their respective barcode sections sections

        Yields
        ------
        tuple[BuildingBlockBarcodeSection, BuildingBlockSet]
        """
        for bb_section, bb_set in zip(self.barcode_schema.building_block_sections, self.bb_sets):
            yield bb_section, bb_set


LibType = TypeVar("LibType", bound=Library)


class LibraryCollection(Generic[LibType]):
    """
    base class for any class that holds a group of DEL libraries
    """

    def __init__(self, libraries: Sequence[LibType]):
        """
        Initialize a DELibrarySchemaGroup object

        Parameters
        ----------
        libraries: List[Library]
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
                raise LibraryBuildError(
                    f"multiple libraries share identical `library_id` '{_library.library_id}'"
                )
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
        Library

        Raises
        ------
        KeyError
            if `library_id` not in the collection
        """
        try:
            return self._library_map[library_id]
        except KeyError as e:
            raise KeyError(
                KeyError(f"cannot find library with id '{library_id}' in collection")
            ) from e

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
