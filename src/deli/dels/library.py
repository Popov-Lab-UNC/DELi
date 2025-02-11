"""defines DEL library functions and classes"""

import json
from typing import Any, Iterator, List, Optional, Self, Tuple, Union, overload

from Levenshtein import distance

from deli.configure import DeliDataLoadable, accept_deli_data_name
from deli.utils.mol_utils import check_valid_smiles

from .barcode import BarcodeSchema, BarcodeSection
from .building_block import BuildingBlockSet
from .reaction import Reaction, ReactionStep, ReactionWorkflow


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
    """

    def __init__(
        self,
        library_id: str,
        library_dna_tag: str,
        barcode_schema: BarcodeSchema,
        bb_sets: List[BuildingBlockSet],
        library_reaction_workflow: ReactionWorkflow,
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
        self.library_tag = library_dna_tag
        self.barcode_schema = barcode_schema
        self.bb_sets = bb_sets
        self.library_reaction_workflow = library_reaction_workflow
        self.dna_barcode_on = dna_barcode_on
        self.scaffold = scaffold

        if len(self.bb_sets) < 2:
            raise LibraryBuildError(
                f"Library requires at least 2 cycle (bb_sets);" f" found {len(self.bb_sets)}"
            )
        self.num_cycles = len(self.bb_sets)

        # not sure this check makes sense anymore
        # if self.num_cycles + (0 if self.scaffold is None else -1) > len(self.reactions):
        #     raise LibraryBuildError(
        #         f"Library requires at least N-1 reactions for N cycles+scaffolds; "
        #         f"found {len(self.bb_sets)} cycles and {len(self.reactions)} reactions"
        #     )

        if self.num_cycles != barcode_schema.num_cycles:
            raise LibraryBuildError(
                f"Number of library cycles does not match barcode schema cycles; "
                f"got {self.num_cycles} and {barcode_schema.num_cycles}"
            )

        self.library_size = sum([len(bb_set) for bb_set in self.bb_sets])

        # handle the dna bases location
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
          deli_config: DeliConfig

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

        if "scaffold" not in data.keys():
            data["scaffold"] = None

        # load bb sets (needed for reaction setup)
        bb_sets: list[BuildingBlockSet] = [BuildingBlockSet.load(bb) for bb in data["bb_sets"]]
        bb_set_ids = [bb_set.bb_set_id for bb_set in bb_sets]
        # handle loading the reaction data

        rxn_steps: list[ReactionStep] = list()
        for rxn_data in data["reactions"]:
            step_id = rxn_data["step_id"]
            rxn_smarts = rxn_data["rxn_smarts"]
            reactant_ids = rxn_data["reactants"]

            # check if a reactant is static
            # static reactants are those that are valid SMILES and do not map to a bb_set id
            static_reactants: list[str] = list()
            variable_reactants: list[str] = list()
            for reactant_id in reactant_ids:
                if (reactant_id not in bb_set_ids) and (check_valid_smiles(reactant_id)):
                    static_reactants.append(reactant_id)
                else:
                    variable_reactants.append(reactant_id)
            rxn_steps.append(
                ReactionStep(
                    step_id=step_id,
                    reaction=Reaction(rxn_smarts),
                    variable_reactant=variable_reactants,
                    static_reactants=static_reactants,
                )
            )

        return cls(
            library_id=data["id"],
            library_dna_tag=data["library_tag"],
            dna_barcode_on=data["dna_barcode_on"],
            barcode_schema=BarcodeSchema.load(data["barcode_schema"]),
            bb_sets=bb_sets,
            library_reaction_workflow=ReactionWorkflow(rxn_steps),
            scaffold=data["scaffold"],
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

    def get_library_barcode_pattern(self) -> str:
        """Get the part of the barcode that is needed to decode the library"""
        return self.barcode_schema.full_barcode


class BaseDELibraryGroup:
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
        self._check_validity()

    def __len__(self) -> int:
        """Return number of libraries in the mega library"""
        return len(self.libraries)

    def __iter__(self) -> Iterator[DELibrary]:
        """Iterate through all libraries in the mega library"""
        return iter(self.libraries)

    @overload
    def __getitem__(self, index: int) -> DELibrary: ...

    @overload
    def __getitem__(self, index: slice) -> List[DELibrary]: ...

    def __getitem__(self, index: Union[int, slice]) -> Union[DELibrary, List[DELibrary]]:
        """Get the Library(s) at the passed index from the set"""
        return self.libraries[index]

    def __add__(self, other: Any) -> Self:
        """Add two library groups together (will fail if they overlap)"""
        if not isinstance(other, self.__class__):
            raise TypeError(
                f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'"
            )
        else:
            return self.__class__(libraries=self.libraries + other.libraries)

    def _check_validity(self):
        """Checks that there are no duplicate or conflicts in del mega library"""
        _ids = []
        _tags = []
        for _library in self.libraries:
            # check id uniqueness
            if _library.library_id in _ids:
                if _library.library_tag == self.libraries[_ids.index(_library.library_id)]:
                    raise LibraryBuildError(
                        f"identical library ids found in mega library: {_library.library_id}"
                    )
                else:
                    raise LibraryBuildError(
                        f"multiple indexes have the same id {_library.library_id}; "
                        f"index_ids must be unique"
                    )
            else:
                _ids.append(_library.library_id)

            # check the bases uniqueness
            if _library.library_tag in _tags:
                _idx = _tags.index(_library.library_tag)
                raise LibraryBuildError(
                    f"library {_library.library_id} and library {self.libraries[_idx]} "
                    f"have the same dna bases: {_library.library_tag}"
                )
            else:
                _tags.append(_library.library_tag)

    def has_library_with_name(self, lib_name: str) -> bool:
        """
        Returns True if the Group has a library with the passed ID/name

        Parameters
        ----------
        lib_name: str
            name/id of the library

        Returns
        -------
        bool
        """
        return any([lib.library_id == lib_name for lib in self.libraries])

    def get_library_with_name(self, lib_name: str) -> DELibrary:
        """
        Returns True if the Group has a library with the passed ID/name

        Parameters
        ----------
        lib_name: str
            name/id of the library

        Returns
        -------
        DELibrary
        """
        for lib in self.libraries:
            if lib.library_id == lib_name:
                return lib
        raise KeyError(f"cannot find library with name '{lib_name}' in LibraryGroup")


class DELibrarySchemaGroup(BaseDELibraryGroup):
    """
    A group of DELibraries where each group member has a calling compatible schema

    Notes
    -----
    Used to make sure calling of libraries can be done

    A schema group has all libraries share the same
    index, primer and library regions of their
    respective barcode schemas

    This way the caller can still
    figure out which library a read is in even if the barcodes
    have variable numbers of building blocks, different BB bases sizes
    or different UMI regions

    Raises
    ------
    LibraryBuildError:
        if SchemaGroup hold incompatible libraries (schemas do not match)
    """

    def __init__(self, libraries: List[DELibrary]):
        """
        Initialize a DELibrarySchemaGroup object

        Parameters
        ----------
        libraries: List[DELibrary]
            libraries to include in the library schema group
        """
        super().__init__(libraries=libraries)

        self.library_call_schema = self._generate_library_call_schema()
        self.requires_multistep_calling = self._requires_multistep_calling()

    def _check_validity(self):
        super()._check_validity()
        _lib_1 = self.libraries[0]
        for _library in self.libraries[1:]:
            if not _library.barcode_schema.is_library_compatible(_lib_1.barcode_schema):
                raise LibraryBuildError(
                    "all libraries in library schema group must have "
                    "the same index, primer and library regions"
                )

    def _generate_library_call_schema(self) -> BarcodeSchema:
        """
        Get the index (if present), primer and library dna bases for calling

        Returns
        -------
        str
        """
        _sections: List[Tuple[BarcodeSection, int]] = list()
        __index = self.libraries[0].barcode_schema.barcode_sections.get("index")
        if __index is None:
            pass
        else:
            _sections.append((__index, self.libraries[0].barcode_schema.get_position("index")))

        _sections.append(
            (
                self.libraries[0].barcode_schema.barcode_sections["library"],
                self.libraries[0].barcode_schema.get_position("library"),
            )
        )
        _sections.append(
            (
                self.libraries[0].barcode_schema.barcode_sections["primer"],
                self.libraries[0].barcode_schema.get_position("primer"),
            )
        )

        # sort by order
        _sections.sort(key=lambda x: x[1])

        return BarcodeSchema(
            schema_id="TMP_DO_NOT_SAVE",
            barcode_sections=[_[0] for _ in _sections],
            override_required=True,
        )

    def _requires_multistep_calling(self) -> bool:
        """
        Returns True if the SchemaGroup requires the multistep calling mode

        Notes
        -----
        Multistep calling is required when the barcodes of the
        group are not all identical.

        While a schema group requires the primer, index and library regions are
        all identical, the BB regions, UMI regions or any other can vary.
        If there are libraries in this group that have different schemas
        this prevents a single alignment to a single barcode schema.

        This function determines if this is the case and if the multistep
        calling mode is needed

        Also, this value is cached to improve performance

        Returns
        -------
        bool
        """
        # checks that all libraries and indexes can be merged into one experiment
        _lib_1 = self.libraries[0]
        for _library in self.libraries[1:]:
            if not _library.barcode_schema == _lib_1.barcode_schema:
                return True
        return False


class DELibraryGroup(BaseDELibraryGroup):
    """A set of many DELibraries"""

    # @classmethod
    # @accept_deli_data_name("libraries", "json")
    # def load(cls, path: str) -> Self:
    #     """
    #     Load a mega library from the DELi data directory
    #
    #     Notes
    #     -----
    #     This is decorated by `accept_deli_data`
    #     which makes this function actually take
    #       path_or_name: str
    #       deli_config: DeliConfig
    #
    #     `path_or_name` can be the full path to the file
    #     or it can be the name of the object to load
    #
    #     See `Storing DEL info` in docs for more details
    #
    #
    #     Parameters
    #     ----------
    #     path: str
    #         path of the mega library to load
    #
    #     Returns
    #     -------
    #     DELibraryGroup
    #     """
    #     _cls = cls.read_json(path)
    #     _cls.loaded_from = path
    #     return _cls
    #
    # @classmethod
    # def read_json(cls, path: str) -> Self:
    #     """
    #     Read a mega library from a json file
    #
    #     Parameters
    #     ----------
    #     path: str
    #         path to mega library json
    #
    #     Returns
    #     -------
    #     DELibraryGroup
    #     """
    #     data = json.load(open(path))
    #
    #     return cls(libraries=[DELibrary.from_dict(d) for d in data])

    def break_into_schema_groups(self) -> List[DELibrarySchemaGroup]:
        """
        Separates all libraries into groups based on schema compatability

        Notes
        -----
        This is done so that down stream calling can take place even if
        there barcode schemas are not all identical

        The groups mean that every library in that sub group has
        a barcode that shares the same 'index' (if used), 'primer'
        and 'library' regions. This way the caller can still
        figure out which library a read is in even if the barcodes
        have variable numbers of building blocks, different BB bases sizes
        or different UMI regions

        Returns
        -------
        List[DELLibrarySchemaGroup]
        """
        _groups: List[Tuple[BarcodeSchema, List[DELibrary]]] = list()

        for _library in self.libraries:
            _found_group = False
            for _schema, _libs in _groups:
                if _library.barcode_schema.is_library_compatible(_schema):
                    _libs.append(_library)
                    _found_group = True
                break
            if not _found_group:
                _groups.append((_library.barcode_schema, [_library]))

        return [DELibrarySchemaGroup(_libs[1]) for _libs in _groups]


def get_min_library_tag_distance(
    included_libraries: Optional[Union[list[DELibrary], BaseDELibraryGroup]],
) -> int:
    """
    Determine the minimum Levenshtein distance between all library tags

    Parameters
    ----------
    included_libraries: list[DELibrary] or BaseDELibraryGroup
        the libraries to use

    Returns
    -------
    min_distance: int
        the minimum Levenshtein distance between the passed library's bases
    """
    # if only one index just give it a constant threshold
    if included_libraries is None or len(included_libraries) <= 1:
        return 0
    _lib_dna_sequences = [i.library_tag for i in included_libraries]
    return min(
        [distance(s1, s2) for s1 in _lib_dna_sequences for s2 in _lib_dna_sequences if s1 != s2]
    )
