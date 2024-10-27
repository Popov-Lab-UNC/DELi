"""defines DEL library functions and classes"""

import json
from typing import Dict, Iterator, List, Optional, Self, Union, overload

from Levenshtein import distance

from deli.configure import accept_deli_data_name

from .barcode import BarcodeSchema
from .base import DeliDataLoadableMixin
from .building_block import BuildingBlockSet


class LibraryBuildError(Exception):
    """error raised when a library build fails"""

    pass


class Reaction:
    """struct to contain info on DEL reactions"""

    def __init__(self, cycle_id_1: str, cycle_id_2: str, reaction: str):
        self.cycle_id_1 = cycle_id_1
        self.bb_set_id_2 = cycle_id_2
        self.reaction = reaction


class DELibrary(DeliDataLoadableMixin):
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
    @accept_deli_data_name("libraries", ".json")
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
        return cls.read_json(path)

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
            library_id=lib_dict["id"],
            library_dna_tag=lib_dict["library_tag"],
            dna_barcode_on=lib_dict["dna_barcode_on"],
            barcode_schema=BarcodeSchema.load_from_json(lib_dict["barcode_schema"]),
            bb_sets=[
                BuildingBlockSet.load_from_csv(bb, no_smiles=True) for bb in lib_dict["bb_sets"]
            ],
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
        data = json.load(open(path))

        if "scaffold" not in data.keys():
            data["scaffold"] = None

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

    def get_library_barcode_pattern(self) -> str:
        """Get the part of the barcode that is needed to decode the library"""
        return self.barcode_schema.full_barcode


class DELibraryGroupMixin:
    """base class for any class that olds a group of DEL libraries"""

    libraries: List[DELibrary]

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


class DELibrarySchemaGroup(DELibraryGroupMixin):
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
    have variable numbers of building blocks, different BB tag sizes
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
        self.libraries = libraries
        self._check_validity()

        self.requires_multistep_calling = self._requires_multistep_calling()

    def _check_validity(self):
        _lib_1 = self.libraries[0]
        for _library in self.libraries[1:]:
            if not _library.barcode_schema.is_library_compatible(_lib_1.barcode_schema):
                raise LibraryBuildError(
                    "all libraries in library schema group must have "
                    "the same index, primer and library regions"
                )

    def get_library_call_tag(self) -> str:
        """
        Get the index (if present), primer and library dna tag for calling

        Returns
        -------
        str
        """
        __index = self.libraries[0].barcode_schema.barcode_sections.get("index")
        if __index is None:
            _index = ""
        else:
            _index = __index.section_tag

        _library = self.libraries[0].barcode_schema.barcode_sections["library"].section_tag
        _primer = self.libraries[0].barcode_schema.barcode_sections["primer"].section_tag

        return _index + _library + _primer

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
                return False
        return True


class DELibraryGroup(DELibraryGroupMixin):
    """A set of many DELibraries"""

    def __init__(self, libraries: List[DELibrary]):
        """
        Initialize a DELibraryGroup object

        Parameters
        ----------
        libraries: List[DELibrary]
            libraries to include in the mega library
        """
        self.libraries = libraries
        self._check_validity()

    @classmethod
    @accept_deli_data_name("libraries", ".json")
    def load(cls, path: str) -> Self:
        """
        Load a mega library from the DELi data directory

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
            path of the mega library to load

        Returns
        -------
        DELibraryGroup
        """
        return cls.read_json(path)

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
        DELibraryGroup
        """
        data = json.load(open(path))

        return cls(libraries=[DELibrary.from_dict(d) for d in data])

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

            # check the tag uniqueness
            if _library.library_tag in _tags:
                _idx = _tags.index(_library.library_tag)
                raise LibraryBuildError(
                    f"library {_library.library_id} and library {self.libraries[_idx]} "
                    f"have the same dna tag: {_library.library_tag}"
                )
            else:
                _tags.append(_library.library_tag)

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
        have variable numbers of building blocks, different BB tag sizes
        or different UMI regions

        Returns
        -------
        List[DELLibrarySchemaGroup]
        """
        _groups: Dict[BarcodeSchema, List[DELibrary]] = dict()

        for _library in self.libraries:
            _found_group = False
            for key, val in _groups.items():
                if _library.barcode_schema.is_library_compatible(key):
                    val.append(_library)
                    _found_group = True
                break
            if _found_group:
                _groups[_library.barcode_schema] = [_library]

        return [DELibrarySchemaGroup(_libs) for _libs in _groups.values()]


def get_min_library_tag_distance(
    included_libraries: Optional[Union[list[DELibrary], DELibraryGroupMixin]],
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
        the minimum Levenshtein distance between the passed library's tag
    """
    # if only one index just give it a constant threshold
    if included_libraries is None or len(included_libraries) <= 1:
        return 0
    _lib_dna_sequences = [i.library_tag for i in included_libraries]
    return min(
        [distance(s1, s2) for s1 in _lib_dna_sequences for s2 in _lib_dna_sequences if s1 != s2]
    )
