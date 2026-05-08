"""define building block classes"""
import warnings
from collections import defaultdict
from csv import DictReader
from io import StringIO
from typing import Generic, Literal, Optional, Sequence, TypeVar, overload

from deli.config import get_deli_config
from deli.utils.mol_utils import SmilesMixin, check_valid_smiles

from ._base import DecodableObject, DecodableObjectSet, TaggedDecodableObject, TaggedDecodableObjectSet
from .parsing import Parsable


def _parse_bb_row_dict(
    row_dict, load_smiles: bool = True, validate_smiles: bool = True
) -> tuple[str | None, str | None, str | None, str | None, bool]:
    """Helper func to parse a row dict from a building block file."""
    config = get_deli_config()

    _id = row_dict.get(config.bb_id_column)
    _tag = row_dict.get(config.bb_tag_column)
    _smiles = row_dict.get(config.bb_smiles_column)
    _subset_id = row_dict.get(config.bb_subset_id_column)
    _is_null_str = row_dict.get(config.bb_is_null_column)

    # convert to a boolean if the is_null column is present, else None
    _is_null_bool = _is_null_str.lower() in ["true", "y", "yes", "1"] if _is_null_str is not None else None

    if _smiles and validate_smiles:
        _valid_smiles = check_valid_smiles(_smiles)
        if _is_null_bool is not None:
            if _valid_smiles == _is_null_bool:
                if _valid_smiles:
                    msg = (
                        f"building block {_id} is explicitly marked "
                        f"as null but appears to have a valid SMILES '{_smiles}'"
                    )
                else:
                    msg = (
                        f"building block {_id} is explicitly marked as not-null "
                        f"but appears to have an invalid SMILES '{_smiles}'"
                    )
                warnings.warn(msg, stacklevel=2)
        else:
            _is_null_bool = not _valid_smiles
    else:
        _is_null_bool = _is_null_bool or False  # default to false if not specified

    return _id, _tag, _smiles, _subset_id, _is_null_bool


@overload
def _parse_building_block_file_contents(
    contents: str,
    include_tags: Literal[False],
    delimiter: str = ",",
    load_smiles: bool = True,
    validate_chemicals: bool = True
) -> "list[BuildingBlock]": ...

@overload
def _parse_building_block_file_contents(
    contents: str,
    include_tags: Literal[True],
    delimiter: str = ",",
    load_smiles: bool = True,
    validate_chemicals: bool = True
) -> "list[TaggedBuildingBlock]": ...

@overload
def _parse_building_block_file_contents(
    contents: str,
    include_tags: bool,
    delimiter: str = ",",
    load_smiles: bool = True,
    validate_chemicals: bool = True
) -> "list[BuildingBlock] | list[TaggedBuildingBlock]": ...

def _parse_building_block_file_contents(
    contents: str,
    include_tags: bool,
    delimiter: str = ",",
    load_smiles: bool = True,
    validate_chemicals: bool = True
) -> "list[BuildingBlock] | list[TaggedBuildingBlock]":
    """
    Parse the contents of a building block set file and return a BuildingBlockSet or TaggedBuildingBlockSet

    Will determine whether to return a BuildingBlockSet or TaggedBuildingBlockSet based on presence of 'tag' column
    in the building block file. The behavior can be control explicitly using the `include_tags` parameter.

    Parameters
    ----------
    contents: str
        contents of the building block set file
    include_tags: bool
        whether to include tags in the returned building block objects. If `True`,
        will return a TaggedBuildingBlockSet, else will return a BuildingBlockSet.
        Note that if `include_tags` is `False` but the file contains a 'tag' column,
        the tags will be ignored.
    delimiter: str, default = ","
        delimiter used in the building block set file (default is comma for csv files)
    load_smiles: bool, default = True
        whether to attempt to load SMILES strings for the building blocks. If `True`,
        will attempt to load SMILES from the csv file. Will *not* fail if SMILES are
        missing
    validate_chemicals: bool, default = True
        whether to validate the chemicals in the building block set file. If `True`,
        will check if the chemicals are valid. This can be computationally expensive,
        so it would be best to set this to `False` unless you want the extra validation.

    Raises
    ------
    BuildingBlockSetError
        If the file is missing required columns or has invalid formatting
    """
    if load_smiles is False and validate_chemicals is True:
       warnings.warn(
        "cannot validate chemicals if SMILES loading is disabled; skipping chemical validation", stacklevel=2
    )

    config = get_deli_config()
    reader = DictReader(StringIO(contents), delimiter=delimiter)

    if reader.fieldnames is None:
        raise BuildingBlockSetError("building block set file is missing header row")
    reader.fieldnames = [field.strip() for field in reader.fieldnames]  # strip whitespace from header fields

    if config.bb_id_column not in reader.fieldnames:
        warnings.warn(
            f"building block set file is missing '{config.bb_id_column}' column; "
            f"auto-generating ids based on row index",
            stacklevel=2
        )
    if include_tags and config.bb_tag_column not in reader.fieldnames:
        raise BuildingBlockSetError(
            f"building block set file is missing required column '{config.bb_tag_column}' "
            f"for tagged building block set"
        )

    _building_blocks: list[BuildingBlock] = list()
    _building_block_map: dict[str, BuildingBlock] = dict()

    # loop through all blocks and condense into unique building blocks and check for duplicates
    for i, row in enumerate(reader):
        bb_id, bb_tag, bb_smiles, bb_subset_id, bb_is_null = _parse_bb_row_dict(row, validate_smiles=validate_chemicals)
        if bb_id is None:
            bb_id = str(i + 1)  # autogenerate an id if not provided, using the row index to ensure uniqueness
        if include_tags and bb_tag is None:
            raise BuildingBlockSetError(
                f"building block set file is missing required column '{config.bb_tag_column}' "
                f"for tagged building block set in row {i+2}"
            )

        # this building block has been seen before
        if bb_id in _building_block_map.keys():
            # check if smiles is the same
            if load_smiles:
                if (bb_smiles is not None) and (bb_smiles != _building_block_map[bb_id].smi):
                    if validate_chemicals:
                        raise BuildingBlockSetError(
                            f"duplicate building block id '{bb_id}' with conflicting SMILES in row {i+2}: "
                            f"saw '{bb_smiles}', expected '{_building_block_map[bb_id].smi}' "
                        )
                    else:
                        warnings.warn(
                            f"duplicate building block id '{bb_id}' with non-identical SMILES in row {i+2}. "
                            f"Without SMILES validation enabled SMILES could encode same chemical. "
                            f"Will treat blocks as potetnial duplicates ignoring SMILES conflict."
                            f"saw '{bb_smiles}', expected '{_building_block_map[bb_id].smi}'\n",
                            stacklevel=2,
                        )
            # check if subset_id is the same
            if (bb_subset_id is not None) and (bb_subset_id != _building_block_map[bb_id].subset_id):
                raise BuildingBlockSetError(
                    f"duplicate building block id '{bb_id}' with conflicting subset_ids in row {i+2}: "
                    f"'saw {bb_subset_id}', expected '{_building_block_map[bb_id].subset_id}' "
                )  # skip duplicate with same smiles
            if (bb_is_null != _building_block_map[bb_id].is_null()):
                raise BuildingBlockSetError(
                    f"duplicate building block id '{bb_id}' with conflicting null status in row {i+2}: "
                    f"'{bb_is_null}', expected '{_building_block_map[bb_id].is_null()}' "
                )
            # if these are tagged register the new tag if differnt otherwise ignore the duplicate
            if include_tags and (bb_tag not in _building_block_map[bb_id].tags):  # ty:ignore[unresolved-attribute]  # if include tag will have tags attribute
                _building_block_map[bb_id].tags.append(bb_tag)  # ty:ignore[unresolved-attribute]
            else:
                warnings.warn(
                    f"duplicate building block id '{bb_id}' detected in row {i+2} "
                    f"with non-conflicting information; ignoring",
                    category=UserWarning,
                    stacklevel=2,
                )
        else:
            if bb_is_null:
                _building_blocks.append(
                    TaggedNullBuildingBlock(bb_id=bb_id, smiles=bb_smiles, tag=[bb_tag], subset_id=bb_subset_id)  # ty:ignore[invalid-argument-type]
                    if include_tags else NullBuildingBlock(bb_id=bb_id, smiles=bb_smiles, subset_id=bb_subset_id)
                )
            else:
                _building_blocks.append(
                    TaggedBuildingBlock(bb_id=bb_id, smiles=bb_smiles, tags=[bb_tag], subset_id=bb_subset_id)  # ty:ignore[invalid-argument-type]
                    if include_tags else BuildingBlock(bb_id=bb_id, smiles=bb_smiles, subset_id=bb_subset_id)
                )
    return _building_blocks


class BuildingBlockSetError(Exception):
    """raised when there is an issue with building block set definition"""

    pass


class BuildingBlock(SmilesMixin, DecodableObject):
    """
    Define building block that has a chemical structure associated with it

    Parameters
    ----------
    bb_id: str
        building block id
    smiles: Optional[str], default = None
        SMILES of the building block
    subset_id: str | None
        if building block belongs to a subset, the subset id
        else None
    """

    def __init__(self, bb_id: str, smiles: Optional[str] = None, subset_id: Optional[str] = None, **kwargs):
        """Initialize the object"""
        super().__init__(obj_id=bb_id, smiles=smiles, **kwargs)
        self.subset_id = str(subset_id) if subset_id else None  # stringify

        get_deli_config().check_if_bb_id_is_reserved(bb_id)

    @property
    def bb_id(self) -> str:
        """Alias for obj_id to help with readablity"""
        return self._obj_id

    def __eq__(self, other):
        """Two building blocks are equal if their ID is equal"""
        return self.bb_id == other.bb_id if isinstance(other, BuildingBlock) else False

    def __str__(self):
        """Cast to str as the bb_id"""
        return self.bb_id

    def __repr__(self):
        """Return a string representation of the object."""
        return f"BuildingBlock(bb_id='{self.bb_id}')"

    @classmethod
    def from_dict(cls, data: dict[str, str]):
        """
        Initialize the building block from a dictionary.

        Optional fields can be missing from the dictionary.

        Parameters
        ----------
        data: dict[str, str]
            Dictionary of building block properties

        Returns
        -------
        BuildingBlock
        """
        return cls(data["id"], smiles=data.get("smiles", None), subset_id=data.get("subset_id", None))

    def has_same_cpd_as_block(self, other: "BuildingBlock") -> bool:
        """
        Check if two building blocks are have the same chemical structure based on their SMILES

        Treats `None` and a unique chemical, so if two compounds both have `None` for a SMILES
        they will *not* be considered to have the same compound.

        Parameters
        ----------
        other: BuildingBlock
            the building block to compare to

        Returns
        -------
        bool
            True if the building blocks are identical in all attributes, else False
        """
        if not (self._smiles is None or other._smiles is None):
            return self.inchi_key == other.inchi_key
        return False

    def in_subset(self) -> bool:
        """
        Check if building block belongs to a subset

        Returns
        -------
        bool
        """
        return self.subset_id is not None

    def get_subset_id(self) -> str:
        """
        Get the subset id of the building block

        Returns
        -------
        str
            subset id of the building block

        Raises
        ------
        ValueError
            if the building block does not belong to a subset
        """
        if self.subset_id is None:
            raise RuntimeError(f"building block '{self.bb_id}' does not belong to a subset")
        return self.subset_id

    @staticmethod
    def is_null() -> bool:
        """
        Check if building block is a null building block

        Returns
        -------
        bool
        """
        return False

    @staticmethod
    def is_real() -> bool:
        """
        Check if building block is a real (non-fake) building block

        Returns
        -------
        bool
        """
        return True


class NullBuildingBlock(BuildingBlock):
    """
    Define a null building block

    Null building blocks are building blocks that are not involved in the combinatorial
    synthesis of a larger building block set.
    They can have SMILES or no SMILES, but they will not be considered in any
    combinatorial enumeration

    In practice, these are common in DELs where a building block position is skipped
    or capped during synthesis to enable QC.

    Notes
    -----
    These are nearly equivalent to BuildingBlock objects, but will
    always return `True` for `is_null()` method.

    Parameters
    ----------
    bb_id: str
        building block id
    smiles: Optional[str], default = None
        SMILES of the building block
    subset_id: Optional[str], default = None
        if building block belongs to a subset, the subset id
    """

    def __init__(self, bb_id: str, smiles: Optional[str] = None, subset_id: Optional[str] = None):
        """Initialize the null building block"""
        super().__init__(bb_id=bb_id, smiles=smiles, subset_id=subset_id)

    @staticmethod
    def is_null() -> bool:
        """
        Check if building block is a null building block

        Returns
        -------
        bool
        """
        return True


class TaggedBuildingBlock(BuildingBlock, TaggedDecodableObject):
    """
    Define building block that has a chemical structure and a DNA tag associated with it

    Attributes
    ----------
    tags: list[str]
        the dna tag(s) associated with this building block

    Parameters
    ----------
    bb_id: str
        building block id
    tag: str | list[str]
        the dna tag associated with this building block
    smiles: str, default = None
        SMILES of the building block
    subset_id: str, default = None
        if building block belongs to a subset, the subset id
    """

    def __init__(
        self,
        bb_id: str,
        tags: str | list[str],
        smiles: Optional[str] = None,
        subset_id: Optional[str] = None,
    ):
        """Initialize the tagged building block"""
        super().__init__(bb_id=bb_id, smiles=smiles, subset_id=subset_id, tags=tags)

    @classmethod
    def from_dict(cls, data: dict[str, str]):
        """
        Initialize the building block from a dictionary

        Parameters
        ----------
        data: dict
            Dictionary of building block properties

        Returns
        -------
        BuildingBlock
        """
        return cls(bb_id=data["id"], tags=data["tag"], smiles=data.get("smiles", None))


class TaggedNullBuildingBlock(TaggedBuildingBlock):
    """
    Define a null tagged building block

    Null building blocks are building blocks that are not involved in the combinatorial
    synthesis of a larger building block set.
    They can have SMILES or no SMILES, but they will not be considered in any
    combinatorial enumeration

    In practice, these are common in DELs where a building block position is skipped
    or capped during synthesis to enable QC.

    These are tagged to enable decoding from DEL selections.
    Unlike the chemical synthesis, tagged null building blocks will have their
    possible full tag sequences enumerated during decoding (since DNA ligation
    occurs even if no chemical synthesis occurs at that position).

    Notes
    -----
    These are nearly equivalent to BuildingBlock objects, but will
    always return `True` for `is_null()` method.

    Parameters
    ----------
    bb_id: str
        building block id
    tag: str | list[str]
        the dna tag associated with this building block
    smiles: str, default = None
        SMILES of the building block
    subset_id: str, default = None
        if building block belongs to a subset, the subset id
    """

    def __init__(
        self,
        bb_id: str,
        tag: str | list[str],
        smiles: Optional[str] = None,
        subset_id: Optional[str] = None,
    ):
        super().__init__(bb_id=bb_id, tags=tag, smiles=smiles, subset_id=subset_id)

    @staticmethod
    def is_null() -> bool:
        """
        Check if building block is a null building block

        Returns
        -------
        bool
        """
        return True


class TaggedFakeBuildingBlock(TaggedNullBuildingBlock):
    """
    Define a fake tagged building block

    Fake tagged building blocks are a special type of Null building block, one that is
    artificially created for decoding purposes only. Like Null blocks, they are not used
    for enumeration, but unlike Null blocks, they are not used during synthesis.

    This is useful for representing compounds within a DEL that follow the same barcode tag
    design, but actually map to compounds that were not synthesized within the DEL and added in
    after. For example, spiked-in controls or known binders added to the pool after synthesis.

    These only exist in a tagged form, since they are not part of the chemical synthesis.
    Thus using them outside a DNA decoding context is meaningless.

    Parameters
    ----------
    bb_id: str
        building block id
    tag: str | list[str]
        the dna tag associated with this building block
    """

    def __init__(self, bb_id: str, tag: str | list[str]):
        """Initialize the fake tagged building block"""
        super().__init__(bb_id=bb_id, tag=tag)

    @staticmethod
    def is_real() -> bool:
        """
        Check if building block is a real (non-fake) building block

        Returns
        -------
        bool
        """
        return False


T = TypeVar("T", bound=BuildingBlock)

class BuildingBlockSet(Parsable, DecodableObjectSet[T], Generic[T], sub_dir="building_blocks"):
    """
    Holds a set of building blocks

    Can hold both real and null building blocks.
    The only required restirction for a set of building blocks
    is that they have unqiue IDs.

    Chemical validation is optional. If enabled DELi will also check that:
    - all non-null and real building blocks have valid SMILES (raises exception)
    - all non-null and real building blocks have a unique chemical structure (raises warning)

    Attributes
    ----------
    has_smiles: bool
        True if all BuildingBlocks in this set have SMILES
        does not guarantee SMILES are valid

    Parameters
    ----------
    building_blocks: Sequence[BuildingBlock]
        list of building blocks in set
    has_unique_smiles: bool, default = False
        whether all building blocks in the set have unique SMILES.
        Will also fail if any (non-null) blocks are missing SMILES

    Notes
    -----
    If validating chemicals, Mol objects will be cached for each building block,
    which can inflate memory usage. You can clear the entire set with `clear_cached_mols()`.
    """

    def __init__(self, building_blocks: Sequence[T], has_unique_smiles: bool = False):
        """Initialize a BuildingBlockSet"""
        self.building_blocks = building_blocks
        self.num_building_blocks = len([bb for bb in self.building_blocks if bb.is_real()])

        self.has_smiles = all([bb.has_smiles() for bb in self.building_blocks if (not bb.is_null()) and bb.is_real()])

        if has_unique_smiles:
            if not self.has_smiles:
                raise BuildingBlockSetError(
                    "building block set is missing SMILES for some building blocks"
                )
            self._check_for_unique_chemicals()

        self._bb_subset_lookup_table: dict[str, list[T]] | None = self._build_subset_lookup_table()

    @classmethod
    def load(
        cls,
        uri_or_id: str | None = None,
        uri: str | None = None,
        id: str | None = None,
        load_smiles: bool = True,
        validate_smiles: bool = True,
        has_unique_smiles: bool = False,
        **kwargs
    ) -> "BuildingBlockSet":
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
        has_unique_smiles: bool
            Whether to check that all building blocks in the set have unique SMILES.
            Will also fail if any (non-null) blocks are missing SMILES
        """
        if load_smiles is False and validate_smiles is True:
            warnings.warn(
                "cannot validate chemicals if SMILES loading is disabled; skipping chemical validation", stacklevel=2
            )
            validate_smiles = False
        if kwargs:
            warnings.warn(
                f"unexpected keyword arguments {list(kwargs.keys())}; these will be ignored",
                stacklevel=1,
            )
        return super().load(
            uri_or_id=uri_or_id,
            uri=uri,
            id=id,
            validate_chemicals=validate_smiles,
            load_smiles=load_smiles,
            has_unique_smiles=has_unique_smiles,
        )

    @classmethod
    def _loads(cls, data: str, **kwargs) -> "BuildingBlockSet":
        """Parse the decoded text contents of a building block set file"""
        return cls(
            building_blocks=_parse_building_block_file_contents(
                contents=data,
                include_tags=False,
                delimiter=",",
                **kwargs
            ),
        )

    def _check_for_unique_chemicals(self):
        """Check that all building blocks in the set have unique chemical structures"""
        _inchi_key_set: set[str] = set()
        for bb in self.building_blocks:
            if bb.inchi_key in _inchi_key_set:
                raise BuildingBlockSetError(
                    f"building block set contains multiple building "
                    f"blocks with the same chemical structure (inchi key '{bb.inchi_key}')"
                )
                break  # only need to warn once about duplicates
            _inchi_key_set.add(bb.inchi_key)

    def _build_subset_lookup_table(self) -> dict[str, list[T]] | None:
        """Build a lookup table from subset id to list of building blocks in that subset"""
        _no_subset_id: bool = self.building_blocks[0].subset_id is None
        _bb_subset_lookup = defaultdict(list)
        for bb in self.building_blocks:
            if (bb.subset_id is None and not _no_subset_id) or (bb.subset_id is not None and _no_subset_id):
                raise BuildingBlockSetError(
                    "Building block set has inconsistent subset id usage"
                )
            elif bb.subset_id is not None:
                _bb_subset_lookup[bb.subset_id].append(bb)

        if _no_subset_id:
            return None
        subset_lookup: dict[str, list[T]] = dict(_bb_subset_lookup)  # cast to dict from default dict
        return subset_lookup

    def __len__(self):
        """Get the number of *real* building blocks in the set"""
        return self.num_building_blocks

    def __iter__(self):
        """Iterate over the building blocks"""
        return iter(self.building_blocks)

    def get_bb_by_id(self, query: str) -> BuildingBlock | None:
        """
        Alias for `get_object_by_id` to improve readability when working with building blocks

        Notes
        -----
        If you need better control over the behavior when a building block
        is not found, you can use `get_object_by_id`

        Parameters
        ----------
        query: str
            bb_id to query

        Returns
        -------
        BuildingBlock | None
            will be `None` if no matching building block is found
            else the matching BuildingBlock object
        """
        return self.get_obj_by_id(query)

    def _get_subset_lookup_table(self) -> dict[str, list[T]]:
        """Helper function to check if subset lookup table is initialized"""
        if self._bb_subset_lookup_table is None:
            raise ValueError("Building block set was not initialized with building blocks with subset ids")
        return self._bb_subset_lookup_table

    def get_bb_subset(self, subset_id: str, drop_null: bool = False) -> list[T]:
        """
        Given a subset id, return all building blocks in that subset

        Notes
        -----
        Order of returned blocks will match order of blocks in the set

        Parameters
        ----------
        subset_id: str
            the subset id to query
        drop_null: bool, default = False
            if `True` will drop null building blocks from the returned list

        Returns
        -------
        list[BuildingBlock]
            list of building blocks in that subset

        Raises
        ------
        ValueError
            if the BuildingBlockSet lacks building blocks with subset ids
        KeyError
            if the subset_id is not found in the BuildingBlockSet
        """
        try:
            bb_list = self._get_subset_lookup_table()[subset_id]
        except KeyError as e:
            raise KeyError(f"Building block set has no building blocks with subset id '{subset_id}'") from e

        if drop_null:
            return [bb for bb in bb_list if not bb.is_null()]
        else:
            return bb_list

    def get_bb_subsets(self, drop_null: bool = False) -> dict[str, list[T]]:
        """
        Get all building block subsets in the building block set

        Will map the subset id to the list of building blocks in that subset

        Returns
        -------
        dict[str, list[BuildingBlock]]
            mapping from building block subset id to list of building blocks in that subset
        drop_null: bool, default = False
            if `True` will drop null building blocks from the returned lists

        Raises
        ------
        ValueError
            if the BuildingBlockSet lacks building blocks with subset ids
        """
        # this is a little inefficient since it will check if the subset table exists each time
        # but this function isn't called in any performance critical functions so it not a big deal
        return {
            subset_id: self.get_bb_subset(subset_id, drop_null=drop_null)
            for subset_id in self._get_subset_lookup_table()
        }

    def get_subset_with_bb_id(self, bb_id: str) -> str:
        """
        Return the subset id that contains the given building block id for this set

        Parameters
        ----------
        bb_id: str
            building block id to check

        Returns
        -------
        str
            the subset id that contains the building block

        Raises
        ------
        KeyError
            if the building block id is not found in any subset
        ValueError
            if the BuildingBlockSet lacks building blocks with subset ids
        """
        bb = self.get_obj_by_id(bb_id, fail_on_missing=True)
        if bb.subset_id is not None:
            return bb.subset_id
        raise ValueError(f"building block '{bb_id}' does not belong to a subset")

    def get_null_building_blocks(self) -> list[T]:
        """
        Get all null building blocks in the set

        Returns
        -------
        list[BuildingBlock]
            list of null building blocks in the set
        """
        return [bb for bb in self.building_blocks if bb.is_null()]

    def get_non_null_building_blocks(self) -> list[T]:
        """
        Get all the non-null building blocks in the set

        Returns
        -------
        list[BuildingBlock]
            list of non-null building blocks in the set
        """
        return [bb for bb in self.building_blocks if not bb.is_null()]

class TaggedBuildingBlockSet(
    BuildingBlockSet[TaggedBuildingBlock],
    TaggedDecodableObjectSet[TaggedBuildingBlock],
    sub_dir="building_blocks",
):
    """
    Define a set of building blocks

    Attributes
    ----------
    has_smiles: bool
        True if all BuildingBlocks in this set have a non-None `smiles` attribute

    Parameters
    ----------
    building_blocks: Sequence[TaggedBuildingBlock]
        list of building blocks in set
    has_unique_smiles: bool, default = False
        whether all building blocks in the set have unique SMILES.
        Will also fail if any (non-null) blocks are missing SMILES
    """

    def __init__(
        self,
        building_blocks: Sequence[TaggedBuildingBlock],
        has_unique_smiles: bool = False,
    ):
        """Initialize the tagged building block set"""
        super().__init__(building_blocks=building_blocks, has_unique_smiles=has_unique_smiles)

        self.tag_length = len(self.building_blocks[0].tags[0])
        for i, _bb in enumerate(self.building_blocks[1:]):
            for _tag in _bb.tags:
                if self.tag_length != len(_tag):
                    raise BuildingBlockSetError(
                        f"expected all tags to have length {self.tag_length}, "
                        f"but tag '{_bb.bb_id}' at row {i} has length {len(_tag)}"
                    )

    @classmethod
    def _loads(cls, data: str, *args, **kwargs) -> "BuildingBlockSet":
        """Parse the decoded text contents of a building block set file"""
        return cls(
            building_blocks=_parse_building_block_file_contents(
                contents=data,
                include_tags=True,
                delimiter=",",
                load_smiles=kwargs.get("load_smiles", True),
                validate_chemicals=kwargs.get("validate_smiles", True),
            ),
            has_unique_smiles=kwargs.get("has_unique_smiles", False),
        )

    def get_bb_by_tag(self, query: str) -> TaggedBuildingBlock | None:
        """
        Alias for `get_obj_by_tag`; returns the building block with the given tag if it exists, else `None`

        Parameters
        ----------
        query: str
            DNA bases to query

        Returns
        -------
        TaggedBuildingBlock | None
            will be `None` if no matching building block is found
            else the matching BuildingBlock object
        """
        return self.get_obj_by_tag(query=query)
