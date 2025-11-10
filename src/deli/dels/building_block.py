"""define building block classes"""

import abc
import os
from collections import defaultdict
from typing import Literal, Optional, Sequence, overload

from deli.configure import DeliDataLoadable, accept_deli_data_name, get_deli_config
from deli.utils.mol_utils import SmilesMixin


# possible column headers for building block csv files
BB_FILE_ID_COLUMN = "id"
BB_FILE_SMILES_COLUMN = "smiles"
BB_FILE_TAG_COLUMN = "tag"
BB_FILE_SUBSET_ID_COLUMN = "subset_id"


def load_bb_set_from_csv_file(path: str) -> "BuildingBlockSet | TaggedBuildingBlockSet":
    """
    Load a building block set from a file

    Will determine if the building block set is tagged or not based on presence of 'tag' column
    in the building block file

    Parameters
    ----------
    path: str
        path to building block set file

    Returns
    -------
    BuildingBlockSet | TaggedBuildingBlockSet
        the loaded building block set
        will be TaggedBuildingBlockSet if 'tag' column is present
    """
    with open(path, "r") as f:
        header = f.readline().strip().split(",")

        if BB_FILE_TAG_COLUMN in header:
            return TaggedBuildingBlockSet.load_from_csv(path)
        else:
            return BuildingBlockSet.load_from_csv(path)


def _validate_bb_set_id(bb_set_id: str) -> bool:
    """Check that a building block set id is valid (not in subset id format)"""
    try:
        parse_building_block_subset_id(bb_set_id)
        return False
    except ValueError:
        return True  # true means the id is good


def _validate_bb_set_file_header(
    header: list[str], required_columns: set[str]
) -> tuple[int, int | None, int | None, dict[str, int]]:
    """
    Locates all required/optional columns in a building block set file header

    Returns the indices of the required columns and any extra required columns

    Returns
    -------
    tuple[int, int | None, int | None, dict[str, int]]
        tuple of (bb_id_col_idx, bb_smiles_col_idx, bb_subset_id_col, extra_cols)
        where extra_cols is a dict mapping from extra required column name to its index
    """
    if BB_FILE_ID_COLUMN not in header:
        raise BuildingBlockSetError(f"missing required column '{BB_FILE_ID_COLUMN}'")
    else:
        bb_col_idx = header.index(BB_FILE_ID_COLUMN)

    bb_smiles_col_idx: int | None
    if BB_FILE_SMILES_COLUMN in header:
        bb_smiles_col_idx = header.index(BB_FILE_SMILES_COLUMN)
    else:
        if BB_FILE_SMILES_COLUMN in required_columns:
            raise BuildingBlockSetError(f"missing required column '{BB_FILE_SMILES_COLUMN}'")
        bb_smiles_col_idx = None

    bb_subset_id_col_idx: int | None
    if BB_FILE_SUBSET_ID_COLUMN in header:
        bb_subset_id_col_idx = header.index(BB_FILE_SUBSET_ID_COLUMN)
    else:
        if BB_FILE_SUBSET_ID_COLUMN in required_columns:
            raise BuildingBlockSetError(f"missing required column '{BB_FILE_SUBSET_ID_COLUMN}'")
        bb_subset_id_col_idx = None

    # handle extra required columns
    required_columns = set(required_columns) - {
        BB_FILE_ID_COLUMN,
        BB_FILE_SMILES_COLUMN,
        BB_FILE_SUBSET_ID_COLUMN,
    }
    extra_cols: dict[str, int] = {}
    for required_column in required_columns:
        if required_column not in header:
            raise BuildingBlockSetError(f"missing required column '{required_column}'")
        else:
            extra_cols[required_column] = header.index(required_column)

    return bb_col_idx, bb_smiles_col_idx, bb_subset_id_col_idx, extra_cols


def generate_building_block_subset_id(bb_set_id: str, subset_id: str) -> str:
    """
    Generate a building block subset id given the building block set id and subset id

    Building block subset ids are formatted as '{bb_set_id}:::{subset_id}'

    Notes
    -----
    This function should always be used to generate building block subset ids to
    ensure consistent formatting.

    Parameters
    ----------
    bb_set_id: str
        building block set id
    subset_id: str
        building block subset id

    Returns
    -------
    str
        building block subset id in the format '{bb_set_id}:::{subset_id}'
    """
    return f"{bb_set_id}:::{subset_id}"


def parse_building_block_subset_id(bb_subset_id: str) -> tuple[str, str]:
    """
    Parse a building block subset id into its building block set id and subset id components

    Building block subset ids are formatted as '{bb_set_id}:::{bb_subset_id}'

    Parameters
    ----------
    bb_subset_id: str
        building block subset id in the format '{bb_set_id}:::{bb_subset_id}'

    Returns
    -------
    tuple[str, str]
        tuple of (bb_set_id, bb_subset_id)

    Raises
    ------
    ValueError
        if the bb_subset_id is not in the correct format
    """
    splits = bb_subset_id.split(":::")
    if len(splits) != 2:
        raise ValueError(
            f"building block subset id '{bb_subset_id}' is not in the correct format "
            f"'<bb_set_id>:::<subset_id>'"
        )
    return splits[0], splits[1]


class BuildingBlockSetError(Exception):
    """raised when there is an issue with building block set definition"""

    pass


class BaseBuildingBlock(abc.ABC):
    """Base class for all BuildingBlock objects"""

    def __init__(self):
        self.bb_id = None

    @abc.abstractmethod
    def is_mask(self) -> bool:
        """
        Check if building block is a masked building block or not

        Returns
        -------
        bool
        """
        raise NotImplementedError

    def is_real(self) -> bool:
        """
        Check if building block is a real (not null or mask) building block or not

        Returns
        -------
        bool
        """
        return not self.is_mask()

    @abc.abstractmethod
    def __eq__(self, other):
        """Determine if two building blocks are equal"""
        raise NotImplementedError

    def __str__(self):
        """Cast to str as the bb_id"""
        return self.bb_id


class MaskedBuildingBlock(BaseBuildingBlock):
    """
    Defines a Masked Building Block

    Notes
    -----
    A masked building block is a placeholder used in synthon-compounds
    where a specific building block that was present and real is masked
    out for this analysis

    Attributes
    ----------
    bb_id : str = BB_MASK
        always set to the constant BB_MASK
    """

    def __init__(self):
        """Initialize the object"""
        super().__init__()
        self.bb_id = get_deli_config()["BB_MASK"]

    def is_mask(self) -> bool:
        """Masked BBs are always masks"""
        return True

    # def is_null(self) -> bool:
    #     """Masked BBs are never null"""
    #     return False

    def __eq__(self, other):
        """Masked BBs can never be equal to anything"""
        return False

    def __copy__(self):
        """Return a new masked BB object"""
        return self.__class__()


class BuildingBlock(BaseBuildingBlock, SmilesMixin):
    """
    Define building block that has a chemical structure associated with it
    """

    def __init__(self, bb_id: str, smiles: Optional[str] = None, subset_id: Optional[str] = None):
        """
        Initialize the building block that has an associated DNA tag

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
        super().__init__()
        self.bb_id = bb_id
        self._smiles = smiles
        self.subset_id = str(subset_id) if subset_id else None  # stringify

    def __eq__(self, other):
        """Two building blocks are equal if their ID is equal"""
        return self.bb_id == other.bb_id if isinstance(other, BuildingBlock) else False

    @classmethod
    def from_dict(cls, data: dict[str, str]):
        """
        Initialize the building block from a dictionary

        Parameters
        ----------
        data: dict[str, str]
            Dictionary of building block properties

        Returns
        -------
        BuildingBlock
        """
        return cls(
            data["id"], smiles=data.get("smiles", None), subset_id=data.get("subset_id", None)
        )

    def is_mask(self) -> bool:
        """Masked BBs are never masks"""
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


class TaggedBuildingBlock(BuildingBlock):
    """
    Define building block that has a chemical structure and a DNA tag associated with it

    Attributes
    ----------
    tags: list[str]
        the dna tag(s) associated with this building block
    """

    def __init__(
        self,
        bb_id: str,
        tag: str | list[str],
        smiles: Optional[str] = None,
        subset_id: Optional[str] = None,
    ):
        """
        Initialize the building block

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
        super().__init__(bb_id=bb_id, smiles=smiles, subset_id=subset_id)
        self.tags = tag if not isinstance(tag, str) else [tag]

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
        return cls(bb_id=data["id"], tag=data["tag"], smiles=data.get("smiles", None))


class BuildingBlockSet(DeliDataLoadable):
    """
    Holds a set of building blocks

    Attributes
    ----------
    has_smiles: bool
        True if all BuildingBlocks in this set have SMILES
        does not guarantee SMILES are valid

    """

    def __init__(self, bb_set_id: str, building_blocks: Sequence[BuildingBlock]):
        """
        Initialize the building blocks

        Parameters
        ----------
        bb_set_id: str
            building block set id/name
        building_blocks: Sequence[BuildingBlock]
            list of building blocks in set
        """
        self.bb_set_id = bb_set_id

        # validate bb_set_id
        if not _validate_bb_set_id(bb_set_id):
            raise BuildingBlockSetError(
                f"building block set id '{bb_set_id}' is not valid, "
                f"contains characters reserved for building block subset id format. "
                f"See the docs for more details."
            )

        self.building_blocks = building_blocks
        self._validate_bb_set()

        self._bb_lookup_table = {bb.bb_id: i for i, bb in enumerate(self.building_blocks)}
        self._bb_subset_lookup_table: dict[str, list[BuildingBlock]] | None = (
            self._build_subset_lookup_table()
        )

        self.has_smiles = all([bb.has_smiles() for bb in self.building_blocks])

    def _validate_bb_set(self):
        """Validate the building block set"""
        _bb_id_set: set[str] = set()
        for bb in self.building_blocks:
            if bb.bb_id in _bb_id_set:
                raise BuildingBlockSetError(
                    f"duplicate building block id '{bb.bb_id}'; "
                    f"BuildingBlockSets cannot contain multiple building blocks "
                    f"with the same id"
                )
            _bb_id_set.add(bb.bb_id)

    def _build_subset_lookup_table(self) -> dict[str, list[BuildingBlock]] | None:
        """Build a lookup table from subset id to list of building blocks in that subset"""
        _no_subset_id: bool = self.building_blocks[0].subset_id is None
        _bb_subset_lookup = defaultdict(list)
        for bb in self.building_blocks:
            if bb.subset_id is None:
                if _no_subset_id:
                    continue
                else:
                    raise BuildingBlockSetError(
                        f"building block {bb.bb_id} in set '{self.bb_set_id}' is missing "
                        f"a subset id when other building blocks have them. All building "
                        f"blocks in a set must either have subset ids or none of them "
                        f"can have subset ids."
                    )
            elif _no_subset_id:
                raise BuildingBlockSetError(
                    f"building block {bb.bb_id} in set '{self.bb_set_id}' has a "
                    f"a subset id when other building blocks lack them. All building "
                    f"blocks in a set must either have subset ids or none of them "
                    f"can have subset ids."
                )
            else:
                _bb_subset_lookup[bb.subset_id].append(bb)

        if _no_subset_id:
            return None
        else:
            subset_lookup = dict(_bb_subset_lookup)  # cast to dict from default dict

            # assert building block set cannot be the same as subset ids
            if self.bb_set_id in _bb_subset_lookup.keys():
                raise BuildingBlockSetError(
                    f"building block set id '{self.bb_set_id}' cannot be the same as a subset id"
                )

            return subset_lookup

    @classmethod
    @accept_deli_data_name(sub_dir="building_blocks", extension="csv")
    def load(cls, path: str, check_for_smiles: bool = False) -> "BuildingBlockSet":
        """
        Load a building block set from the DELi data directory

        Notes
        -----
        This is decorated by `accept_deli_data`
        which makes allows to load by build block set name or path

        Parameters
        ----------
        path: str
            path of the building block set to load
        check_for_smiles: bool, default = False
            if `True` will check that the building block file has a SMILES column
            *will not* check that all SMILES are valid or present for all compounds

        Returns
        -------
        BuildingBlockSet
        """
        _cls = cls.load_from_csv(
            path, set_id=os.path.basename(path).split(".")[0], check_for_smiles=check_for_smiles
        )
        _cls.loaded_from = path
        return _cls

    @classmethod
    def load_from_csv(
        cls, path: str, set_id: Optional[str] = None, check_for_smiles: bool = False
    ) -> "BuildingBlockSet":
        """
        Read a building block set from a csv file

        Parameters
        ----------
        path: str
            path to csv file
        set_id: str, default = None
            An ID for the building block set
            if set_id will be the basename of the file if not passed
        check_for_smiles: bool, default = False
            if `True` will check that the building block file has a SMILES column
            *will not* check that all SMILES are valid or present for all compounds

        Returns
        -------
        BuildingBlockSet
        """
        # get set id name from file name if None
        if set_id is None:
            _set_id = os.path.basename(path).split(".")[0]
        else:
            _set_id = set_id

        _building_blocks: list[BuildingBlock] = list()
        _building_block_map: dict[str, BuildingBlock] = dict()
        with open(path, "r") as f:
            header = f.readline().strip().split(",")

            try:
                _id_col_idx, _smi_col_idx, _subset_id_col_idx, _ = _validate_bb_set_file_header(
                    header=header,
                    required_columns={BB_FILE_SMILES_COLUMN} if check_for_smiles else set(),
                )
            except BuildingBlockSetError as e:
                raise BuildingBlockSetError(
                    f"missing required columns in building block set '{_set_id}': {str(e)}"
                ) from e

            for line in f:
                splits = line.strip().split(",")
                _id = splits[_id_col_idx]
                _smiles = splits[_smi_col_idx] if _smi_col_idx is not None else None
                _subset = splits[_subset_id_col_idx] if _subset_id_col_idx is not None else None

                if _id in _building_block_map.keys():
                    if _smiles != _building_block_map[_id].smi:
                        raise BuildingBlockSetError(
                            f"duplicate building block id '{_id}' with conflicting SMILES: "
                            f"'{_smiles}', '{_building_block_map[_id].smi}' "
                        )
                    elif _subset == _building_block_map[_id].subset_id:
                        raise BuildingBlockSetError(
                            f"duplicate building block id '{_id}' with conflicting subset_ids: "
                            f"'{_subset}', '{_building_block_map[_id].subset_id}' "
                        )  # skip duplicate with same smiles
                    else:
                        continue  # skip since it's a duplicate
                else:
                    _building_blocks.append(
                        BuildingBlock(bb_id=_id, smiles=_smiles, subset_id=_subset)
                    )
        return cls(_set_id, _building_blocks)

    def __len__(self):
        """Get the number of building blocks in the set"""
        return len(self.building_blocks)

    def __iter__(self):
        """Iterate over the building blocks"""
        return iter(self.building_blocks)

    @overload
    def get_bb_by_id(
        self, query: str, fail_on_missing: Literal[False]
    ) -> BuildingBlock | None: ...

    @overload
    def get_bb_by_id(self, query: str, fail_on_missing: Literal[True]) -> BuildingBlock: ...

    def get_bb_by_id(self, query: str, fail_on_missing: bool = False) -> BuildingBlock | None:
        """
        Given a bb_id, search for corresponding BB for that ID

        Notes
        -----
        Will return `None` if no matching building block is found

        Parameters
        ----------
        query: str
            bb_id to query
        fail_on_missing: bool, default False
            if `True` raise a KeyError is no match is found
            else return `None`

        Returns
        -------
        Optional[BuildingBlock]
            will be `None` if no matching building block is found
            else the matching BuildingBlock object

        """
        _idx = self._bb_lookup_table.get(query, None)
        if _idx is None:
            if fail_on_missing:
                raise KeyError(
                    f"BuildingBlock id '{query}' not found in BuildingBlockSet '{self.bb_set_id}'"
                )
            return None
        else:
            return self.building_blocks[_idx]

    def _get_subset_lookup_table(self) -> dict[str, list[BuildingBlock]]:
        """Get the subset lookup table if it exists, else raise an error"""
        if self._bb_subset_lookup_table is None:
            raise ValueError(
                f"BuildingBlockSet '{self.bb_set_id}' was not initialized with a subset_id_map"
            )
        else:
            return self._bb_subset_lookup_table

    def get_bb_subset(self, subset_id: str) -> list[BuildingBlock]:
        """
        Given a subset id, return all building blocks in that subset

        Notes
        -----
        Will return all building blocks in the set if the subset_id matches the bb_set_id
        Can take subset ids in the full building block subset id format or as just the
        subset id

        Parameters
        ----------
        subset_id: str
            the subset id to query

        Returns
        -------
        list[BuildingBlock]
            list of building blocks in that subset

        Raises
        ------
        ValueError
            if the BuildingBlockSet was not initialized with a subset_id_map
        KeyError
            if the subset_id is not found in the BuildingBlockSet
        """
        try:
            bb_set_id, subset_id = parse_building_block_subset_id(subset_id)
            if bb_set_id != self.bb_set_id:
                raise KeyError(
                    f"Building block subset id '{subset_id}' does not belong "
                    f"to BuildingBlockSet '{self.bb_set_id}'"
                )
        except ValueError:
            # not in bb_subset_id format, assume just subset_id
            pass

        if self._bb_subset_lookup_table is None:
            if subset_id == self.bb_set_id:
                return list(self.building_blocks)
        try:
            bb_list = self._get_subset_lookup_table()[subset_id]
        except KeyError as e:
            raise KeyError(
                f"Subset id '{subset_id}' not found in BuildingBlockSet '{self.bb_set_id}'"
            ) from e
        return bb_list

    def get_bb_subsets(self) -> dict[str, list[BuildingBlock]]:
        """
        Get all building block subsets in the set

        Will map the full building block subset id (which includes the bb_set_id,
        not just the subset id) to the list of building blocks in that subset

        Returns
        -------
        dict[str, list[BuildingBlock]]
            mapping from building block subset id to list of building blocks in that subset

        See Also
        --------
        generate_building_block_subset_id
        """
        return {
            generate_building_block_subset_id(self.bb_set_id, subset_id): subset
            for subset_id, subset in self._get_subset_lookup_table().items()
        }

    def get_subset_with_bb(self, bb: BuildingBlock, as_bb_subset_id: bool = False) -> str:
        """
        Return the subset id that contains the given building block

        Parameters
        ----------
        bb: BuildingBlock
            building block to check
        as_bb_subset_id: bool, default = False
            if `True` will return the building block subset id format

        Returns
        -------
        str
            the subset id that contains the building block

        Raises
        ------
        KeyError
            if the building block is not found in any subset
        """
        return self.get_subset_with_bb_id(bb.bb_id, as_bb_subset_id=as_bb_subset_id)

    def get_subset_with_bb_id(self, bb_id: str, as_bb_subset_id: bool = False) -> str:
        """
        Return the subset id that contains the given building block id

        Notes
        -----
        If building block set has no subsets, will return the building block set id
        if the building block is found in the set, otherwise it will raise an error

        Parameters
        ----------
        bb_id: str
            building block id to check
        as_bb_subset_id: bool, default = False
            if `True` will return the building block subset id format

        Returns
        -------
        str
            the subset id that contains the building block

        Raises
        ------
        KeyError
            if the building block id is not found in any subset
        """
        if self._bb_subset_lookup_table is not None:
            for subset_id, bbs in self._bb_subset_lookup_table.items():
                for bb in bbs:
                    if bb.bb_id == bb_id:
                        if as_bb_subset_id:
                            return generate_building_block_subset_id(self.bb_set_id, subset_id)
                        else:
                            return subset_id
        else:
            if any(bb_id == _bb.bb_id for _bb in self.building_blocks):
                return self.bb_set_id
        raise KeyError(
            f"Building block '{bb_id}' not found in any subset "
            f"of BuildingBlockSet '{self.bb_set_id}'"
        )


class TaggedBuildingBlockSet(BuildingBlockSet):
    """
    Define a set of building blocks

    Attributes
    ----------
    has_smiles: bool
        True if all BuildingBlocks in this set have a non-None `smiles` attribute
    """

    def __init__(
        self,
        bb_set_id: str,
        building_blocks: Sequence[TaggedBuildingBlock],
    ):
        """
        Initialize the building blocks

        Parameters
        ----------
        bb_set_id: str
            building block set id/name
        building_blocks: Sequence[TaggedBuildingBlock]
            list of building blocks in set

        Notes
        -----
        Building blocks that are attached to DNA should mark the atom
        that bind DNA with a "X" in the SMILES
        """
        super().__init__(bb_set_id, building_blocks)
        self.building_blocks: Sequence[TaggedBuildingBlock] = building_blocks  # for type checker

        self.tag_length = len(self.building_blocks[0].tags[0])
        for i, _bb in enumerate(self.building_blocks[1:]):
            for _tag in _bb.tags:
                if self.tag_length != len(_tag):
                    raise BuildingBlockSetError(
                        f"expected all tags to have length {self.tag_length}, "
                        f"but tag '{_bb.bb_id}' at row {i} has length {len(_tag)}"
                    )

        self._dna_lookup_table = {
            tag: i for i, bb in enumerate(self.building_blocks) for tag in bb.tags
        }
        self._bb_lookup_table = {bb.bb_id: i for i, bb in enumerate(self.building_blocks)}

    def _validate_bb_set(self):
        """Validate the building block set"""
        _bb_id_set: set[str] = set()
        _tag_set: set[str] = set()
        for bb in self.building_blocks:
            if bb.bb_id in _bb_id_set:
                raise BuildingBlockSetError(
                    f"duplicate building block id '{bb.bb_id}'; "
                    f"BuildingBlockSets cannot contain multiple building blocks "
                    f"with the same id"
                )
            _bb_id_set.add(bb.bb_id)

            for tag in bb.tags:
                if tag in _tag_set:
                    raise BuildingBlockSetError(
                        f"duplicate building block tag '{tag}'; "
                        f"BuildingBlockSets cannot contain multiple building blocks "
                        f"with the same tag"
                    )
                _tag_set.add(tag)

    @classmethod
    @accept_deli_data_name(sub_dir="building_blocks", extension="csv")
    def load(cls, path: str, check_for_smiles: bool = False) -> "TaggedBuildingBlockSet":
        """
        Load a tagged building block set from the DELi data directory

        Notes
        -----
        This is decorated by `accept_deli_data`
        which makes allows to load by build block set name or path

        Parameters
        ----------
        path: str
            path of the building block set to load
        check_for_smiles: bool, default = False
            if `True` will check that the building block file has a SMILES column
            *will not* check that all SMILES are valid or present for all compounds

        Returns
        -------
        TaggedBuildingBlockSet
        """
        _cls = cls.load_from_csv(path, set_id=os.path.basename(path).split(".")[0])
        _cls.loaded_from = path
        return _cls

    @classmethod
    def load_from_csv(
        cls, path: str, set_id: Optional[str] = None, check_for_smiles: bool = False
    ) -> "TaggedBuildingBlockSet":
        """
        Read a building block set from a csv file

        Parameters
        ----------
        path: str
            path to csv file
        set_id: str, default = None
            An ID for the building block set
            if set_id will be the basename of the file if not passed
        check_for_smiles: bool, default = False
            if `True` will check that the building block file has a SMILES column
            *will not* check that all SMILES are valid or present for all compounds

        Returns
        -------
        TaggedBuildingBlockSet
        """
        # get set id name from file name if None
        if set_id is None:
            _set_id = os.path.basename(path).split(".")[0]
        else:
            _set_id = set_id

        _building_blocks: list[TaggedBuildingBlock] = list()
        _building_block_map: dict[str, TaggedBuildingBlock] = dict()
        with open(path, "r") as f:
            header = f.readline().strip().split(",")

            try:
                _id_col_idx, _smi_col_idx, _subset_id_col_idx, extra_cols = (
                    _validate_bb_set_file_header(
                        header=header,
                        required_columns={BB_FILE_SMILES_COLUMN, BB_FILE_TAG_COLUMN}
                        if check_for_smiles
                        else {BB_FILE_TAG_COLUMN},
                    )
                )
            except BuildingBlockSetError as e:
                raise BuildingBlockSetError(
                    f"missing required columns in building block set '{_set_id}': {str(e)}"
                ) from e

            for line in f:
                splits = line.strip().split(",")
                _id = splits[_id_col_idx]
                _smiles = splits[_smi_col_idx] if _smi_col_idx is not None else None
                _subset = splits[_subset_id_col_idx] if _subset_id_col_idx is not None else None
                _tag = splits[extra_cols[BB_FILE_TAG_COLUMN]]

                if _id in _building_block_map.keys():
                    if _smiles != _building_block_map[_id].smi:
                        raise BuildingBlockSetError(
                            f"duplicate building block id '{_id}' with conflicting SMILES: "
                            f"'{_smiles}', '{_building_block_map[_id].smi}' "
                        )
                    elif _subset != _building_block_map[_id].subset_id:
                        raise BuildingBlockSetError(
                            f"duplicate building block id '{_id}' with conflicting subset_ids: "
                            f"'{_subset}', '{_building_block_map[_id].subset_id}' "
                        )  # skip duplicate with same smiles
                    else:
                        _building_block_map[_id].tags.append(_tag)  # skip since it's a duplicate
                else:
                    _building_blocks.append(
                        TaggedBuildingBlock(bb_id=_id, smiles=_smiles, tag=_tag, subset_id=_subset)
                    )
        return cls(_set_id, _building_blocks)

    @overload
    def search_tags(
        self, query: str, fail_on_missing: Literal[False]
    ) -> TaggedBuildingBlock | None: ...

    @overload
    def search_tags(self, query: str, fail_on_missing: Literal[True]) -> TaggedBuildingBlock: ...

    def search_tags(self, query: str, fail_on_missing: bool = False) -> TaggedBuildingBlock | None:
        """
        Given a query DNA bases, search for corresponding BB with that bases

        Notes
        -----
        Will return `None` if no matching building block is found

        Parameters
        ----------
        query: str
            DNA bases to query
        fail_on_missing: bool, default False
            if `True` raise a KeyError is no match is found
            else return `None`

        Returns
        -------
        TaggedBuildingBlock | None
            will be `None` if no matching building block is found
            else the matching BuildingBlock object

        """
        _idx = self._dna_lookup_table.get(query, None)
        if _idx is None:
            if fail_on_missing:
                raise KeyError(
                    f"BuildingBlock DNA tag '{query}' not found "
                    f"in BuildingBlockSet '{self.bb_set_id}'"
                )
            return None
        else:
            return self.building_blocks[_idx]
