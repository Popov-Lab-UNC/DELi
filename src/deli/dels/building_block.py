"""define building block classes"""

import abc
import os
from typing import Literal, Optional, Sequence, overload

from deli.configure import DeliDataLoadable, accept_deli_data_name, get_deli_config
from deli.utils.mol_utils import SmilesMixin


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
    A masked building block is a place holder used in synthon-compounds
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
    """Define building block that has a chemical structure associated with it"""

    def __init__(self, bb_id: str, smiles: Optional[str] = None):
        """
        Initialize the building block that has an associated DNA tag

        Parameters
        ----------
        bb_id: str
            building block id
        smiles: Optional[str], default = None
            SMILES of the building block
        """
        super().__init__()
        self.bb_id = bb_id
        self._smiles = smiles

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
        return cls(data["id"], smiles=data.get("smiles", None))

    def is_mask(self) -> bool:
        """Masked BBs are never masks"""
        return False


class TaggedBuildingBlock(BuildingBlock):
    """Define building block that has a chemical structure and a DNA tag associated with it"""

    def __init__(self, bb_id: str, tag: str, smiles: Optional[str] = None):
        """
        Initialize the building block

        Parameters
        ----------
        bb_id: str
            building block id
        tag: str
            the dna tag associated with this building block
        smiles: Optional[str], default = None
            SMILES of the building block
        """
        super().__init__(bb_id=bb_id, smiles=smiles)
        self.tag = tag

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


class BuildingBlockSet(DeliDataLoadable, abc.ABC):
    """Holds a set of building blocks"""

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
        self.building_blocks = building_blocks
        self._bb_lookup_table = {bb.bb_id: i for i, bb in enumerate(self.building_blocks)}

        self.has_smiles = all([bb.has_smiles() for bb in self.building_blocks])

    @classmethod
    @accept_deli_data_name(sub_dir="building_blocks", extension="csv")
    def load(cls, path: str) -> "BuildingBlockSet":
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

        Returns
        -------
        BuildingBlockSet
        """
        _cls = cls.load_from_csv(path, set_id=os.path.basename(path).split(".")[0])
        _cls.loaded_from = path
        return _cls

    @classmethod
    def load_from_csv(cls, path: str, set_id: Optional[str] = None) -> "BuildingBlockSet":
        """
        Read a building block set from a csv file

        Parameters
        ----------
        path: str
            path to csv file
        set_id: str, default = None
            An ID for the building block set
            if set_id will be the basename of the file if not passed

        Returns
        -------
        BuildingBlockSet
        """
        # get set id name from file name if None
        if set_id is None:
            _set_id = os.path.basename(path).split(".")[0]
        else:
            _set_id = set_id

        _building_blocks = []
        with open(path, "r") as f:
            header = f.readline().strip().split(",")

            try:
                _id_col_idx = header.index("id")
            except ValueError as e:
                raise BuildingBlockSetError(
                    f"missing column 'id' of building block set '{_set_id}'"
                ) from e

            _smi_col_idx = header.index("smiles") if "smiles" in header else None

            for line in f:
                splits = line.strip().split(",")
                _id = splits[_id_col_idx]
                _smiles = splits[_smi_col_idx] if _smi_col_idx is not None else None
                _building_blocks.append(BuildingBlock(bb_id=_id, smiles=_smiles))

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


class TaggedBuildingBlockSet(BuildingBlockSet):
    """
    Define a set of building blocks

    Attributes
    ----------
    has_smiles: bool
        True if all BuildingBlocks in this set have a non-None `smiles` attribute
    """

    def __init__(self, bb_set_id: str, building_blocks: Sequence[TaggedBuildingBlock]):
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

        self.tag_length = len(self.building_blocks[0].tag)
        for i, _bb in enumerate(self.building_blocks[1:]):
            if self.tag_length != len(_bb.tag):
                raise BuildingBlockSetError(
                    f"expected all tags to have length {self.tag_length}, "
                    f"but tag '{_bb.bb_id}' at row {i} has length {len(_bb.tag)}"
                )

        self._dna_lookup_table = {bb.tag: i for i, bb in enumerate(self.building_blocks)}
        self._bb_lookup_table = {bb.bb_id: i for i, bb in enumerate(self.building_blocks)}

    @classmethod
    @accept_deli_data_name(sub_dir="building_blocks", extension="csv")
    def load(cls, path: str) -> "TaggedBuildingBlockSet":
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

        Returns
        -------
        TaggedBuildingBlockSet
        """
        _cls = cls.load_from_csv(path, set_id=os.path.basename(path).split(".")[0])
        _cls.loaded_from = path
        return _cls

    @classmethod
    def load_from_csv(cls, path: str, set_id: Optional[str] = None) -> "TaggedBuildingBlockSet":
        """
        Read a building block set from a csv file

        Parameters
        ----------
        path: str
            path to csv file
        set_id: str, default = None
            An ID for the building block set
            if set_id will be the basename of the file if not passed

        Returns
        -------
        TaggedBuildingBlockSet
        """
        # get set id name from file name if None
        if set_id is None:
            _set_id = os.path.basename(path).split(".")[0]
        else:
            _set_id = set_id

        _building_blocks = []
        with open(path, "r") as f:
            header = f.readline().strip().split(",")

            try:
                _id_col_idx = header.index("id")
            except ValueError as e:
                raise BuildingBlockSetError(
                    f"missing column 'id' of building block set '{_set_id}'"
                ) from e

            try:
                _dna_col_idx = header.index("tag")
            except ValueError as e:
                raise BuildingBlockSetError(
                    f"missing column 'tag' of tagged building block set '{_set_id}'"
                ) from e

            _smi_col_idx = header.index("smiles") if "smiles" in header else None

            for line in f:
                splits = line.strip().split(",")
                _id = splits[_id_col_idx]
                _smiles = splits[_smi_col_idx] if _smi_col_idx is not None else None
                _tag = splits[_dna_col_idx]
                _building_blocks.append(TaggedBuildingBlock(bb_id=_id, tag=_tag, smiles=_smiles))

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
