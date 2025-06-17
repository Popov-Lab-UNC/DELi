"""define building block classes"""

import abc
import os
import warnings
from typing import List, Literal, Optional, Self, overload

from deli.configure import DeliDataLoadable, accept_deli_data_name, get_deli_config
from deli.utils.mol_utils import to_mol


class BuildingBlockSetError(Exception):
    """raised when there is an issue with building block set definition"""

    pass


class BaseBuildingBlock(abc.ABC):
    """Base class for all BuildingBlock objects"""

    def __init__(self):
        self.bb_id = None
        self.smiles = None

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


class BuildingBlock(BaseBuildingBlock):
    """
    Define building block class

    Attributes
    ----------
    mol: Chem.Mol
        the rdkit mol for the building block
        if smiles is None will be None
    """

    def __init__(self, bb_id: str, tag: str, smiles: Optional[str] = None, load_mol: bool = False):
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
        load_mol: bool, default = False
            load the rdkit mol for this smiles
            requires smiles is not None
        """
        super().__init__()
        self.bb_id = bb_id
        self.smiles = smiles
        self.tag = tag

        self._mol = None
        if load_mol:
            if self.smiles is not None:
                self._mol = to_mol(self.smiles, fail_on_error=True)
            else:
                warnings.warn(
                    f"cannot load mol for building block {self.bb_id} because smiles is None",
                    stacklevel=1,
                )

    @property
    def mol(self):
        """The rdkit mol for the building block if smiles is None will be None"""
        if self._mol is None and self.smiles is not None:
            self._mol = to_mol(self.smiles, fail_on_error=True)
        return self._mol

    @mol.deleter
    def mol(self):
        """Set the mol property to None"""
        self._mol = None

    @classmethod
    def from_dict(cls, data: dict):
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
        return cls(data["id"], data["tag"], data.get("smiles"))

    def is_mask(self) -> bool:
        """Masked BBs are never masks"""
        return False

    def __eq__(self, other):
        """BBs are equal if they share the same id and SMILES"""
        if not isinstance(other, self.__class__):
            return False
        else:
            return (self.bb_id == other.bb_id) and (
                self.smiles == other.smiles if self.smiles else True
            )

    def __copy__(self):
        """Makes a new BB object with the same id and smiles"""
        return self.__class__(self.bb_id, self.smiles)


class BuildingBlockSet(DeliDataLoadable):
    """
    Define a set of building blocks

    Attributes
    ----------
    has_smiles: bool
        True if all BuildingBlocks in this set have a non-None `smiles` attribute
    """

    def __init__(self, bb_set_id: str, building_blocks: List[BuildingBlock]):
        """
        Initialize the building blocks

        Parameters
        ----------
        bb_set_id: str
            building block set id/name
        building_blocks: List[BuildingBlock]
            list of building blocks in set

        Notes
        -----
        Building blocks that are attached to DNA should mark the atom
        that bind DNA with a "X" in the SMILES
        """
        self.bb_set_id = bb_set_id
        self.building_blocks = building_blocks
        self.tag_length = len(self.building_blocks[0].tag)
        for i, _bb in enumerate(building_blocks[1:]):
            if self.tag_length != len(_bb.tag):
                raise BuildingBlockSetError(
                    f"expected all tags to have length {self.tag_length}, "
                    f"but tag '{_bb.bb_id}' at row {i} has length {len(_bb.tag)}"
                )

        self._dna_lookup_table = {bb.tag: i for i, bb in enumerate(self.building_blocks)}
        self._bb_lookup_table = {bb.bb_id: i for i, bb in enumerate(self.building_blocks)}

        self.has_smiles = False
        if all([bb.smiles for bb in self.building_blocks]):
            self.has_smiles = True

    @classmethod
    @accept_deli_data_name(sub_dir="building_blocks", extension="csv")
    def load(
        cls,
        path: str,
        set_id: Optional[str] = None,
        no_smiles: bool = False,
        load_mols: bool = False,
    ) -> Self:
        """
        Load a building block set from the DELi data directory

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
            path of the building block set to load
        set_id: Optional[str]
            name of the building block set
        no_smiles: bool, default False
            set to true if bb_set lacks a smiles column
        load_mols: bool, default False
            set to true if you want to load the mols for each building block
            this is useful when you expect to be operating on the BBs chemically
            for example, doing enumeration

        Returns
        -------
        BuildingBlockSet
        """
        _cls = cls.load_from_csv(path, set_id=set_id, no_smiles=no_smiles, load_mols=load_mols)
        _cls.loaded_from = path
        return _cls

    @classmethod
    def load_from_csv(
        cls,
        path: str,
        set_id: Optional[str] = None,
        no_smiles: bool = False,
        load_mols: bool = False,
    ) -> Self:
        """
        Read a building block set from a csv file

        Notes
        -----
        if `set_id` is not passed, will use the file name as the set id

        Parameters
        ----------
        path: str
            path to csv file
        set_id: Optional[str]
            name of the building block set
        no_smiles: bool, default False
            set to true if bb_set lacks a smiles column
        load_mols: bool, default False
            set to true if you want to load the mols for each building block
            this is useful when you expect to be operating on the BBs chemically
            for example, doing enumeration

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
            try:
                _dna_col_idx = header.index("tag")
            except ValueError as e:
                raise BuildingBlockSetError(
                    f"missing column 'tag' of building block set '{_set_id}'"
                ) from e

            if not no_smiles:
                if "smiles" not in header:
                    warnings.warn(
                        f"building block file {path} missing 'smiles' column in header; "
                        f"set `no_smiles` to `True` to turn off this message",
                        stacklevel=0,
                    )
                    no_smiles = True

            _smi_col_idx = header.index("smiles") if not no_smiles else None

            for line in f:
                splits = line.strip().split(",")
                _id = splits[_id_col_idx]
                _smiles = splits[_smi_col_idx] if _smi_col_idx is not None else None
                _dna = splits[_dna_col_idx]
                _building_blocks.append(
                    BuildingBlock(bb_id=_id, smiles=_smiles, tag=_dna, load_mol=load_mols)
                )
        return cls(_set_id, _building_blocks)

    def __len__(self):
        """Get the number of building blocks in the set"""
        return len(self.building_blocks)

    def __iter__(self):
        """Iterate over the building blocks"""
        return iter(self.building_blocks)

    @overload
    def search_tags(
        self, query: str, fail_on_missing: Literal[False]
    ) -> Optional[BuildingBlock]: ...

    @overload
    def search_tags(self, query: str, fail_on_missing: Literal[True]) -> BuildingBlock: ...

    def search_tags(self, query: str, fail_on_missing: bool = False) -> Optional[BuildingBlock]:
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
        Optional[BuildingBlock]
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

    @overload
    def get_bb_by_id(
        self, query: str, fail_on_missing: Literal[False]
    ) -> Optional[BuildingBlock]: ...

    @overload
    def get_bb_by_id(self, query: str, fail_on_missing: Literal[True]) -> BuildingBlock: ...

    def get_bb_by_id(self, query: str, fail_on_missing: bool = False) -> Optional[BuildingBlock]:
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
