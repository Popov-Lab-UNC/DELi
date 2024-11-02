"""define building block classes"""

import abc
import os
import warnings
from typing import List, Optional, Self

from deli.configure import DeliDataLoadable, accept_deli_data_name
from deli.constants import BB_MASK


class BaseBuildingBlock(abc.ABC):
    """Base class for all BuildingBlock objects"""

    def __init__(self):
        self.bb_id = None
        self.smiles = None

    @abc.abstractmethod
    def is_mask(self) -> bool:
        """
        check if building block is a masked building block or not

        Returns
        -------
        bool
        """
        raise NotImplementedError

    def is_real(self) -> bool:
        """
        check if building block is a real (not null or mask) building block or not

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
        self.bb_id = BB_MASK

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
    """define building block class"""

    def __init__(self, bb_id: str, smiles: Optional[str] = None, tag: Optional[str] = None):
        """
        Initialize the building block

        Parameters
        ----------
        bb_id: str
            building block id
        smiles: Optional[str], default = None
            SMILES of the building block
        """
        super().__init__()
        self.bb_id = bb_id
        self.smiles = smiles
        self.tag = tag

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
        return cls(data["id"], data.get("smiles"), data.get("bases"))

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

        self._dna_lookup_table = {bb.tag: i for i, bb in enumerate(self.building_blocks)}

    @classmethod
    @accept_deli_data_name(sub_dir="building_blocks", extension="csv")
    def load(cls, path: str) -> Self:
        """
        Load a building block set from the DELi data directory

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
            path of the building block set to load

        Returns
        -------
        BuildingBlockSet
        """
        return cls.load_from_csv(path)

    @classmethod
    def load_from_csv(cls, path: str, set_id: Optional[str] = None, no_smiles: bool = False):
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
            _id_col_idx = header.index("id")
            _dna_col_idx = header.index("tag")

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
                _building_blocks.append(BuildingBlock(bb_id=_id, smiles=_smiles, tag=_dna))
        return cls(_set_id, _building_blocks)

    def __len__(self):
        """Get the number of building blocks in the set"""
        return len(self.building_blocks)

    def __iter__(self):
        """Iterate over the building blocks"""
        return iter(self.building_blocks)

    def search_tags(self, query: str) -> Optional[BuildingBlock]:
        """
        Given a query DNA bases, search for corresponding BB with that bases

        Notes
        -----
        Will return `None` if no matching building block is found

        Parameters
        ----------
        query: str
            DNA bases to query

        Returns
        -------
        Optional[BuildingBlock]
            will be `None` if no matching building block is found
            else the matching BuildingBlock object

        """
        _idx = self._dna_lookup_table.get(query, None)
        if _idx is None:
            return None
        else:
            return self.building_blocks[_idx]
