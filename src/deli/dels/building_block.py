"""define building block classes"""

import abc
from os import PathLike
from typing import List, Optional, Union

from deli.constants import BB_MASK
from deli.dels.configure import validate_file_path


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

    # @abc.abstractmethod
    # def is_null(self) -> bool:
    #     """
    #     check if building block is a null building block or not
    #
    #     Returns
    #     -------
    #     bool
    #     """
    #     raise NotImplementedError

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


# class NullBuildingBlock(BaseBuildingBlock):
#     """
#     Defines a Null Building Block
#
#     Notes
#     -----
#     A masked building block is a place holder used to specify that
#     there was no building block here. This is used when trying to
#     force a 2 cycle library into a 3 cycle library
#
#     Attributes
#     ----------
#     bb_id : str = BB_NULL
#         always set to the constant BB_NULL
#     """
#
#     def __init__(self):
#         super().__init__()
#         self.bb_id = BB_NULL
#
#     def is_mask(self) -> bool:
#         """Masked BBs are never masks"""
#         return False
#
#     def is_null(self) -> bool:
#         """Masked BBs are always null"""
#         return True
#
#     def __eq__(self, other):
#         """Null BB can never be equal to another BB"""
#         return False
#
#     def __copy__(self):
#         """Return a new null BB object"""
#         return self.__class__()


class BuildingBlock(BaseBuildingBlock):
    """define building block class"""

    def __init__(self, bb_id: str, smiles: Optional[str] = None,
                 tag: Optional[str] = None):
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
        return cls(data["id"], data.get("smiles"), data.get("tag"))

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


class BuildingBlockSet:
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

    # @classmethod
    # def from_json(cls, path: str):
    #     """
    #     Initialize the building block set from a json file
    #
    #     Parameters
    #     ----------
    #     path: str
    #         path to json file
    #
    #     Returns
    #     -------
    #     BuildingBlockSet
    #     """
    #     data = json.load(open(path))
    #     return cls(
    #         data["id"], [BuildingBlock.from_dict(bb_data) for bb_data in data["building_blocks"]]
    #     )

    @classmethod
    @validate_file_path(sub_dir='building_blocks')
    def from_csv(cls, file_path: Union[str, PathLike], set_id: Optional[str] = None):
        _building_blocks = []
        with open(file_path, "r") as f:
            header = f.readline()
            _id_col_idx = header.index("id")
            _smi_col_idx = header.index("smiles")
            _dna_col_idx = header.index("tag")
            for line in f:
                splits = line.strip().split(",")
                _id = splits[_id_col_idx]
                _smiles = splits[_smi_col_idx]
                _dna = splits[_dna_col_idx]
                _building_blocks.append(BuildingBlock(
                    bb_id = _id,
                    smiles = _smiles,
                    tag = _dna
                ))
        return cls(set_id, _building_blocks)

    def __len__(self):
        """Get the number of building blocks in the set"""
        return len(self.building_blocks)

    def __iter__(self):
        """Iterate over the building blocks"""
        return iter(self.building_blocks)
