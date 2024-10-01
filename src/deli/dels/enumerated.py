"""enumerated compound classes"""

from typing import List, Optional

from .building_block import BaseBuildingBlock
from .member import SelectableDELMember
from .synthon import HasDisynthonMixin


class EnumeratedCompoundError(Exception):
    """error raised when a compound cannot be made"""

    pass


class FullyEnumeratedCompound(SelectableDELMember, HasDisynthonMixin):
    """
    A full (no masked building block) compound from a DEL

    Notes
    -----
    Similar to the synthons, except for a complete compound with no masked
    building blocks
    """

    def __init__(self, building_blocks: List[BaseBuildingBlock], library_id: Optional[str] = None):
        """
        Initialize the object

        Parameters
        ----------
         building_blocks: List[BaseBuildingBlock]
            the building blocks that make up the synthon
            should be MaskedBuildingBlock for cycles were the BB info isn't needed
        library_id: str, optional
            the library ID for the DEL the synthon came from
            only necessary if BuildingBlock ids are not unique across libraries
        """
        super().__init__(building_blocks, library_id)

        if not all([bb.is_real() for bb in building_blocks]):
            raise EnumeratedCompoundError(
                "Fully Enumerated compounds can only have real building blocks"
            )

        self.compound_id = "-".join([building_block.bb_id for building_block in building_blocks])

        if self.library_id is not None:
            self.compound_id = self.library_id + "-" + self.compound_id

    def __eq__(self, other):
        """True if two FullyEnumeratedCompound objects have the same id"""
        if isinstance(other, FullyEnumeratedCompound):
            return self.compound_id == other.compound_id
        return False
