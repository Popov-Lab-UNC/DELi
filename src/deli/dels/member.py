"""generic abstract classes for chemicals/synthons in a DEL library"""

import abc
from typing import List, Optional

from deli.dels.building_block import BaseBuildingBlock


class SelectableMixIn:
    """mixin for adding selection data to object"""

    selection_count: Optional[float] = None
    control_count: Optional[float] = None
    competitive_count: Optional[float] = None
    selected = False

    def select(
        self,
        selection_count: Optional[float],
        control_count: Optional[float],
        competitive_count: Optional[float],
    ):
        """
        Add selection data to this object

        Notes
        -----
        once called at least once will set "selected" to True
        initializes all selection counts to -1 (impossible measurement)

        Parameters
        ----------
        selection_count: float, optional
            added a selection condition count
        control_count: float, optional
            added a control condition count
        competitive_count: float, optional
            added a competitive condition count
        """
        if selection_count is not None:
            self.selection_count = selection_count
        if control_count is not None:
            self.control_count = control_count
        if competitive_count is not None:
            self.competitive_count = competitive_count
        self.selected = True

    def been_selected(self) -> bool:
        """Return True if object has selection data attached"""
        return self.selected


class DELMember(abc.ABC):
    """
    Abstract class for all objects that derive from compounds in the DEL

    Notes
    -----
    This class includes fully enumerated compounds AND synthons as children

    Attributes
    ----------
    real_cycles: list[int]
        the cycle indexes with real building blocks
    """

    def __init__(
        self,
        building_blocks: List[BaseBuildingBlock],
        library_id: Optional[str] = None,
    ):
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
        self.building_blocks = building_blocks
        self.library_id = library_id
        self.real_cycles = [i for i, bb in enumerate(self.building_blocks) if bb.is_real()]

        if len(self.building_blocks) < 2:
            raise ValueError(
                "All DEL chemicals must have at least 2 building_blocks"
                " (including Masked building blocks)"
            )

    @abc.abstractmethod
    def __eq__(self, other):
        """Return true if other has have the same ID"""
        raise NotImplementedError


class SelectableDELMember(DELMember, SelectableMixIn, abc.ABC):
    """
    Abstract class for a DEL member that can have selected data attached to it
    """

    def __init__(
        self,
        building_blocks: List[BaseBuildingBlock],
        library_id: Optional[str] = None,
    ):
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

    def __add__(self, other):
        """Add the selection data of two identical del members together"""
        if not isinstance(other, self.__class__):
            raise TypeError(f"cannot add {self.__class__.__name__} to other.__class__.__name__")
        if self.building_blocks != other.building_blocks:
            raise ValueError("cannot add two DEL chemicals together if do not share the same ID")
        _new_obj = self.__class__(self.building_blocks, self.library_id)

        _new_obj.selection_count = self.selection_count + other.selection_count
        _new_obj.control_count = self.control_count + other.control_count
        _new_obj.competitive_count = self.competitive_count + other.competitive_count

        return _new_obj
