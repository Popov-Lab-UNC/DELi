"""synthon functionality"""

import abc
from copy import deepcopy
from itertools import combinations
from typing import List, Optional

from .building_block import BaseBuildingBlock, MaskedBuildingBlock
from .member import SelectableDELMember


class SynthonError(Exception):
    """Error for issues related to creating synthons"""

    pass


class Synthon(SelectableDELMember, abc.ABC):
    """
    abstract class for all synthon objects

    Notes
    -----
    Sythons are partial fully enumerated compounds, they only include 1 or 2
    Building blocks and the rest are masked out. This is different than just
    having a building block or a smaller fully enumerated DEL since synthons
    exist in the context of the original DEL library they come from

    For example: a Monosynthon vs a Building Block. Both only have information
    about a single building block, but the Monosynthon retains information about
    what cycle that building block was in the DEL. Further, Two monosythons
    with the same building block but at different cycles are NOT equal, but if
    they were building block objects we consider them equal

    We make this very minor distinction to help with analysis. Synthons are not
    "real"; they are not in the library. However, analyzing the result in a
    synthon based approach (for example looking for enrichment in disynthons
    over fully enumerated compounds) has shown to be an effective way to precess
    the selection data. We want to make sure these objects are treated differently
    because they are not real chemicals, unlike the fully enumerated compounds and
    building blocks, which are physically present

    Attributes
    ----------
    synthon_id: str
        the unique ID for the synthon
        based on passed building blocks
    real_cycles: str
        the cycles in the synthon with real building blocks
    num_cycles: int
        the number of cycles in the DEL the synthon came from
    real_cycles: list[int]
        the cycle indexes with real building blocks
    """

    def __init__(
        self,
        building_blocks: List[BaseBuildingBlock],
        library_id: Optional[str] = None,
    ):
        """
        Initialize a synthon

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

        # _last_non_null_bb = max([i for i, bb in enumerate(self.building_blocks)
        # if not bb.is_null()])
        # _first_null_bb = max([i for i, bb in enumerate(self.building_blocks)
        # if bb.is_null()])
        #
        # if _first_null_bb <= _last_non_null_bb:
        #     raise SynthonError(f"null building blocks can only occur after all
        #     non-null building blocks;"
        #                            f" first null bb at cycle {_first_null_bb},"
        #                            f" last non-null bb at cycle {_last_non_null_bb}")

        self.synthon_id = "-".join([building_block.bb_id for building_block in building_blocks])

        if self.library_id is not None:
            self.synthon_id = self.library_id + "-" + self.synthon_id

        self.num_cycles = len(self.building_blocks)

    def __contains__(self, item):
        """Return true if synthon contains a given building block"""
        if isinstance(item, BaseBuildingBlock):
            return any([_bb == item for _bb in self.building_blocks])
        return False

    def __eq__(self, other):
        """Return true if two synthons share the same ID"""
        if isinstance(other, Synthon):
            return self.synthon_id == other.synthon_id
        return False


class Monosynthon(Synthon):
    """
    Monosynthons are synthons with 1 Building block retained and all others masked out
    """

    def __init__(self, building_blocks: List[BaseBuildingBlock], library_id: Optional[str] = None):
        """
        Initialize a Monosynthon

        Parameters
        ----------
        building_blocks: List[BaseBuildingBlock]
            the building blocks that make up the synthon
            should be MaskedBuildingBlock for cycles were the BB info isn't needed
        library_id: str, optional
            the library ID for the DEL the synthon came from
            only necessary if BuildingBlock ids are not unique across libraries

        Raises
        ------
        SynthonError
            if Monosynthon has more than 1 real building block
        """
        super().__init__(building_blocks, library_id)
        num_real_bbs = sum([building_block.is_real() for building_block in self.building_blocks])
        if num_real_bbs != 1:
            raise SynthonError(
                f"Monosynthons can only have 1 real building_block, found {num_real_bbs}"
            )

        self.mono_bb = self.building_blocks[self.real_cycles[0]]

    def select(
        self,
        selection_count: Optional[float],
        control_count: Optional[float],
        competitive_count: Optional[float],
    ) -> None:
        """
        Add selection data to this Monosynthon

        Parameters
        ----------
        selection_count: float, optional
            added a selection condition count
        control_count: float, optional
            added a control condition count
        competitive_count: float, optional
            added a competitive condition count
        """
        pass


class HasMonosynthonMixin:
    """
    A mixin to give an object the ability to be decomposed into Monosynthons
    """

    library_id: Optional[str]
    building_blocks: List[BaseBuildingBlock]
    real_cycles: List[int]

    def get_monosynthons(self, skip: Optional[List[int]] = None) -> List[Monosynthon]:
        """
        Decompose object into all possible Monosynthons

        Parameters
        ----------
        skip: List[int], optional
            if passed, the cycle idx that you want to exclude from
            using when making monosynthons

        Returns
        -------
        List[Monosynthon]
        """
        _monosynthons = []
        if skip is None:
            skip = []
        for i in self.real_cycles:
            if i in skip:  # if asked to skip this cycle skip it
                continue
            if not self.building_blocks[i].is_real():  # skip bbs that aren't real
                continue
            _bbs: List[BaseBuildingBlock] = [MaskedBuildingBlock() for _ in self.building_blocks]
            _bbs[i] = deepcopy(self.building_blocks[i])
            _monosynthons.append(Monosynthon(building_blocks=_bbs, library_id=self.library_id))
        return _monosynthons

    def __contains__(self, item: object) -> bool:
        """Return True if object contains the passed BaseBuildingBlock OR Monosynthon"""
        if isinstance(item, BaseBuildingBlock):
            return any([_bb == item for _bb in self.building_blocks])
        else:
            if isinstance(item, Monosynthon):
                return (self.library_id == item.library_id) and (
                    item.mono_bb in self.building_blocks
                )
            return False


class Disynthon(Synthon, HasMonosynthonMixin):
    """
    Disynthons are synthons with 2 Building blocks retained and all others masked out
    """

    def __init__(self, building_blocks: List[BaseBuildingBlock], library_id: Optional[str] = None):
        """
        Initialize a Disynthon

        Parameters
        ----------
        building_blocks: List[BaseBuildingBlock]
            the building blocks that make up the synthon
            should be MaskedBuildingBlock for cycles were the BB info isn't needed
        library_id: str, optional
            the library ID for the DEL the synthon came from
            only necessary if BuildingBlock ids are not unique across libraries

        Raises
        ------
        SynthonError
            if Disynthon does not have 2 real building blocks
        """
        super().__init__(building_blocks, library_id)
        num_real_bbs = sum([building_block.is_real() for building_block in self.building_blocks])
        if num_real_bbs != 2:
            raise SynthonError(
                f"Disynthons can only have 2 real building_block, found {num_real_bbs}"
            )

        self.di_bbs = [self.building_blocks[i] for i in self.real_cycles]

    def select(
        self,
        selection_count: Optional[float],
        control_count: Optional[float],
        competitive_count: Optional[float],
    ) -> None:
        """
        Add selection data to this Disynthon

        Parameters
        ----------
        selection_count: float, optional
            added a selection condition count
        control_count: float, optional
            added a control condition count
        competitive_count: float, optional
            added a competitive condition count
        """
        pass


class HasDisynthonMixin(HasMonosynthonMixin):
    """Mixin for adding functionality to decompose into disynthon(s)"""

    def get_disynthons(self, skip: Optional[List[int]] = None) -> List[Disynthon]:
        """
        Decompose object into all possible disynthons

        Parameters
        ----------
        skip: List[int], optional
            if passed, the cycle idx that you want to exclude from
            using when making monosynthons

        Returns
        -------
        List[Disynthon]
        """
        _disynthons = []
        if skip is None:
            skip = []
        for i, j in combinations(self.real_cycles, 2):
            if (i in skip) or (j in skip):  # if asked to skip this cycle skip it
                continue
            # skip bbs that aren't real
            if (not self.building_blocks[i].is_real()) and (not self.building_blocks[j].is_real()):
                continue
            _bbs: List[BaseBuildingBlock] = [MaskedBuildingBlock() for _ in self.building_blocks]
            _bbs[i] = deepcopy(self.building_blocks[i])
            _bbs[j] = deepcopy(self.building_blocks[j])

            _disynthons.append(Disynthon(building_blocks=_bbs, library_id=self.library_id))
        return _disynthons

    def __contains__(self, item: object) -> bool:
        """Return True if object contains the passed BaseBuildingBlock OR Mono/Disynthon"""
        if item in self:
            return True
        else:
            if isinstance(item, Disynthon):
                return (self.real_cycles == item.synthon_id) and (
                    all([di in self.building_blocks for di in item.di_bbs])
                )
            return False
