"""classes for enumerating combinatorial libraries based on reaction trees"""
import warnings
from itertools import product as iter_product
from typing import Sequence, overload, Literal, Iterator, Optional, TYPE_CHECKING

from rdkit.Chem import Mol
from tqdm import tqdm

from deli.dels.building_block import BuildingBlockSet, BuildingBlock
from deli.dels.compound import DELCompound, DELCompoundRaw
from deli.utils import SmilesMixin

from .reaction import ReactionTree, ReactionVial, ReactionError

if TYPE_CHECKING:
    from deli.dels.library import Library


class EnumerationRunError(Exception):
    """error raised when an enumeration fails"""

    pass


class EnumeratedDELCompoundRaw(DELCompoundRaw, SmilesMixin):
    """
    Low memory DEL enumerated compound class

    This class is used to represent a DEL compound that has been fully enumerated
    and has a SMILES string mapped to the full compound, represented as string IDs.
    It inherits from DELCompoundRaw and implements the SmilesMixin.

    If Mol is not provided, it will be generated from the SMILES string
    and cached when first accessed.
    """

    def __init__(
        self,
        library_id: str,
        building_blocks_ids: list[str],
        smiles: str,
        mol: Optional[Mol] = None,
    ):
        """
        Initialize the EnumeratedDELCompoundRaw object.

        Parameters
        ----------
        library_id: str
            The ID of the library.
        building_blocks_ids: list[str]
            The IDs of the building blocks.
        smiles: str
            The SMILES string of the compound.
        mol: Optional[Mol], default=None
            The RDKit Mol object for the compound.
            If None, it will be generated from the SMILES string and
            cached when first accessed.
            NOTE: it is not recommended to set this directly unless
                  you are confident the mol object
                  is the results of the provided SMILES
        """
        super().__init__(library_id=library_id, building_blocks_ids=building_blocks_ids)
        self._smiles = smiles
        self._mol = mol


class EnumeratedDELCompound(DELCompound, SmilesMixin):
    """
    DEL enumerated compound class

    This class is used to represent a DEL compound that has been fully enumerated
    and has a SMILES string mapped to the full compound.
    It inherits from DELCompound and implements the SmilesMixin.

    If Mol is not provided, it will be generated from the SMILES string
    and cached when first accessed.
    """

    def __init__(
        self,
        library: "Library",
        building_blocks: list["BuildingBlock"],
        smiles: str,
        mol: Optional[Mol] = None,
    ):
        """
        Initialize the EnumeratedDELCompound object.

        Parameters
        ----------
        library: DELibrary
            The DELibrary object that this compound belongs to.
        building_blocks: dict[str, BuildingBlock]
            The building blocks of the compound
            should be in cycle order (cycle 1 first, then cycle 2, etc.)
        smiles: str
            The SMILES string of the compound.
        mol: Optional[Mol], default=None
            The RDKit Mol object for the compound.
            If None, it will be generated from the SMILES string
            and cached when first accessed.
            NOTE: it is not recommended to set this directly
                  unless you are confident the mol object
                  is the results of the provided SMILES
        """
        super().__init__(library=library, building_blocks=building_blocks)
        self._smiles = smiles
        self._mol = mol


class Enumerator:
    def __init__(self, building_block_sets: Sequence[BuildingBlockSet],  reaction_tree: ReactionTree):
        self.building_block_sets = building_block_sets
        self.reaction_tree = reaction_tree

    def get_bb_set(self, bb_set_id: str) -> BuildingBlockSet:
        """
        Get a building block set by its ID
        Parameters
        ----------
        bb_set_id: str
            The ID of the building block set to retrieve
        Returns
        -------
        BuildingBlockSet
            The building block set with the specified ID
        Raises
        ------
        KeyError
            If no building block set with the specified ID exists
        """
        for bb_set in self.building_block_sets:
            if bb_set.bb_set_id == bb_set_id:
                return bb_set
        raise KeyError(f"No building block set with ID {bb_set_id} found")

    def get_enumeration_size(self) -> int:
        """
        Get the total number of possible compounds in the enumeration

        This could be a different number than just a product of the building block set sizes,
        because not all building block sets may be used in the reaction tree.

        Returns
        -------
        int
        """
        total_size = 0
        for rxn_thread in self.reaction_tree.threads:
            thread_size = 1
            for bb_set_id, bb_set_reactant_id in zip(rxn_thread.bb_set_ids, rxn_thread.bb_set_reactant_ids):
                thread_size *= len(self.get_bb_set(bb_set_id).get_bb_subset(bb_set_reactant_id))
            total_size += thread_size
        return total_size

    @overload
    def enumerate(
        self, dropped_failed: Literal[True], fail_on_error: bool, use_tqdm: bool
    ) -> Iterator[tuple[list[BuildingBlock], Mol]]: ...

    @overload
    def enumerate(
        self, dropped_failed: bool, fail_on_error: Literal[True], use_tqdm: bool
    ) -> Iterator[tuple[list[BuildingBlock], Mol]]: ...

    @overload
    def enumerate(
        self, dropped_failed: Literal[False], fail_on_error: Literal[False], use_tqdm: bool
    ) -> Iterator[tuple[list[BuildingBlock], Mol | None]]: ...

    @overload
    def enumerate(
        self, dropped_failed: bool, fail_on_error: bool, use_tqdm: bool
    ) -> Iterator[tuple[list[BuildingBlock], Mol | None]]: ...

    def enumerate(
        self, dropped_failed: bool = False, fail_on_error: bool = False, use_tqdm: bool = False
    ) -> Iterator[tuple[list[BuildingBlock], Mol | None]]:
        """
        Run reactions for all possible products of the building block set

        The default behavior is to yield all compounds,
        whether they failed enumeration or not.
        This behavior can be controlled with the `dropped_failed` and `fail_on_error` parameters.

        Parameters
        ----------
        dropped_failed: bool, default = False
            skip yielding compounds that failed enumeration
            will only be utilized if `fail_on_error` is False
            otherwise an error will be raised if any compound fails enumeration
        fail_on_error: bool, default = False
            if True, will raise an error if ANY compound fails enumeration
            Unless you are confident all compounds will enumerate successfully,
            better to leave as False and instead yield the failed compounds without SMILES
        use_tqdm : bool
            turn on tqdm progress bar

        Yields
        ------
        tuple[list[BuildingBlock], Mol | None]
            A tuple of (building_block_ids, enumerated_mol)
            if dropped_failed is False, will yield Mol, not Mol | None

        Raises
        ------
        EnumerationRunError
            if `fail_on_error` is `True` and enumeration fails for any compound
            will raise an error with the failed compound's bb_id_map

        Warnings
        --------
        If `dropped_failed` is set to True and `fail_on_error` is also set to True,
        `fail_on_error` will take precedence
        """
        if dropped_failed and fail_on_error:
            warnings.warn(
                "'dropped_failed' is set to True, "
                "but 'fail_on_error' is also set to True; "
                "ignoring 'dropped_failed' and will raise "
                "an error if any compound fails enumeration",
                stacklevel=1,
            )

        # loop through all the threads in the reaction tree
        p_bar = tqdm(total=self.get_enumeration_size(), disable=not use_tqdm)
        for rxn_thread in self.reaction_tree.threads:
            bb_set_dict = {
                bb_set_reactant_id: self.get_bb_set(bb_set_id).get_bb_subset(bb_set_reactant_id)
                for bb_set_id, bb_set_reactant_id in zip(rxn_thread.bb_set_ids, rxn_thread.bb_set_reactant_ids)
            }
            bb_combos: iter_product[tuple[tuple[str, BuildingBlock], ...]] = iter_product(
                *[[(bb_set_id, bb) for bb in bb_set] for bb_set_id, bb_set in bb_set_dict.items()]
            )

            for bb_combo in bb_combos:
                try:
                    enumerated_mol = rxn_thread.run_thread(
                        vail=ReactionVial({bb[0]: bb[1].mol for bb in bb_combo})
                    )
                    yield [bb[1] for bb in bb_combo], enumerated_mol
                except ReactionError as e:
                    if fail_on_error:
                        raise EnumerationRunError(
                            f"failed to enumerate compound with building blocks: "
                            f"{[bb[0] for bb in bb_combo]}"
                        ) from e
                    elif dropped_failed:
                        continue
                    else:
                        yield [bb[1] for bb in bb_combo], None
                p_bar.update(1)

    def enumerate_by_bbs(self, bbs: list[BuildingBlock]) -> Mol:
        """
        Enumerate a single compound by providing the building blocks

        Parameters
        ----------
        bbs: list[BuildingBlock]
            list of building block ids to enumerate
            The should match the order of the building block sets in the library

        Returns
        -------
        EnumeratedDELCompound
            the enumerated compound with the provided building blocks

        Raises
        ------
        EnumerationRunError
            if enumeration fails
        """

        try:
            bb_subset_id_map = {bb_set.get_subset_with_bb(bb, as_bb_subset_id=True): bb.mol for bb_set, bb in
                                zip(self.building_block_sets, bbs)}
            rxn_thread = self.reaction_tree.get_thread_for_bb_subset_ids(frozenset(bb_subset_id_map.keys()))
        except (KeyError, ReactionError) as e:
            raise EnumerationRunError(
                f"failed to enumerate compound with building blocks: "
                f"{[bb.bb_id for bb in bbs]}"
            ) from e
        return rxn_thread.run_thread(ReactionVial(bb_subset_id_map))

    def enumerate_by_bb_ids(self, bb_ids: list[str]) -> Mol:
        """
        Enumerate a single compound by providing the building block IDs

        Parameters
        ----------
        bb_ids: list[str]
            list of building block ids to enumerate
            The should match the order of the building block sets in the library

        Returns
        -------
        Mol
            the enumerated Mol with the provided building blocks

        Raises
        ------
        EnumerationRunError
            if enumeration fails
        """
        bbs = []
        for bb_set, bb_id in zip(self.building_block_sets, bb_ids):
            try:
                bbs.append(bb_set.get_bb_by_id(bb_id, fail_on_missing=True))
            except KeyError as e:
                raise EnumerationRunError(
                    f"failed to find building block {bb_id} in building block set {bb_set.bb_set_id}"
                ) from e
        return self.enumerate_by_bbs(bbs)
