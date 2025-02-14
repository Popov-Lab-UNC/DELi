"""defines DEL enumerator functions and classes"""

import json
from collections.abc import Iterator
from functools import reduce
from itertools import product as iter_product
from operator import mul
from pathlib import Path
from typing import Callable, List, Optional, Union

import pandas as pd
from rdkit import Chem
from tqdm import tqdm

from deli.utils.mol_utils import to_smi

from ..configure import DeliDataLoadable, accept_deli_data_name
from .building_block import BuildingBlock, BuildingBlockSet
from .reaction import ReactionError, ReactionVial, ReactionWorkflow


FAILED_ENUMERATION_STR = "ENUMERATION_FAILED"


class EnumeratedDELCompound:
    """
    A fully enumerated DEL compound

    Will have the rdkit Mol of the chemical, the SMILES
    and all teh building blocks in at as attributes
    named after the building block set id
    """

    def __init__(self, mol: Chem.rdchem.Mol, building_block_id_map: dict[str, BuildingBlock]):
        """
        initializes the enumerated DEL compound

        Parameters
        ----------
        mol: Chem.Mol
            the enumerated compounds as a mol object
        building_block_id_map: dict[str, BuildingBlock]
            the mapping from building block set IDs to building blocks in the compound
        """
        self.mol: Chem.rdchem.Mol = mol
        for bb_set_id, bb_id in building_block_id_map.items():
            self.__setattr__(bb_set_id, bb_id)
        self.smi: str = to_smi(self.mol)


class FailedEnumeratedDELCompound(EnumeratedDELCompound):
    """
    special case of `EnumeratedDELCompound` when enumeration fails

    Will set the `smi` to `FAILED_ENUMERATION_STR` and `mol` to an empty mol object
    """

    def __init__(self, building_block_id_map: dict[str, BuildingBlock]):
        """
        Initializes the failed enumerated DEL compound

        Parameters
        ----------
        building_block_id_map: dict[str, BuildingBlock]
            the mapping from building block set IDs to building blocks in the compound
        """
        super().__init__(Chem.Mol(), building_block_id_map)
        self.smi: str = FAILED_ENUMERATION_STR


class EnumeratorBuildError(Exception):
    """error raised when a enumerator build fails"""

    pass


class EnumeratorRunError(Exception):
    """error raised when a enumerator run fails"""

    pass


class DELEnumerator(DeliDataLoadable):
    """Handles DNA-encoded library enumeration using library reaction workflows"""

    def __init__(
        self,
        reaction_workflow: ReactionWorkflow,
        bb_sets: List[BuildingBlockSet],
        scaffold: Optional[str] = None,
    ):
        """
        Initialize enumerator with building blocks and reaction workflow.

        Parameters
        ----------
        bb_sets : List[BuildingBlockSet]
            Available building block sets
        reaction_workflow : ReactionWorkflow
            Reaction workflow definition
        scaffold : Optional[str]
            Optional scaffold SMILES string
        """
        self.bb_sets: dict[str, BuildingBlockSet] = {
            bb_set.bb_set_id: bb_set for bb_set in bb_sets
        }
        self.reaction_workflow = reaction_workflow
        self.scaffold = scaffold

        ### validation checks ###

        # check that there are SMILES
        if not all([bb_set.has_smiles for bb_set in self.bb_sets.values()]):
            for bb_set in self.bb_sets.values():
                if not bb_set.has_smiles:
                    raise EnumeratorBuildError(
                        "BuildingBlockSet '{}' lacks SMILES for all compounds".format(
                            bb_set.bb_set_id
                        )
                    )

        # check that the enumerator reaction workflow requires the bb_sets used
        _required_bb_set_ids = set(
            [
                _
                for step in self.reaction_workflow.rxn_steps
                for _ in step.requires
                if not _.startswith("product_")
            ]
        )
        _found_bb_set_ids = set(list(self.bb_sets.keys()))

        if (_required_bb_set_ids - _found_bb_set_ids) != set():
            raise EnumeratorBuildError(
                f"missing required building block set(s) "
                f"'{_required_bb_set_ids - _found_bb_set_ids}'"
            )
        if (_found_bb_set_ids - _required_bb_set_ids) != set():
            raise EnumeratorBuildError(
                f"found extra building block set(s) "
                f"'{_found_bb_set_ids - _required_bb_set_ids}'"
            )

        # check that scaffold is not None is scaffold is required for reaction
        if "scaffold" in _required_bb_set_ids and self.scaffold is None:
            raise EnumeratorBuildError(
                "enumeration reaction requires a scaffold but no scaffold was specified"
            )

    @classmethod
    @accept_deli_data_name("libraries", "json")
    def load(cls, path: str):
        """Load in a enumerator from a library json file"""
        return cls.from_json(path)

    @classmethod
    def from_json(cls, path: str):
        """Load in a enumerator from a library json file"""
        data = json.load(open(path))

        # load bb sets (needed for reaction setup)
        bb_sets: list[BuildingBlockSet] = [BuildingBlockSet.load(bb) for bb in data["bb_sets"]]
        bb_set_ids = set([bb_set.bb_set_id for bb_set in bb_sets] + ["scaffold"])

        reaction_workflow = ReactionWorkflow.load_from_json_list(data["reactions"], bb_set_ids)

        return cls(
            bb_sets=bb_sets,
            reaction_workflow=reaction_workflow,
            scaffold=data.get("scaffold"),
        )

    def __repr__(self) -> str:
        """Represent the enumerator as number of BB sets and reactions in it"""
        return (
            f"<DELEnumerator: {len(self.bb_sets)} BB sets, "
            f"{len(self.reaction_workflow)} reactions>"
        )

    def enumerate(self, use_tqdm: bool = False) -> Iterator[EnumeratedDELCompound]:
        """
        Run reactions for all possible products of the building block set

        Parameters
        ----------
        use_tqdm : bool
            turn on tqdm progress bar

        Yields
        ------
        EnumeratedDELCompound
            the enumerated DEL compound
            if enumeration failed will be a special child of
            EnumeratedDELCompound called FailedEnumeratedDELCompound
        """
        bb_combos: iter_product[tuple[tuple[str, BuildingBlock], ...]] = iter_product(
            *[[(key, bb) for bb in val.building_blocks] for key, val in self.bb_sets.items()]
        )
        num_bb_combos = reduce(mul, [len(bb_set) for bb_set in self.bb_sets.values()])
        for bb_combo in tqdm(bb_combos, total=num_bb_combos, disable=not use_tqdm):
            bb_id_map = {bb[0]: bb[1] for bb in bb_combo}
            try:
                enumerated_mol = self.reaction_workflow.run_workflow(
                    ReactionVial({bb[0]: bb[1].mol for bb in bb_combo})
                )
                yield EnumeratedDELCompound(enumerated_mol, bb_id_map)
            except ReactionError:
                yield FailedEnumeratedDELCompound(bb_id_map)

    def __iter__(self) -> Iterator[EnumeratedDELCompound]:
        """Enumerate all the compounds; alias for self.enumerate()"""
        return self.enumerate()

    def get_enumerated_compound_from_bb_ids(self, building_block_id_map: dict[str, str]):
        """
        Given a mapping of bb_set_ids to bb_ids, enumerate the single compound

        Parameters
        ----------
        building_block_id_map: dict[str, str]
            mapping of bb_set_ids to bb_ids for the compound to enumerate

        Returns
        -------
        EnumeratedDELCompound
            the enumerated DEL compound
            if enumeration failed will be a special child of
            EnumeratedDELCompound called FailedEnumeratedDELCompound
        """
        vial = ReactionVial()
        bb_set_id_map: dict[str, BuildingBlock] = dict()
        for bb_set_id, bb_id in building_block_id_map.items():
            try:
                bb = self.bb_sets[bb_set_id].get_bb_by_id(bb_id, fail_on_missing=True)
            except KeyError as e:
                raise EnumeratorRunError(
                    f"unrecognized building block set id '{bb_set_id}'"
                ) from e

            if bb is None:
                raise EnumeratorRunError(
                    f"building block id '{bb_set_id}' not found in BB set '{bb_id}'"
                )
            bb_set_id_map[bb_set_id] = bb
            vial[bb_set_id] = bb.mol

        try:
            return EnumeratedDELCompound(self.reaction_workflow.run_workflow(vial), bb_set_id_map)
        except ReactionError:
            return FailedEnumeratedDELCompound(bb_set_id_map)

    def enumerate_to_csv_file(
        self,
        out_path: Union[str, Path],
        compound_id_function: Optional[Callable] = None,
        use_tqdm: bool = False,
    ) -> None:
        """
        Write the enumerated compounds to a CSV file.

        Writing happens on the fly, so file will be written to as
        the enumeration progresses

        File header will be "SMILES,[BB_SET_IDS]"
        If compound_id_function is not None, "CompoundID" column will
        appear before "SMILES" column

        Parameters
        ----------
        out_path: Union[str, Path]
            the path to the csv file to write
        compound_id_function: Optional[Callable]
            a function that takes a mapping of bb_set_ids to BuildingBlock objects
            and generate a compound id from the mapping
            if left as None will exclude the "CompoundID" column
        use_tqdm: bool, default False
            whether to use tqdm progress bar
        """
        bb_sets_ids = list(self.bb_sets.keys())
        header = "CompoundID," if compound_id_function else ""
        header += "SMILES," + ",".join(bb_sets_ids) + "\n"
        with open(out_path, "w") as f:
            f.write(header)
            for enumerated_compound in self.enumerate(use_tqdm=use_tqdm):
                row = (
                    f"{compound_id_function(enumerated_compound.__dict__)},"
                    if compound_id_function
                    else ""
                )
                row += f"{enumerated_compound.smi}"
                for bb_sets_id in bb_sets_ids:
                    row += f",{str(enumerated_compound.__getattribute__(bb_sets_id))}"
                row += "\n"
                f.write(row)

    def enumerate_to_pandas(self) -> pd.DataFrame:
        """
        Enumerate all compounds into pandas DataFrame

        WARNING: this could be extremely memory hungry for large libraries, use with caution

        Columns will be
        - "smi" with smiles
        - "mol" with rdkit mol object of compound
        - "<BB_SET_ID>" with the corresponding BuildingBlock object for that set

        Returns
        -------
        pd.DataFrame
        """
        return pd.DataFrame([comp.__dict__ for comp in self.enumerate()])
