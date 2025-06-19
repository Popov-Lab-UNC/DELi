"""define reaction classes"""

import abc
import warnings
from typing import List, Tuple

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from deli.utils.mol_utils import to_mol


class ReactionError(Exception):
    """error for generic reaction failure"""

    pass


class ReactionWarning(UserWarning):
    """warning for abnormal reaction behavior that is not a failure"""

    pass


class Reaction:
    """Class for describing a singular chemical reaction"""

    def __init__(self, rxn_smarts: str):
        """
        Initialize the reaction object

        Parameters
        ----------
        rxn_smarts: str
            SMARTS or SMIRKS to define a reaction
        """
        self.rxn = rdChemReactions.ReactionFromSmarts(rxn_smarts)
        self.num_reactants = self.rxn.GetNumReactantTemplates()

        if self.rxn.GetNumProductTemplates() != 1:
            raise ReactionError(
                "reaction must produce exactly one product; found {}".format(
                    self.rxn.GetNumProductTemplates()
                )
            )

    def _order_reactants(self, *args: Chem.Mol) -> Tuple[Chem.Mol, ...]:
        """RDKit needs reactants in the right order, this will do that"""
        matches: List[List[Chem.Mol]] = list()
        for i in range(self.num_reactants):
            react_matches: List[Chem.Mol] = list()
            _react_temp = self.rxn.GetReactantTemplate(i)
            for mol in args:
                if mol.HasSubstructMatch(_react_temp):
                    react_matches.append(mol)
            matches.append(react_matches)

        reactant_tuple: List[Chem.Mol] = []
        for i, match_set in enumerate(matches):
            if len(match_set) > 1:
                warnings.warn(
                    f"found multiple matches for reactant template {i} in reaction "
                    f"{rdChemReactions.ReactionToSmarts(self.rxn)}. "
                    f"selecting first match by default; could cause invalid product",
                    ReactionWarning,
                    stacklevel=1,
                )
            elif len(match_set) == 0:
                raise ReactionError(
                    f"failed to match reactants to reactant template {i} in reaction "
                    f"{rdChemReactions.ReactionToSmarts(self.rxn)}"
                )
            reactant_tuple.append(match_set[0])

        return tuple(reactant_tuple)

    def react(self, *args: Chem.Mol) -> Chem.Mol:
        """
        Run reaction on input molecules/SMILES

        Parameters
        ----------
        *args: Chem.Mol
            the molecules to use as reactants

        Returns
        -------
        product: Chem.Mol
            the product of the reaction
        """
        # Validate input count
        if len(args) != self.num_reactants:
            raise ReactionError(f"Expected {self.num_reactants} reactants, got {len(args)}")

        # Run reaction
        products = self.rxn.RunReactants(args)
        if not products:
            raise ReactionError("Failed to generate product from reactants")

        if len(products) > 1:
            reactant_smiles = [Chem.MolToSmiles(arg) for arg in args]
            warnings.warn(
                f"multiple reaction pathways found for reaction "
                f"{rdChemReactions.ReactionToSmarts(self.rxn)} with "
                f"reactants {reactant_smiles}, "
                f"selecting first product by default",
                ReactionWarning,
                stacklevel=1,
            )

        # Take the first outcome's product
        product = products[0][0]

        # Sanitize the product before return
        try:
            product.UpdatePropertyCache()
            return product
        except Exception as e:
            raise RuntimeError(f"Failed to sanitize reaction product: {e}") from e


class ReactionVial(dict[str, Chem.Mol]):
    """holds the chemicals that are current exist as a reaction progresses"""

    def add_product(self, product_id: str, product: Chem.Mol):
        """
        Add a product to the reaction "soup"

        Notes
        -----
        This is no different than just assigning the key value pair like a normal dict.
        This is just to help with readability

        Parameters
        ----------
        product_id: str
            the id of the product
        product: Chem.Mol
            the product molecule
        """
        self[product_id] = product


class Reactant(abc.ABC):
    """abstract class for all reactants"""

    @abc.abstractmethod
    def get_from_vial(self, vial: ReactionVial) -> Chem.Mol:
        """Return a mol of the reactant from the pass vial"""
        raise NotImplementedError


class StaticReactant(Reactant):
    """A reactant that is the same for all DEL collections"""

    def __init__(self, mol: Chem.Mol):
        """
        Initialize a static reactant

        Parameters
        ----------
        mol: Chem.Mol
            the molecule of the static reactant
        """
        self.mol = mol

    def get_from_vial(self, vial: ReactionVial) -> Chem.Mol:
        """
        Returns the static molecule

        The vial is ignored, as static reactants never change.
        Argument is for compatability

        Parameters
        ----------
        vial: ReactionVial
            ignored; for compatability

        Returns
        -------
        Chem.Mol
        """
        return self.mol


class BBSetReactant(Reactant):
    """A reactant that varies based on the BB_ID from the BB_Set it comes from"""

    def __init__(self, bb_set_id: str):
        """
        Initialize a BBSet reactant

        Parameters
        ----------
        bb_set_id: the BuildingBlockSet ID for the reactant to be pulled from
        """
        self.bb_set_id = bb_set_id

    def get_from_vial(self, vial: ReactionVial) -> Chem.Mol:
        """
        Returns the variable molecule from the BB set with matching ID from the vial

        Parameters
        ----------
        vial: ReactionVial
            vial to search in

        Returns
        -------
        Chem.Mol

        Raises
        ------
        ReactionError
            if there is no matching BB set id in the reaction vial
        """
        reactant = vial.get(self.bb_set_id, None)
        if reactant is None:
            raise ReactionError(
                f"failed to find a reactant from BB Set {self.bb_set_id} "
                f"in current reactant set {list(vial.keys())}"
            )
        return reactant


class ReactionStep:
    """
    Contains required information to carry out a reaction step

    Reaction steps are abstract guidelines for how a reaction is carried out.
    It include the reaction itself, and info about the reactants that should be used.

    Attributes
    ----------
    product_id: str
        the id of the product produced by the reaction
        by default this is always "product_<step_id>"
    """

    def __init__(self, step_id: int, reaction: Reaction, reactants: List[Reactant]):
        """
        Initialize a reaction step

        Parameters
        ----------
        step_id: int
            the id of the step
        reaction: Reaction
            the reaction to be carried out in the step
        reactants: List[Reactants]
            the reactants to be used for this step
            must be passed in order that the rxn smarts uses them
            See https://github.com/rdkit/rdkit/issues/4361
        """
        self.step_id = step_id
        self.reaction = reaction
        self.reactants = reactants
        self.product_id = f"product_{self.step_id}"

    def run_step(self, possible_reactants: ReactionVial):
        """
        Run the reaction step

        This will update the ReactionVial inplace, adding the product to the vial with
        with the steps `product_id`

        Parameters
        ----------
        possible_reactants: ReactionVial
            the pool of possible reactants to be used during the step
        """
        possible_reactants.add_product(
            self.product_id,
            self.reaction.react(
                *[reactant.get_from_vial(possible_reactants) for reactant in self.reactants]
            ),
        )


class ReactionWorkflowError(Exception):
    """error raised if a reaction workflow fails (during build or execution)"""

    pass


class ReactionWorkflow:
    """
    Merge a series of ReactionSteps into a single workflow
    """

    def __init__(self, rxn_steps: List[ReactionStep]):
        """
        Initialize the reaction workflow

        Notes
        -----
        The passed reaction steps must have unique ids and grow
        sequentially from 1 (e.g. 1, 2, 3, 4, ...).
        `0` is not a valid ID.
        Steps do not need to be passed in order, they will be sorted based on their ID
        (lower ids come first).

        Parameters
        ----------
        rxn_steps: List[ReactionStep]
            the list of reaction steps
        """
        # check that step ids are unique and sequential
        _step_ids = sorted([step.step_id for step in rxn_steps])
        if list(range(1, len(_step_ids) + 1)) != _step_ids:
            raise ReactionWorkflowError(
                f"step ids must start from 1 and increase sequentially; "
                f"expected ids {list(range(1, len(_step_ids) + 1))}, "
                f"found ids {_step_ids}"
            )

        # sort the steps
        self.rxn_steps = sorted(rxn_steps, key=lambda x: x.step_id)
        self._final_product_id = f"product_{self.rxn_steps[-1].step_id}"

    def __len__(self) -> int:
        """Number of reaction steps in workflow"""
        return len(self.rxn_steps)

    def run_workflow(self, starting_reactants: ReactionVial) -> Chem.Mol:
        """
        Execute the reaction workflow for a given ReactionVial

        The ReactionVial should only contain the starting reactants
        pulled from the DEL Building Block libraries

        Parameters
        ----------
        starting_reactants: ReactionVial

        Returns
        -------
        product: Chem.Mol
            the final enumerated DEL
        """
        for step in self.rxn_steps:
            step.run_step(starting_reactants)
        return starting_reactants[self._final_product_id]

    @classmethod
    def load_from_json_list(cls, rxn_list: list[dict], bb_set_ids: set[str]) -> "ReactionWorkflow":
        """
        Load a reaction workflow from a list of reaction steps defined in json format

        This requires the building blocks sets ids to determine what is
        a static reactant

        Parameters
        ----------
        rxn_list: list[dict]
            the list of reaction steps defined in json dict format
        bb_set_ids: set[str]
            the set of building block set ids

        Returns
        -------
        ReactionWorkflow
        """
        rxn_steps: list[ReactionStep] = list()
        for rxn_data in rxn_list:
            step_id = rxn_data["step"]
            rxn_smarts = rxn_data["rxn_smarts"]
            reactant_ids = rxn_data["reactants"]

            bb_set_ids.add(f"product_{step_id}")  # to make it skip products

            # check if a reactant is static
            # static reactants are those that are valid SMILES and do not map to a bb_set id
            reactants: list[Reactant] = list()
            for reactant_id in reactant_ids:
                if reactant_id not in bb_set_ids:
                    reactants.append(StaticReactant(to_mol(reactant_id, fail_on_error=True)))
                else:
                    reactants.append(BBSetReactant(reactant_id))
            rxn_steps.append(
                ReactionStep(step_id=step_id, reaction=Reaction(rxn_smarts), reactants=reactants)
            )
        return cls(rxn_steps)
