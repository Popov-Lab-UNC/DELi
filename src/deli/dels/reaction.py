"""define reaction classes"""

from typing import List, Tuple

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from deli.utils.mol_utils import Molable, to_mol


class ReactionError(Exception):
    """error for generic reaction failure"""

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

        if self.rxn.GetNumReactants() != 1:
            raise ReactionError(
                "reaction must produce exactly one product; " "found {}".format(
                    self.rxn.GetNumReactants()
                )
            )

    def react(self, *args: Tuple[Chem.Mol, ...]) -> Chem.Mol:
        """
        Run reaction on input molecules/SMILES

        Parameters
        ----------
        *args: Tuple[Chem.Mol, ...]
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

        # Take the first outcome's product
        product = products[0][0]

        # Sanitize the product before return
        try:
            return Chem.SanitizeMol(product)
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


class ReactionStep:
    """
    Contains required information to carry out a reaction step

    Reaction steps are abstract guidelines for how a reaction is carried out.
    It include the reaction itself, and info about the reactants that should be used.

    This separates "static" reactants from "variable" ones.
    "variable" reactants can have chemical structures vary between ReactantVials,
    and are instead identified by their reactant ID (from the ReactantVial).
    "static" reactants are always the same chemical.
    These are identified by their SMILES rather than an ID

    Attributes
    ----------
    product_id: str
        the id of the product produced by the reaction
        by default this is always "product_<step_id>"
    requires: List[str]
        the ids of the reactants that are necessary for the reaction
        same as the "variable_reactions" parameter that is passed
    """

    def __init__(
        self,
        step_id: int,
        reaction: Reaction,
        variable_reactant: List[str],
        static_reactants: List[Molable],
    ):
        """
        initialize a reaction step

        Parameters
        ----------
        step_id: int
            the id of the step
        reaction: Reaction
            the reaction to be carried out in the step
        variable_reactant: List[str]
            the ids of the reactants to be used during the step
        static_reactants: List[Molable]
            any chemicals that are always present during the reaction
        """
        self.step_id = step_id
        self.reaction = reaction
        self.requires = set(variable_reactant)
        self.static_reactants = [
            to_mol(static_reactant, fail_on_error=True) for static_reactant in static_reactants
        ]
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
        reactants = [possible_reactants[req] for req in self.requires] + self.static_reactants
        possible_reactants.add_product(self.product_id, self.reaction.react(*reactants))


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
