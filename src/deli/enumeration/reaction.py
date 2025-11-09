"""define reaction classes"""

import abc
import warnings
from typing import Optional, Sequence, no_type_check

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from deli.dels.building_block import parse_building_block_subset_id
from deli.utils.mol_utils import to_mol, to_smi


class ReactionError(Exception):
    """error for generic reaction failure"""

    pass


class ReactionParsingError(Exception):
    """error for issues with reaction information"""

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

    def _order_reactants(self, *args: Chem.Mol) -> tuple[Chem.Mol, ...]:
        """RDKit needs reactants in the right order, this will do that"""
        matches: list[list[Chem.Mol]] = list()
        for i in range(self.num_reactants):
            react_matches: list[Chem.Mol] = list()
            _react_temp = self.rxn.GetReactantTemplate(i)
            for mol in args:
                if mol.HasSubstructMatch(_react_temp):
                    react_matches.append(mol)
            matches.append(react_matches)

        reactant_tuple: list[Chem.Mol] = []
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

    def to_smarts(self) -> str:
        """
        Return the reaction SMARTS/SMIRKS string

        Returns
        -------
        str
            the reaction SMARTS/SMIRKS
        """
        return rdChemReactions.ReactionToSmarts(self.rxn)

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
                f"{self.to_smarts()} with reactants {reactant_smiles}, "
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
    """Holds the chemicals that currently exist as a reaction progresses"""

    def add_product(self, product_id: str, product: Chem.Mol):
        """
        Add a product to the reaction "vial"

        Notes
        -----
        This is no different from just assigning the key value pair like a normal dict.
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

    def __str__(self):
        """String representation of the static reactant"""
        return to_smi(self.mol)

    def __repr__(self):
        """Representation of the static reactant"""
        return f"StaticReactant({to_smi(self.mol)})"

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

    def __init__(self, bb_set_reactant_id: str):
        """
        Initialize a BBSet reactant

        Parameters
        ----------
        bb_set_reactant_id: str
            This can be either the bb_set_id or a bb_subset_id
            for the building block set (or subset) that this reactant comes from
        """
        self.bb_set_reactant_id = bb_set_reactant_id

    def __str__(self):
        """String representation of the BB set reactant"""
        return self.bb_set_reactant_id

    def __repr__(self):
        """Representation of the BB set reactant"""
        return f"BBSetReactant({self.bb_set_reactant_id})"

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
        reactant = vial.get(self.bb_set_reactant_id, None)
        # error handling for easier debugging
        if reactant is None:
            try:
                bb_set_id, bb_subset_id = parse_building_block_subset_id(self.bb_set_reactant_id)
            except ValueError:
                bb_set_id = self.bb_set_reactant_id
                bb_subset_id = None
            _error_msg = (
                f"failed to find a reactant with id {self.bb_set_reactant_id} "
                f"from BB Set {bb_set_id}"
            )
            if bb_subset_id is not None:
                _error_msg += f" and subset {bb_subset_id}"
            _error_msg += f" in current reactant set {list(vial.keys())}"
            raise ReactionError(_error_msg)
        return reactant


class ProductReactant(Reactant):
    """A reactant that is the product of a previous reaction step"""

    def __init__(self, step_id: str):
        """
        Initialize a Product reactant

        Parameters
        ----------
        step_id: str
            the id of the reaction step that produced the product to be used as a reactant
        """
        self.step_id = step_id

    def __str__(self):
        """String representation of the product reactant"""
        return self.product_id

    def __repr__(self):
        """Representation of the product reactant"""
        return f"ProductReactant({self.product_id})"

    @property
    def product_id(self):
        """The product id derived from the step id"""
        return f"product_{self.step_id}"

    @product_id.setter
    def product_id(self, value):
        """Prevent modification of product_id property"""
        raise RuntimeError(
            "product_id is a read-only property derived from step_id; cannot be modified"
        )

    def get_from_vial(self, vial: ReactionVial) -> Chem.Mol:
        """
        Returns the product molecule from the vial

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
            if there is no matching product id in the reaction vial
        """
        reactant = vial.get(self.product_id, None)
        if reactant is None:
            raise ReactionError(
                f"failed to find a product from previous step {self.product_id[8:]} "
                f"in current reactant set {list(vial.keys())}; "
                f"make sure you are using the correct `step_id` (not the `step_name`). "
                f"See the reaction definition docs for me details."
            )
        return reactant


class PooledReactant(Reactant):
    """A reactant that is pooled from multiple sources"""

    def __init__(self, reactants: Sequence[Reactant]):
        """
        Initialize a pooled reactant

        Parameters
        ----------
        reactants: Sequence[Reactant]
            the list of reactants to pool
            the first available reactant from the pool list found in the vial will be used
        """
        self.reactants = reactants

        if any([isinstance(reactant, PooledReactant) for reactant in reactants]):
            raise ReactionParsingError("nested PooledReactants are not supported")

    def get_from_vial(self, vial: ReactionVial) -> Chem.Mol:
        """
        Returns the first available reactant from the vial

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
            if there is no matching reactant id in the reaction vial
        """
        for reactant in self.reactants:
            try:
                return reactant.get_from_vial(vial)
            except ReactionError:
                continue  # try the next reactant in the pool
        raise ReactionError(
            f"failed to find any of the pooled reactants {[str(r) for r in self.reactants]} "
            f"in current reactant set {list(vial.keys())}"
        )


class ReactionStep:
    """
    Contains required information to carry out a reaction step

    Reaction steps are abstract guidelines for how a reaction is carried out.
    Includes the reaction itself, and info about the reactants that should be used.

    Attributes
    ----------
    product_id: str
        the id of the product produced by the reaction
        by default this is always "product_<step_id>"
    """

    def __init__(
        self,
        step_name: str,
        reaction: Reaction | Sequence[Reaction],
        reactants: list[Reactant],
        step_id: Optional[str] = None,
    ):
        """
        Initialize a reaction step

        Parameters
        ----------
        step_id: int
            the id of the step
        reaction: Reaction | list[Reaction]
            the reaction to be carried out in the step
            if a list is passed, the first reaction that successfully
            produces a product will be used
        reactants: list[Reactants]
            the reactants to be used for this step
            must be passed in order that the rxn smarts uses them
            See https://github.com/rdkit/rdkit/issues/4361
        """
        self.step_name = step_name
        self.step_id: str = step_id if step_id else step_name
        self.reactions: Sequence[Reaction] = (
            [reaction] if isinstance(reaction, Reaction) else reaction
        )
        self.reactants = reactants
        self.product_id = f"product_{self.step_id}"

        # some validation to check reactant count matches reaction requirements:
        for rxn in self.reactions:
            if rxn.num_reactants != len(self.reactants):
                raise ReactionError(
                    f"reaction {rxn.to_smarts()} in step {self.step_id} "
                    f"expects {rxn.num_reactants} "
                    f"reactants, but {len(self.reactants)} were provided"
                )

    def flatten_reactants(self) -> list[Reactant]:
        """
        Return a flattened list of all reactants used in this step

        This will expand any PooledReactants into their individual reactants

        Returns
        -------
        list[Reactant]
            the flattened list of reactants
        """
        flattened_reactants: list[Reactant] = list()
        for reactant in self.reactants:
            if isinstance(reactant, PooledReactant):
                flattened_reactants.extend(reactant.reactants)
            else:
                flattened_reactants.append(reactant)
        return flattened_reactants

    def run_step(self, possible_reactants: ReactionVial):
        """
        Run the reaction step

        All reactions defined in the step will be attempted in the priority order
        they were passed until one successfully produces a product.

        Notes
        -----
        This will update the ReactionVial inplace, adding the product to the vial with the
        steps `product_id` attribute as the key.

        Parameters
        ----------
        possible_reactants: ReactionVial
            the pool of possible reactants to be used during the step

        Raises
        ------
        ReactionError
            if all reactions fail to produce a product
        """
        for reaction in self.reactions:
            try:
                possible_reactants.add_product(
                    self.product_id,
                    reaction.react(
                        *[
                            reactant.get_from_vial(possible_reactants)
                            for reactant in self.reactants
                        ]
                    ),
                )
                return  # success, exit the method
            except ReactionError:
                continue  # try the next reaction
        # if we reach here, all reactions failed
        raise ReactionError(f"all reactions failed for step {self.step_id}")


class ReactionThread:
    """
    A single linear path through a reaction tree

    It defines the ordered reaction steps that are needed to enumerate
    a given compound in a DEL library. Not all compounds from a DEL will
    utilize any given thread; instead they must be only compatible with
    a single reaction thread.

    Attributes
    ----------
    bb_set_ids: set[str]
        the set of building block set ids used by this reaction thread
    bb_set_reactant_ids: set[str]
        the set of building block set *or subset* ids used by this reaction thread
        if subsets are being used, this can be different from `bb_set_ids`
    """

    def __init__(self, reaction_steps: Sequence["ReactionStep"]):
        """
        Initialize a reaction thread

        Parameters
        ----------
        reaction_steps: tuple[_ReactionNode, ...]
            the ordered reaction steps that make up the thread
            *must* be in order from first to last step
        """
        self.reaction_steps = (
            reaction_steps if isinstance(reaction_steps, tuple) else tuple(reaction_steps)
        )
        self._final_product_id = self.reaction_steps[-1].product_id

        # validate that each BuildingBlockSet used by the thread only appears once
        self.bb_set_ids, self.bb_set_reactant_ids = self._validate_thread()

    def _validate_thread(self) -> tuple[list[str], list[str]]:
        """Validate the reaction thread for correctness"""
        # make sure thread uses each unique BB set once
        _bb_set_ids: list[str] = list()
        _bb_set_reactant_ids: list[str] = list()
        for step in self.reaction_steps:
            for reactant in step.flatten_reactants():
                if isinstance(reactant, BBSetReactant):
                    try:
                        bb_set_id, _ = parse_building_block_subset_id(reactant.bb_set_reactant_id)
                    except ValueError:  # not a building block subset id
                        bb_set_id = reactant.bb_set_reactant_id

                    if bb_set_id in _bb_set_ids:
                        raise ReactionParsingError(
                            f"building block set {bb_set_id} used multiple times "
                            f"in reaction thread {self._get_thread_description()}; "
                            f"each building block set can only "
                            f"be used once per reaction thread"
                        )
                    _bb_set_ids.append(bb_set_id)
                    _bb_set_reactant_ids.append(reactant.bb_set_reactant_id)
        return _bb_set_ids, _bb_set_reactant_ids

    def __len__(self) -> int:
        """Number of reaction steps in workflow"""
        return len(self.reaction_steps)

    def _get_thread_description(self) -> str:
        """
        Return a string description of the reaction thread

        Returns
        -------
        str
            the string description of the reaction thread
        """
        return " -> ".join([step.step_id for step in self.reaction_steps])

    def run_thread(self, starting_reactants: ReactionVial) -> Chem.Mol:
        """
        Run the reaction thread with the given starting reactants

        Returns
        -------
        product: Chem.Mol
            the final product of the reaction thread
        """
        for step in self.reaction_steps:
            step.run_step(starting_reactants)
        return starting_reactants[self._final_product_id]


class _ReactionNode:
    """Internal support class for building reaction trees"""

    def __init__(self, step: ReactionStep):
        self.step = step
        self.parents: list["_ReactionNode"] = list()
        self.children: list["_ReactionNode"] = list()

    def __eq__(self, other):
        """Reaction nodes are equal if their step ids are equal"""
        if not isinstance(other, _ReactionNode):
            return False
        return self.step.step_id == other.step.step_id


class ReactionTree:
    """
    Contains the full reaction scheme for a DEL library

    DEL reactions can involve complex branching steps where not all
    building blocks from a given cycle undergo the same reaction.

    DELi handles this by building a reaction tree, starting with the DNA tag linker
    and progressing down to the final product. Every possible path from the start to end
    is a valid reaction "thread" and specifies which building blocks will utilize this thread
    during enumeration.

    Notes
    -----
    The reaction tree is an acyclic-directed graph (DAG) where each node is a ReactionStep
    and each path from root to leaf is a ReactionThread. This is a limitation of chemistry,
    not of the algorithm. If a cycle is detected during building of the reaction tree,
    an error will be raised since this means the reaction scheme attempting to be parsed is
    not chemically valid.

    """

    def __init__(self, reaction_steps: list[ReactionStep]):
        """
        Initialize the reaction tree

        Parameters
        ----------
        reaction_steps: list[ReactionStep]
            the list of reaction steps that make up the reaction tree
        """
        self.reaction_steps = reaction_steps

        # build the tree
        self._roots, self._nodes = self._build_tree()

        # find all reaction threads in the tree
        self.threads = self._find_threads()

        # validate that all threads have different building block sets/subsets
        _bb_set_reactant_id_thread_sets: set[frozenset[str]] = set()
        for thread in self.threads:
            bb_set_reactant_id_thread_set = frozenset(thread.bb_set_reactant_ids)
            if bb_set_reactant_id_thread_set in _bb_set_reactant_id_thread_sets:
                raise ReactionParsingError(
                    f"multiple reaction threads use the same set of building block sets/subsets "
                    f"{thread.bb_set_reactant_ids}; each reaction thread must use a unique set "
                    f"of building block sets/subsets to avoid ambiguous synthesis"
                )
            _bb_set_reactant_id_thread_sets.add(bb_set_reactant_id_thread_set)

    def _build_tree(self) -> tuple[tuple[_ReactionNode, ...], list[_ReactionNode]]:
        _nodes: dict[str, _ReactionNode] = {
            step.step_id: _ReactionNode(step) for step in self.reaction_steps
        }

        for step_id, node in _nodes.items():
            for reactant in node.step.flatten_reactants():
                # product nodes give us all the info we will need
                if isinstance(reactant, ProductReactant):
                    parent_step_id = reactant.step_id  # no "product_" prefix
                    try:
                        parent_node = _nodes[parent_step_id]
                    except KeyError as e:
                        raise ReactionParsingError(
                            f"reaction step {step_id} references missing "
                            f"parent step {parent_step_id}"
                        ) from e
                    node.parents.append(parent_node)
                    parent_node.children.append(node)

        # check how many root nodes there are (nodes with no parents)
        root_nodes = [node for node in _nodes.values() if not node.parents]
        if len(root_nodes) == 0:
            raise ReactionParsingError(
                "no root node found in reaction tree; possible recation cycle "
                "detected (all step rely on previous steps)"
            )

        return tuple(root_nodes), list(_nodes.values())

    def _find_threads(self) -> tuple[ReactionThread, ...]:
        # do a DFS from each root node and save all possible paths to leaves
        paths: list[list[_ReactionNode]] = list()
        for root_node in self._roots:
            stack = [(root_node, [root_node])]
            while stack:
                current_node, path = stack.pop()
                if not current_node.children:
                    paths.append(path)
                else:
                    for child in current_node.children:
                        if child in path:
                            raise ReactionParsingError(
                                f"step {child.step.step_name} requires itself as a dependency; "
                                f"this is chemically impossible"
                            )
                        stack.append((child, path + [child]))
        return tuple(ReactionThread([node.step for node in path]) for path in paths)

    @no_type_check
    @classmethod
    def load_from_dict(
        cls,
        data: dict,
        possible_bb_set_ids: frozenset[str],
        static_comp_lookup: Optional[dict[str, str | None]] = None,
    ) -> "ReactionTree":
        """
        Load a reaction tree from reaction steps defined in dict format

        Capable of loading a reaction tree from a library definition json
        by passing the list of reaction steps defined in the in "reactions" key.

        Parameters
        ----------
        data: dict
            the dictionary containing reaction step data
        possible_bb_set_ids: frozenset[str]
            the list of possible building block set ids
            needed to determine what is a static reactant
        static_comp_lookup: Optional[dict[str, str | None]], default=None
            optional lookup table for static components
            if provided, reactants that are static components
            will be looked up in this table to get their SMILES *before*
            being checked if they are a valid SMILES string

        Returns
        -------
        ReactionTree
        """
        if static_comp_lookup is None:
            _static_comp_lookup = dict()
        else:
            _static_comp_lookup = static_comp_lookup

        _current_reaction_ids: set[str] = set()
        reaction_steps: list[ReactionStep] = list()
        for rxn_step_name, rxn_step_data in data.items():
            # parse rxn step id
            rxn_step_id = rxn_step_data.get("step_id", rxn_step_name)
            if rxn_step_id in _current_reaction_ids:
                raise ReactionParsingError(
                    f"duplicate reaction step id '{rxn_step_id}' found in reaction definition; "
                    f"make sure all reaction `step_ids` are unique"
                )
            else:
                _current_reaction_ids.add(rxn_step_id)

            # parse reactions
            try:
                rxn_smarts = rxn_step_data["rxn_smarts"]
            except KeyError as e:
                raise ReactionParsingError(
                    f"reaction step {rxn_step_name} missing required 'rxn_smarts' field"
                ) from e

            if isinstance(rxn_smarts, list):
                reactions = [Reaction(smart) for smart in rxn_smarts]
            else:
                reactions = [Reaction(rxn_smarts)]

            # parse reactants
            reactants: list[Reactant] = list()
            try:
                reactant_ids = rxn_step_data["reactants"]
            except KeyError as e:
                raise ReactionParsingError(
                    f"reaction step {rxn_step_name} missing required 'reactants' field"
                ) from e
            if not isinstance(reactant_ids, list):
                raise ReactionParsingError(
                    f"reaction step {rxn_step_name} 'reactants' field "
                    f"must be a list of reactant ids; got type '{type(reactant_ids)}'"
                )
            for reactant_id in reactant_ids:
                if isinstance(reactant_id, list):  # pool reactants:
                    _reactant_pool: list[Reactant] = list()
                    for pooled_reactant_id in reactant_id:
                        if not isinstance(pooled_reactant_id, str):
                            raise ReactionParsingError(
                                f"reaction step '{rxn_step_name}' contains a pooled "
                                f"reactant with an invalid tpye;"
                                f"reactant pools can only have strings; found type "
                                f"'{type(pooled_reactant_id)}' in "
                                f"pool '{reactant_id}'"
                            )
                        _reactant_pool.append(
                            _build_reactant(
                                pooled_reactant_id, possible_bb_set_ids, _static_comp_lookup
                            )
                        )
                    reactants.append(PooledReactant(_reactant_pool))
                elif isinstance(reactant_id, str):
                    reactants.append(
                        _build_reactant(reactant_id, possible_bb_set_ids, _static_comp_lookup)
                    )
                else:
                    raise ReactionParsingError(
                        f"reaction step '{rxn_step_name}' has an invalid "
                        f"reactant id of type '{type(reactant_id)}'; "
                        f"reactant ids must be strings or lists of strings "
                        f"(for pooled reactants)"
                    )

            reaction_steps.append(
                ReactionStep(
                    step_name=rxn_step_name,
                    step_id=rxn_step_id,
                    reaction=reactions if len(reactions) > 1 else reactions[0],
                    reactants=reactants,
                )
            )
        return cls(reaction_steps)

    def get_thread_for_bb_subset_ids(self, bb_subset_ids: frozenset[str]) -> ReactionThread:
        """
        Return the reaction thread that matches the given building block set ids

        Parameters
        ----------
        bb_subset_ids: frozenset[str]
            the set of building block set ids to match

        Returns
        -------
        ReactionThread
            the matching reaction thread

        Raises
        ------
        ReactionError
            if no matching reaction thread is found
        """
        for thread in self.threads:
            if frozenset(thread.bb_set_reactant_ids) == bb_subset_ids:
                return thread
        raise ReactionError(
            f"no reaction thread found matching building block set ids {bb_subset_ids}"
        )


def _build_reactant(
    reactant_id: str, bb_ids: frozenset[str], static_comp_map: dict[str, str | None]
) -> Reactant:
    """
    Build the correct reactant from its id string

    Parameters
    ----------
    reactant_id: str
        the reactant id string
    bb_ids: frozenset[str]
        the list of possible building block set ids
    static_comp_map: dict[str, str | None]
        lookup table for static components

    Returns
    -------
    Reactant
    """
    if reactant_id.startswith("product_"):
        return ProductReactant(step_id=reactant_id[8:])
    else:
        try:
            bb_set_id, _ = parse_building_block_subset_id(reactant_id)
        except ValueError:
            bb_set_id = reactant_id
        if bb_set_id in bb_ids:
            return BBSetReactant(reactant_id)
        elif reactant_id in static_comp_map:
            try:
                _value = static_comp_map[reactant_id]  # for mypy type checking
                if _value is not None:
                    return StaticReactant(to_mol(_value, fail_on_error=True))
            except ValueError as e:
                raise ReactionParsingError(
                    f"static reactant '{reactant_id}' has an invalid SMILES"
                ) from e
            raise ReactionParsingError(
                f"static reactant '{reactant_id}' is missing, but was used in a reaction step"
            )
        else:
            try:
                return StaticReactant(to_mol(reactant_id, fail_on_error=True))
            except ValueError as e:
                raise ReactionParsingError(
                    f"reactant id '{reactant_id}' appears to be a static reactant, "
                    f"but is not a valid SMILES"
                ) from e
