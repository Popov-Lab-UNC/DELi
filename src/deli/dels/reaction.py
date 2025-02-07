"""define reaction classes"""

import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions
from typing import Dict, List, Union, Optional

def clean_smi(smi):
    comps = smi.split('.')
    smi = comps[np.argmax([len(_) for _ in comps])]
    return smi

class Reaction:
    """Class for reaction definition"""
    def __init__(self, rxn_smarts: str):
        """
        Initialize the reaction. Describes a singular chemical reaction.

        Parameters
        ----------
        rxn_smarts: str
            SMARTS or SMIRKS to define a reaction
        """
        self.rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        self.num_reactants = self.rxn.GetNumReactantTemplates()
        
    def run(self, reactants: Union[List[Union[str, Chem.Mol]], Union[str, Chem.Mol]]):
        """
        Run reaction on input molecules/SMILES
        
        Parameters
        ----------
        reactants: list of SMILES/mol or SMILES/mol if only one reactant

        Returns
        -------
        Products
        
        """
        # Handle single reactant case
        if self.num_reactants == 1 and not isinstance(reactants, list):
            reactants = [reactants]
            
        # Validate input count
        if len(reactants) != self.num_reactants:
            raise ValueError(f"Expected {self.num_reactants} reactants, got {len(reactants)}")
            
        # Convert SMILES to Mol objects and validate
        processed = []
        for i, mol in enumerate(reactants):
            if isinstance(mol, str):
                mol = Chem.MolFromSmiles(mol)
                if not mol:
                    raise ValueError(f"Invalid SMILES at position {i}")
            if not isinstance(mol, Chem.Mol):
                raise TypeError(f"Expected SMILES string or Mol object at position {i}")
            processed.append(mol)
            
        # Run reaction
        products = self.rxn.RunReactants(processed)
        if not products:
            return None
            
        # Check each reaction outcome produces exactly one product
        for outcome in products:
            if len(outcome) != 1:
                raise ValueError("Reaction produced an outcome with more than one product.")
        
        # Take the first outcome's product
        product = products[0][0]
        
        # Sanitize the product
        try:
            Chem.SanitizeMol(product)
        except:
            return None
        
        return product

class MultiReactionPipeline:
    """Class for managing and executing a sequence of chemical reactions."""
    def __init__(self):
        self.steps: Dict[int, tuple] = {}  # {step: (Reaction, reactants_list)}

    def __len__(self) -> int:
        """Return the number of reaction steps in the pipeline"""
        return len(self.steps)

    @classmethod
    def from_json(cls, json_path: str) -> "MultiReactionPipeline":
        """Initialize pipeline from JSON file"""
        with open(json_path, "r") as f:
            data = json.load(f)
        return cls.from_dict(data)

    @classmethod
    def from_dict(cls, data: dict) -> "MultiReactionPipeline":
        """Initialize pipeline from dictionary (matches JSON structure)"""
        pipeline = cls()
        
        # Validate and sort reactions
        reactions = sorted(data["reactions"], key=lambda x: x["step"])
        seen_steps = set()
        
        for rxn in reactions:
            step = rxn["step"]
            if step in seen_steps:
                raise ValueError(f"Duplicate step number {step} in reaction definition")
            seen_steps.add(step)
            
            pipeline.add_reaction(
                step=step,
                rxn_smarts=rxn["rxn_smarts"],
                reactants=rxn["reactants"]
            )
        
        return pipeline
        
    def add_reaction(
        self,
        step: int,
        rxn_smarts: str,
        reactants: Union[List[str], str]
    ) -> None:
        """
        Add a reaction step to the pipeline.

        Parameters
        ----------
        step: int
            The step number (execution order).
        rxn_smarts: str
            SMARTS/SMIRKS defining the reaction.
        reactants: list of str or str
            Reactant identifiers (e.g., ["product_0", "reactant_1"]).
        """
        # Ensure reactants is a list
        if not isinstance(reactants, list):
            reactants = [reactants]
        
        # Create Reaction instance and validate reactant count
        reaction = Reaction(rxn_smarts)
        if reaction.num_reactants != len(reactants):
            raise ValueError(
                f"Step {step}: Reaction expects {reaction.num_reactants} reactants, "
                f"but {len(reactants)} were provided."
            )
        
        self.steps[step] = (reaction, reactants)
        
    def run(
        self,
        reactants_dict: Dict[str, Union[str, Chem.Mol]]
    ) -> Dict[int, Optional[Chem.Mol]]:
        """
        Execute the multi-step reaction pipeline.

        Parameters
        ----------
        reactants_dict: dict
            Maps reactant identifiers (e.g., "reactant_1") to SMILES/Chem.Mol.

        Returns
        -------
        dict
            Maps step numbers to their products. Returns `None` if any step fails.
        """
        # Sort steps by step number
        sorted_steps = sorted(self.steps.items(), key=lambda x: x[0])
        products: Dict[int, Optional[Chem.Mol]] = {}
        
        for step, (reaction, reactants) in sorted_steps:
            resolved_reactants = []
            for r in reactants:
                # Resolve reactant references (e.g., "product_0" or "reactant_1")
                if r.startswith("product_"):
                    # Get product from previous step
                    product_step = int(r.split("_")[1])
                    if product_step not in products:
                        raise ValueError(f"Step {step}: Missing product_{product_step}.")
                    mol = products[product_step]
                    if mol is None:
                        return products  # Pipeline failed earlier
                else:
                    # Get reactant from input dictionary
                    if r not in reactants_dict:
                        raise ValueError(f"Step {step}: Reactant '{r}' not provided.")
                    mol = reactants_dict[r]
                
                resolved_reactants.append(mol)
            
            # Run the reaction
            product = reaction.run(resolved_reactants)
            products[step] = product
            
            # Stop if any step fails
            if product is None:
                break
        
        return products