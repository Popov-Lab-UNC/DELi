"""defines DEL enumerator functions and classes"""

from rdkit import Chem
from rdkit.Chem import AllChem
import json, csv
from itertools import product
from typing import Dict, List, Optional, Iterator, Set, Union, Tuple
from deli.configure import DeliDataLoadable, accept_deli_data_name
from .reaction import MultiReactionPipeline
from .building_block import BuildingBlock, BuildingBlockSet

class EnumeratorBuildError(Exception):
    """error raised when a enumerator build fails"""
    pass

class EnumeratorRunError(Exception):
    """error raised when a enumerator run fails"""
    pass

class DELEnumerator(DeliDataLoadable):
    """Handles DNA-encoded library enumeration using multi-step reaction pipelines."""
    
    def __init__(
        self,
        bb_sets: List[BuildingBlockSet],
        reactions: List[dict],
        scaffold: Optional[str] = None,
    ):
        """
        Initialize enumerator with building blocks and reaction scheme.
        
        Parameters
        ----------
        bb_sets : List[BuildingBlockSet]
            Available building block sets
        reactions : List[dict]
            Reaction definitions with
            - Step number: int
            - Reaction SMARTS: str
            - Reactants: list of IDs (str)
        scaffold : Optional[str]
            Optional scaffold SMILES string
        """
        self.bb_sets = bb_sets
        self.reactions = reactions
        self.scaffold = scaffold
        self.static_reactants = self._find_static_reactants() # holds mol objects
        self.pipeline = self._build_pipeline()

        self._validate_initialization()

    def _find_static_reactants(self) -> Set[str]:
        """Identify static SMILES strings used directly in reactions."""
        static = {}
        if self.scaffold: # includes scaffold
            static['scaffold'] = Chem.MolFromSmiles(self.scaffold)

        for rxn in self.reactions:
            reactants = rxn["reactants"]
            if not isinstance(reactants, list):
                reactants = [reactants]
                
            for r in reactants:
                if (r not in self.bb_sets 
                    and r != "scaffold"
                    and not r.startswith("product_")):
                    static[r] = Chem.MolFromSmiles(r)
        return static

    def _build_pipeline(self) -> MultiReactionPipeline:
        """Construct reaction pipeline from JSON definitions."""
        return MultiReactionPipeline.from_dict({"reactions": self.reactions})
    
    def _validate_initialization(self):
        """Validate critical initialization parameters."""
        if not self.bb_sets:
            raise EnumeratorBuildError("No building block sets provided")
            
        if not self.reactions:
            raise EnumeratorBuildError("No reaction definitions provided")
    
    @classmethod
    @accept_deli_data_name("libraries", "json")
    def load(cls, path: str) -> "DELEnumerator":
        """Load enumerator from DELi data directory."""
        data = json.load(open(path))
        return cls.from_dict(data)

    @classmethod
    def from_dict(cls, data: dict) -> "DELEnumerator":
        """Create enumerator from dictionary data."""
        # Load building block sets
        bb_sets = []
        for bb_name in data["bb_sets"]:
            bb_set = BuildingBlockSet.load(bb_name)
            if not bb_set:
                raise EnumeratorBuildError(f"Building block set {bb_name} not found")
            bb_sets.append(bb_set)

        return cls(
            bb_sets=bb_sets,
            reactions=data["reactions"],
            scaffold=data.get("scaffold")
        )

    def __repr__(self) -> str:
        return (f"<DELEnumerator: {len(self.bb_sets)} BB sets, "
                f"{len(self.reactions)} reactions>")

    def _get_final_product_smiles(self, results: dict) -> Optional[str]:
        """Extract final product SMILES from pipeline results"""
        if not results:
            return None
        final_step = max(results.keys())
        product = results.get(final_step)
        return Chem.MolToSmiles(product) if product else None
    
    def enumerate(self, input_bb_ids: Dict[str, str]) -> Optional[str]:
        """
        Run reactions for specific building block IDs and return final product SMILES
        
        Parameters
        ----------
        input_bb_ids : Dict[str, str]
            Dictionary mapping BB set IDs to building block IDs
            Example: {"BB_A": "A1", "BB_B": "B3"}
            
        Returns
        -------
        Optional[str]
            Final product SMILES if successful, None otherwise
        """
        # Validate input contains all required BB sets
        required_bb_sets = {bb_set.bb_set_id for bb_set in self.bb_sets}
        provided_bb_sets = set(input_bb_ids.keys())
        
        if provided_bb_sets != required_bb_sets:
            return None

        # Retrieve building blocks
        building_blocks = {}
        for bb_set in self.bb_sets:
            bb = bb_set.get_bb_by_id(input_bb_ids[bb_set.bb_set_id])
            if not bb:
                return None
            building_blocks[bb_set.bb_set_id] = bb

        # Create reactants dictionary
        reactants = {bb_set_id: Chem.MolFromSmiles(bb) for bb_set_id, bb in building_blocks.items()}
        reactants.update(self.static_reactants)

        # Run reaction pipeline
        results = self.pipeline.run(reactants)
        if not results:
            return None
            
        # Get final product
        final_step = max(results.keys())
        product_mol = results.get(final_step)
        
        if not product_mol:
            return None
            
        try:
            product_mol.UpdatePropertyCache()
            Chem.SanitizeMol(product_mol)
            return Chem.MolToSmiles(product_mol)
        except Exception:
            return None
    
    def _bb_combination_generator(self) -> Iterator[Dict[str, BuildingBlock]]:
        """Generates {bb_set_id: BuildingBlock} combinations"""
        # Get ordered lists of set IDs and their building blocks
        set_ids = [bb_set.bb_set_id for bb_set in self.bb_sets]
        all_blocks = [bb_set.building_blocks for bb_set in self.bb_sets]
        
        # Generate all combinations across BB sets
        for combo in product(*all_blocks):
            yield {set_id: bb for set_id, bb in zip(set_ids, combo)}
    
    def _build_csv_row(self, 
                      bb_combo: Dict[str, BuildingBlock], 
                      product_smiles: Optional[str]) -> List[str]:
        """Build CSV row with BB metadata and product SMILES"""
        row = []
        # Process BB sets in order of self.bb_sets
        for bb_set in self.bb_sets:
            bb = bb_combo[bb_set.bb_set_id]
            row += [bb.bb_id, bb.smiles]
        row.append(product_smiles or "")
        return row

    def enumerate_all(self, output_file: str) -> None:
        """
        Write all combinations and products to a CSV file
        
        Parameters
        ----------
        output_file : str
            Path to output CSV file with columns:
            [BB_A_ID, BB_A_SMILES, BB_B_ID, BB_B_SMILES, ..., Product_SMILES]
        """
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Write header
            header = []
            for bb_set in self.bb_sets:
                header += [f"{bb_set.bb_set_id}_ID", f"{bb_set.bb_set_id}_SMILES"]
            header.append("Product_SMILES")
            writer.writerow(header)
            
            # Process all combinations
            for bb_combo in self._bb_combination_generator():
                try:
                    reactants = {bb_set_id: Chem.MolFromSmiles(bb.smiles) for bb_set_id, bb in bb_combo.items()}
                    reactants.update(self.static_reactants)
                    
                    # Run reaction
                    results = self.pipeline.run(reactants)
                    product_smiles = self._get_product_smiles(results)
                    
                    # Write row to file
                    row = self._build_csv_row(bb_combo, product_smiles)
                    writer.writerow(row)
                    
                except Exception as e:
                    if self.strict_mode:
                        raise EnumeratorRunError(f"Failed combination {bb_combo}: {str(e)}")
                    # Write partial data for failed reactions
                    row = self._build_csv_row(bb_combo, None)
                    writer.writerow(row)