#!/usr/bin/env python3
"""Efficient script to enumerate DEL006 using DELi with custom building blocks and reactions"""

import sys
import os
import pandas as pd

from deli.dels.library import Library
from deli.dels.building_block import BuildingBlockSet, BuildingBlock
from deli.dels.reaction import ReactionWorkflow, ReactionStep, Reaction, StaticReactant, BBSetReactant
from rdkit import Chem

# 1. Build building block sets from CSV files using ID and SMILES
def create_building_blocks_from_csv(csv_file, set_id):
    df = pd.read_csv(csv_file)
    building_blocks = []
    for _, row in df.iterrows():
        # Use the actual SMILES from the CSV file
        bb = BuildingBlock(str(row['ID']), row['SMILES'])
        building_blocks.append(bb)
    return BuildingBlockSet(set_id, building_blocks)

bba_set = create_building_blocks_from_csv('DEL006_A.csv', 'DEL006_A')
bbb_set = create_building_blocks_from_csv('DEL006_B.csv', 'DEL006_B') 
bbc_set = create_building_blocks_from_csv('DEL006_C.csv', 'DEL006_C')

# 2. Define reaction steps using ReactionStep and Reaction (SMARTS, Reactants)
reactions = [
    ReactionStep(1, Reaction("[#6:1]-C(=O)-O.[15N:2]>>[15N:2]-C(=O)-[#6:1]"), [BBSetReactant("DEL006_A"), StaticReactant(Chem.MolFromSmiles("[15N]"))]),
    ReactionStep(2, Reaction("[#8]=[#6](-[#8]-[#6]-[#6]1-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-[#6]2:[#6]-1:[#6]:[#6]:[#6]:[#6]:2)-[#7:1][*:2]>>[#7:1]([H])[*:2]"), [BBSetReactant("product_1")]),
    ReactionStep(3, Reaction("[Nh:2].[#6:3](=O)-O >>[#6:3](=O)-[N:2]"), [BBSetReactant("product_2"), StaticReactant(Chem.MolFromSmiles("OC(C1=CC=C(F)C(N(=O)=O)=C1)=O"))]),
    ReactionStep(4, Reaction("[#6:0]1(-[*:1]):[#6]:[#6]:[#6:2](-[#9]):[#6:4](-[#7+](-[#8-])=[#8]):[#6]:1.[NH2:6][*:3]>>[#6:0]1(-[*:1]):[#6]:[#6]:[#6:2](-[NH:6][*:3]):[#6:4](-[#7+](-[#8-])=[#8]):[#6]:1"), [BBSetReactant("product_3"), BBSetReactant("DEL006_B")]),
    ReactionStep(5, Reaction("[#6:0]1(-[*:1]):[#6]:[#6]:[#6:2](-[#7][*:3]):[#6:4](-[#7+](-[#8-])=[#8]):[#6]:1.[CX3H1:6](=O)[*:5]>>[#6:0]1(-[*:1]):[#6]:[#6]:[#6:2]2-[#7]([*:3])-[#6:6](-[*:5])=[#7]-[#6:4]:2:[#6]:1"), [BBSetReactant("product_4"), BBSetReactant("DEL006_C")])
]

# 3. Create reaction workflow
workflow = ReactionWorkflow(reactions)

# 4. Create library with custom building blocks and reactions
lib = Library('DEL006', [bba_set, bbb_set, bbc_set], workflow)

# 5. You can enumerate the entire library. This might take a while depending on the size of the library
# lib.enumerate_to_file('DEL006_enumerated.csv', use_tqdm=True)

# 6. Or you can enumerate specific combinations
result = lib.enumerate_by_bb_ids(['A10', 'B47', 'C60'])
print(f"SMILES: {result._smiles}")
print(f"Building blocks: {[bb.bb_id for bb in result.building_blocks]}")