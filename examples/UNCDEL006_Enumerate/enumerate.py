#!/usr/bin/env python3
"""Example script to enumerate DEL006 using DELi."""

from pathlib import Path

from deli.configure import set_deli_data_dir
from deli.dels.combinatorial import CombinatorialLibrary

data_dir = Path(__file__).resolve().parent.parent / "example_deli_data_dir"
set_deli_data_dir(data_dir)

lib = CombinatorialLibrary.load("DEL006")

# Enumerate a single compound by building block IDs
result = lib.enumerate_by_bb_ids(["A010", "B047", "C060"])
print(f"SMILES: {result.smi}")
print(f"Building blocks: {[bb.bb_id for bb in result.building_blocks]}")

# To enumerate the full library, uncomment:
# lib.enumerate_to_file("DEL006_enumerated.csv", use_tqdm=True)
