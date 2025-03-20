==========
Enumerator
==========

Overview
========

The DELi enumerator is a specialized component for enumerating DNA-encoded library (DEL) compounds. 
It generates compounds by combining building blocks according to specified reaction workflows, producing SMILES representations of the resulting molecules.

How it Works
===========

The enumerator performs the following:

1. Takes a reaction workflow definition
2. Takes sets of building blocks
3. Takes a scaffold SMILES string (optional)
4. Combines building blocks according to the reaction steps
5. Generates SMILES strings and RDKit mol objects for the resulting compounds

Command-Line Interface
=====================

The enumerator provides a CLI for direct compound generation without Python coding:

.. code-block:: bash

    deli enumerate <path/to/library.json> [OPTIONS]

Required Arguments
------------------

* ``<path/to/library.json>`` - JSON file containing reaction workflow, building block definitions, and optional scaffold

Options
-------

* ``-o``, ``--out_path <file.csv>`` - Write results to CSV file (default: "enumerated_library.csv")
* ``--debug`` - Enable debug mode for detailed output
* ``--no-tqdm`` - Disable progress bar display
* ``--help`` - Show help message and exit

JSON Formatting
----------------

The JSON file should contain this structure:

.. code-block:: json

    {
        "bb_sets": ["DEL004_BBA", "DEL004_BBB"],
        "reactions": [
            {
                "step": 1,
                "rxn_smarts": "...",
                "reactants": ["DEL004_BBA", "scaffold"]
            },
            {
                "step": 2,
                "rxn_smarts": "...", 
                "reactants": ["product_1", "DEL004_BBB"]
            }
        ],
        "scaffold": "Nc1ccccc1"
    }

This is the same JSON format as a typical library JSON. For full details see :doc:`defining_libraries`.
The items in "bb_sets" must be loadable from the specified building block files.
Similarly, the items in "reactions" must contain all necessary information for reaction definition.

Example Commands
----------------

Generate compounds and save to CSV:

.. code-block:: bash

    deli enumerate library.json -o my_compounds.csv

Python Usage
============

To use Python instead of the command-line interface, you need to provide:

1) A reaction workflow
2) A list of building block sets
3) An optional scaffold SMILES string

.. code-block:: python

    from deli.dels.enumerator import DELEnumerator
    from deli.dels.reaction import ReactionWorkflow
    from deli.dels.building_block import BuildingBlock, BuildingBlockSet

    # Define the reaction workflow
    reaction_workflow = ReactionWorkflow.load_from_json_list(
        rxn_list = [
            {"step": 1, "rxn_smarts": "...", "reactants": ["DEL004_BBA", "scaffold"]},
            {"step": 2, "rxn_smarts": "...", "reactants": ["product_1", "DEL004_BBB"]},
            {"step": 3, "rxn_smarts": "...", "reaction": ["product_2", "DEL004_BBC"]}
        ],
        bb_set_ids: ["DEL004_BBA", "DEL004_BBB", "DEL004_BBC"]
    )

    # Define building block sets
    bb_sets = [
        BuildingBlockSet("DEL004_BBA", [BuildingBlock("BB_0", "ATGCTGTA", "CC(=O)Nc1ccc(O)cc1"), BuildingBlock("BB_1", "ATGCAGTA", "CC(=O)Nc1ccccc1")]),
        BuildingBlockSet("DEL004_BBB", [BuildingBlock("BB_2", "CTGCTGTA", "CCNc1ccc(O)cc1"), BuildingBlock("BB_3", "TCAGCAGTA", "CCNc1ccccc1")]),
        BuildingBlockSet("DEL004_BBC", [BuildingBlock("BB_4", "TCGCTGTA", "CCNc1cccc(O)cc1"), BuildingBlock("BB_5", "TTAGCAGTA", "Nc1cccc(O)cc1")])
    ]

    # Optional scaffold
    scaffold = "Nc1ccccc1"

    # Initialize the enumerator
    enumerator = DELEnumerator(reaction_workflow, bb_sets, scaffold)

Initialization from JSON
-----------------------

Like with the command-line interface, the DELEnumerator can also be initialized from the same JSON file:

.. code-block:: python

    import json
    from deli.dels.enumerator import DELEnumerator
    
    enumerator = DELEnumerator.load("library.json")

Enumeration
-------------

The enumerator can generate all possible compounds.
Or, it can also generate a specific compound based on building block IDs.

.. code-block:: python

    # Enumerate all possible compounds
    for compound in enumerator.enumerate():
        print(compound.smi)  # Access SMILES of enumerated compound
        print(compound.mol)  # Access RDKit mol object

    # Enumerate a specific compound
    building_block_id_map = {
        "BB1": "10", # Building block with ID "10" in the BB1 set
        "BB2": "12", # Building block with ID "12" in the BB2 set
        "BB3": "18" # Building block with ID "18" in the BB3 set
    }
    compound = enumerator.get_enumerated_compound_from_bb_ids(building_block_id_map)

Writing Results
-------------

You can write enumeration results to a pandas dataframe or a CSV file:

The pandas dataframe has headers: "smi", "mol", multiples of "[BuildingBlockSet ID]".
WARNING: This could be extremely memory hungry for large libraires.

The CSV has headers: "SMILES", multiples of "[BuildingBlockSet ID]", and an optional "CompoundID".
The "CompoundID" is a unique identifier for each compound, which can be specified by 
a callable function or left as None for default behavior.

.. code-block:: python

    # Write enumeration results to pandas dataframe
    df = enumerator.enumerate_to_pandas()

    # Write enumeration results to CSV
    def compound_id_function(building_block_id_map: dict[str, str]) -> str:
        unique_id = building_block_id_map["BB1"] + "_" + building_block_id_map["BB2"]
        return unique_id

    enumerator.enumerate_to_csv_file(
        out_path="enumerated_compounds.csv",
        compound_id_function=compound_id_function, # can be set to None
        use_tqdm=True  # Show progress bar
    )

Error Handling
============

The enumerator includes robust error handling:

* Failed enumerations return a FailedEnumeratedDELCompound
* Building block validation occurs during initialization

For failed enumerations, the SMILES string will be set to "ENUMERATION_FAILED".
