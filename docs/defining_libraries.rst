==================
Defining Libraries
==================

A DNA encoded library is define as a collection of several parts
    * A DNA Barcode schema
    * Where the Barcode is attached
    * At least two set of building blocks used to synthesize the library
    * The reactions that linked cycles together (not always the same reaction)
    * An optional base scaffold that all compounds in the library will share

Below is a deeper dive into each of these parts, and outlines how to define them so that DELi can
understand and build your library

DNA Barcode Schema
==================
This define the general structure of the DNA barcode that will get attached to every compound in the library.
This is was is later sequences to process the results of the DEL selection experiment.
See Barcode Schema for more details

Building Block Sets
===================
A quick refresher about compounds in DEL:
Each compound is assembled from various building blocks (BB).
The compound is built sequentially: First a build block 1 is
added (often called the cycle1), then building block two is added and a reaction occurs to link them (cycle2)
and so on. Generally, we only see cycles between 2/3 (di/trisynthon libraries).
What allows DELs to have such large size is combinatorics.
If I have 1,000 BBs in cycle1, 1,000 in cycle2 and 1,000 in cycle3 I have 1,000,000,000
(1 billion) unique compounds in the resulting library.

There are some restrictions on what can be in those BB sets.
Since each additional cycle is linked with the same reaction,
all BB in a given set must have the same specific reaction site
to allow a compound from the next cycle to be added.

Thus to define a library we need to have at least 2 BB sets, each set containing
BBs that have the correct reaction sites to allow for the correct reactions to occur.
To define a library you just need to provide the files that contain the information about each set
See Defining Building Block Sets for more information

Scaffold
========
Not all DEL synthesis schemas involve the Building Blocks being linked together.
(e.g. BB1 -- BB2 -- BB3). Some might all share a common initial scaffold that each BB is linked to.
If this is the case, you need to provide this scaffold (as a SMILES).

.. note::
    Scaffold should not have DNA tags in the Barcode. If you're scaffold DOES have a
    section of the Barcode that encodes it, you should treat it like a Building blocks
    (it would be a building block set of size 1)

DNA Barcode Attachment
======================
The DNA Barcode must be linked somehow to the compound, and is generally linked to the first cycle's building block.
However it can be attached to other cycles or to a scaffold. You should specify this when defining the library.
For readability and long term support; you should allows specify it.

.. note::
    if left unspecified, will assume the tag is attached to the "scaffold" if one exist otherwise to BB1

Reactions
=========
Every time a new cycle is added, a reaction occurs to link them together. The reaction that occurs is static between any
two cycles, e.g. cycle1 + cycle2 is always an Amide bond coupling.

Therefore we define a reaction as 3 things:
 * Reactant cycle 1
 * Reactant cycle 2
 * SMARTS/SMIRKS that defines the reactions

We format them as dictionaries
::
    {
        "cycle_id_1": "my_buildblock_set_1"
        "cycle_id_2": "my_buildblock_set_2"
        "reaction": "C(=O)O.OCC>>C(=O)OCC.O"
    }

For a library to be valid there must be at least ``N-1`` reactions where N = num_cycle (plus scaffold if present).
So if I had 3 cycles I need 2 reactions. If I also had a scaffold I would need 3 reactions

.. note::
    DELi supports Macrocyclic DELs (a library where the a the addition of a new cycle will reaction with two previous cycles instead of 1. You define this by simpling also including the required reaction.


You can read more about reaction SMARTS here: https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html

If you have a scaffold that reacts with your building blocks, you can use the keyword "scaffold" instead of the
building block set id/name to specify that reaction
::
    {
        "cycle_id_1": "my_buildblock_set_1"
        "cycle_id_2": "scaffold"
        "reaction": "C(=O)O.OCC>>C(=O)OCC.O"
    }


File Format
===========
To define a library we use a JSON format with several named keys:
 * ``id``; the name/id of the library (i.e. "DEL004"). Should also be the name of the file (``DEL004.json``)
 * ``library_dna_tag``; the DNA tag for this library
 * ``barcode_schema``; the name of the Barcode Schema to use (should have a file with the same name in your ``$BARCODE_SCHEMA_FOLDER`` directory, see Defining Barcode Schemas)
 * ``bb_sets``; A list of bb_set ids to include in the library (should have a file with the same name in your ``$BB_SET_FOLDER`` directory, see Defining Building Block Sets)
 * ``scaffold``; the (optional) scaffold present as SMILES. If not present should be ``null`` for readability (assumes ``null`` if not present)
 * ``reactions``; a list of reactions linking the bb_sets

As an example, I have 5 BB sets, named ``BB_A``, ``BB_B``, ``BB_C``, ``BB_D`` and ``BB_E``, thus I have a folder named
``my_building_blocks`` with
::
    /my_building_blocks
    --- BB_A.json
    --- BB_B.json
    --- BB_C.json
    --- BB_D.json
    --- BB_E.json

I can define a new library with
::
    {
        "id": My_Library_1,
        "bb_sets": [
            "BB_A",
            "BB_C",
            "BB_D"
        ],
        "scaffold": "C1CC(Cl)CCC1"
        "reactions": [
            {
                "cycle_id_1": "BB_A"
                "cycle_id_2": "scaffold"
                "reaction": "[C:1]Cl.[C:2]O>>[C:1]O[C:2]"
            },
            {
                "cycle_id_1": "BB_B"
                "cycle_id_2": "BB_A"
                "reaction": "(C(=O)O).(OCC)>>C(=O)OCC.O"
            },
            {
                "cycle_id_1": "BB_B"
                "cycle_id_2": "BB_C"
                "reaction": "[C:1]Cl.[C:2]O>>[C:1]O[C:2]"
            }
        ],
        "dna_barcode_on": "scaffold"
    }

I can define another one without a scaffold with
::
    {
        "id": My_Library_1,
        "bb_sets": [
            "BB_B",
            "BB_E"
        ],
        "scaffold": null
        "reactions": [
            {
                "cycle_id_1": "BB_B"
                "cycle_id_2": "BB_E"
                "reaction": "(C(=O)O).(OCC)>>C(=O)OCC.O"
            }
        ],
        "dna_barcode_on": "BB_B"
    }
