==================
Defining Libraries
==================
A DNA encoded library is define as a collection of several parts
    * The name of the library
    * The DNA tag associated with the library
    * The building block sets used for each cycle
    * The reactions that linked cycles together
    * An optional base scaffold that all compounds in the library will share
    * A DNA Barcode schema that all barcodes in the library follow
    * Which cycle building block the DNA barcode is attached to

Below is a deeper dive into each of these parts, and outlines
how to define them so that DELi can understand and load your library

Library File Format
==================
**Libraries are defined following JSON format.
Each part of the library has a specific key**

Below is and example of a library schema and
definitions/explanations of each part
::
  {
    "id": "DEL004",
    "library_tag": "GCCAGACC",
    "bb_sets": [
      "DEL004_BBA",
      "DEL004_BBB",
      "DEL004_BBC"
    ],
    "barcode_schema": "MyExampleSchema",
    "scaffold": "NA",
    "reactions": [
      {
          "cycle_id_1": "DEL004_BBA",
          "cycle_id_2": "DEL004_BBB",
          "reaction": "NA"
      },
      {
          "cycle_id_1": "DEL004_BBB",
          "cycle_id_2": "DEL004_BBC",
          "reaction": "NA"
      }
    ],
    "dna_barcode_on": "DEL004_BBA"
  }

- ``id``
This element should hold a string.
This is the ID associated with the library.
It should be the same as the name of the file

- ``library_tag``
This element should be a string.
This should be the full DNA sequence associated with
library

.. note::
    if your DNA barcode has a overhang region for the
    library tag it **should not** be included in the
    ``library_tag``

- ``bb_sets``
This element should be a list of strings.
It should include the path to each Building Block Set
that is used in the library. They should be listed in
the same order as they are synthesized in, e.i. cycle A/1,
then cycle B/2 and so on.

Instead of using file paths, you can also use the name of
building block set file if you have configured the DELi
Data Directory

- ``barcode_schema``
This element should be a string.
It should be the path the Barcode Schema file this library follows.
You can also use the name of the barcode schema if you have
configured the DELi Data Directory

- ``scaffold``
This element should be a string.
In the case that this library has common central
structure that all compounds in the DEL share we can
call that the "scaffold". You can specify the SMILES
of that scaffold here.

If there is not scafold you can leave this section out
OR you can specify 'NA'

.. note::
    Scaffold should not have DNA tags in the Barcode.
    If you're scaffold DOES have a section of the Barcode
    that encodes it, you should treat it like another cycle
    (since it technically is, otherwise it should not be tagged)

- ``reactions``
This element should be a dictionary (see more below).
Every time a new cycle is added, a reaction occurs to link them together.
The reaction that occurs is static between any
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

For a library to be valid there must be at least ``N-1``
reactions where N = num_cycle (plus scaffold if present).
So if I had 3 cycles I need 2 reactions. If I also had a
scaffold I would need 3 reactions

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

- ``dna_barcode_on``
This element should be a string.
It should match the file name of one of the ``bb_sets``.
This tells DELi which building block the DNA is attached to.
DELi uses this during some data analysis methods.

Defining Library Groups
=======================
Most DELs screening involves the screening of a mix
of lirbaries all at once. Often this is a standard mix
that is commonly used. To allow for this, DELi also allows
for the encoding a large set of libraries in one file.
Internally, DELi calls these a ``DELibraryGroup``.

To define a library group in a file, you follow the same definition
for libraries, except you encapsulate them in a list
::
    [
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
        },
        {
            "id": My_Library_2,
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
    ]

Saving in the DELi Data Directory
=================================
if you have configured the DELi Data Directory,
you should save the library files in the
``building_blocks`` sub-directory.
The name of the file should be the same as the
``id`` element in the file
