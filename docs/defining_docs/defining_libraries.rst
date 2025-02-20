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
-------------------
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

``id``
^^^^^^
This element should hold a string.
This is the ID associated with the library.
It should be the same as the name of the file

``library_tag``
^^^^^^^^^^^^^^^
This element should be a string.
This should be the full DNA sequence associated with
library

.. note::
    if your DNA barcode has a overhang region for the
    library tag it **should not** be included in the
    ``library_tag``

``bb_sets``
^^^^^^^^^^^
This element should be a list of strings.
It should include the path to each Building Block Set
that is used in the library. They should be listed in
the same order as they are synthesized in, e.i. cycle A/1,
then cycle B/2 and so on.

Instead of using file paths, you can also use the name of
building block set file if you have configured the :ref:`deli-data-dir-ref`

``barcode_schema``
^^^^^^^^^^^^^^^^^^
This element should be a string.
It should be the path the Barcode Schema file this library follows.
You can also use the name of the barcode schema if you have
configured the DELi Data Directory

.. _scaffold-ref:

``scaffold``
^^^^^^^^^^^^
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

``reactions``
^^^^^^^^^^^^^
All compounds in a library should follow the same reaction scheme to
create the compounds.
This is usually carried out in steps.
For example: Cycle1 is attached to DNA, Cycle2 is attached to Cycle1,
Cycle3 is attached to the Cycle1-2 product and so on.

While not required, if the reaction schema is provided, it will enable
DELi to enumerate out the full SMILES for any compound in the library.

.. note::
    Enumeration is expensive, so it is recommend to run a full library
    enumeration once and the results saved in a database mapping DEL IDs
    to the full SMILES. See Enumerating Libraries for more info

The reaction is defined as a list of dictionaries,
with each dictionary defining a single 'step' of the reaction scheme.
Each step should include the follow key-value pairs:

* ``step``: the order/position of this this step, starting from 1. Each reaction steps should have a unique step position and they should be sequential (i.e. 1, 2, 3 is valid but 1, 3, 4 is not as 2 is missing)

* ``rxn_smarts``: the SMARTS that defines the reaction that will occur.
  for more info on how reaction SMARTS are defined, see `the Daylight docs <https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html>`_

* ``reactants``: this value should be a list, and contain the BB cycle ids or
  SMILES of the reactants. Only use the SMILES if all compounds in the library
  are reacting with the same compound. You can also include 'scaffold' if the
  reaction includes the DEL scaffold defined in :ref:`scaffold-ref`. You can reference the product of any other steps by
  using ``product_<step>``, i.e. ``product_1`` is the product from the first
  reaction step.

  .. warning::
        Reactants in the list **MUST** match the order they are used in the
        reactants part of the reaction SMARTS. This is because RDKit expects this
        For example, for an amide reaction [NH2:1].[C:2](=[O:3])(O)>>[C:2]
        (=[O:3])[NH:1], you need to order the recants as ['amine', 'carbo-acid'.
        If you ordered it as ['carbo-acid', 'amine'] the reaction
        would not be carried out.


An example of a reaction step dictionary for an amide coupling between
a BB set 'BB1' and the scaffold would be
::
    {
        "step": 1
        "rxn_smarts": "[NH2:1].[C:2](=[O:3])(O)>>[C:2](=[O:3])[NH:1]"
        "reactants": ["BB1", "scaffold"]
    }

Reactions are not limited to 2 reactants, it can be any number that matches
the reaction SMARTS.
An example of a reaction step dictionary for a three step reaction is
::
    {
        "step": 1
        "rxn_smarts": "[NH2].[C(=O)O].[OH]>>[C(=O)N].[C(=O)O]"
        "reactants": ["BB1", "BB2", "c1ccccc1[OH]"]
    }

.. note::
    reaction steps do not need to be ordered in the list as long as the ``step``
    key still provides the correct order. However, for readability providing the
    steps in the order they occur is preferred

``dna_barcode_on``
^^^^^^^^^^^^^^^^^^

This element should be a string.
It should match the BB Set ID of one of the ``bb_sets``.
This tells DELi which building block the DNA is attached to.
DELi uses this during some data analysis methods.


Saving in the DELi Data Directory
=================================
if you have configured the DELi Data Directory,
you should save the library files in the
``libraries`` sub-directory.
The name of the file should be the same as the
``id`` element in the json file
