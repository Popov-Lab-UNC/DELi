.. _defining_libraries:

==================
Defining Libraries
==================

A DNA encoded library (DEL) is defined as a collection of several parts, some of which are
optional depending on what you want to do with the library:
    * The name/id of the library (REQUIRED)
    * The building block sets used for each cycle (REQUIRED)
    * The DNA tag barcode design used for the library (OPTIONAL)
    * Which building block cycle the DNA barcode is attached to (OPTIONAL)
    * The reactions that linked cycles together (OPTIONAL)
    * A scaffold that all compounds in the library will share (OPTIONAL)

Below is a deeper dive into each of these parts, and outlines
how to define them so that DELi can understand and load your DEL


Library File Format
-------------------
**Libraries are defined following JSON format**

Below is and example of a library json file
::
    {
        "barcode_schema": {
            "library": {
                "tag":  "CCTTGGCACCCGAGAATTCCAATCGCTGA",
                "overhang": "CTA"
            },
            "bb1": {
                "tag": "NNNNNNNN",
                "overhang": "ACG",
                "error_correction": "levenshtein_dist:1,asymmetrical"
            },
            "bb2": {
                "tag": "NNNNNNNN",
                "overhang": "GAT",
                "error_correction": "hamming_dist:1"
            },
            "bb3": {
                "tag": "NNNNNNNN",
                "overhang": "TGC",
                "error_correction": "hamming_code:8_4"
            },
            "pre-umi": {
                "tag": "AATGCCAGTACG"
            },
            "umi": {
                "tag": "NNNNNNNNNNN"
            }
        },
        "bb_sets": [
            {
                "cycle": 1,
                "bb_set_name": "DEL006_BBA"
            },
            {
                "cycle": 2,
                "bb_set_name": "DEL006_BBB"
            },
            {
                "cycle": 3,
                "bb_set_name": "DEL006_BBC"
            }
        ],
        "reactions": [
            {
                "step": 1,
                "rxn_smarts": "[#6:1]-C(=O)-O.[15N:2]>>[15N:2]-C(=O)-[#6:1]",
                "reactants": ["DEL006_BBA", "[15N]"]
            },
            {
                "step": 2,
                "rxn_smarts": "[#8]=[#6](-[#8]-[#6]-[#6]1-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-[#6]2:[#6]-1:[#6]:[#6]:[#6]:[#6]:2)-[#7:1][*:2]>>[#7:1]([H])[*:2]",
                "reactants": ["product_1"]
            },
            {
                "step": 3,
                "rxn_smarts": "[Nh:2].[#6:3](=O)-O >>[#6:3](=O)-[N:2]",
                "reactants": ["product_2", "OC(C1=CC=C(F)C(N(=O)=O)=C1)=O"]
            },
            {
                "step": 4,
                "rxn_smarts": "[#6:0]1(-[*:1]):[#6]:[#6]:[#6:2](-[#9]):[#6:4](-[#7+](-[#8-])=[#8]):[#6]:1.[NH2:6][*:3]>>[#6:0]1(-[*:1]):[#6]:[#6]:[#6:2](-[NH:6][*:3]):[#6:4](-[#7+](-[#8-])=[#8]):[#6]:1",
                "reactants": ["product_3", "DEL006_BBB"]
            },
            {
                "step": 5,
                "rxn_smarts": "[#6:0]1(-[*:1]):[#6]:[#6]:[#6:2](-[#7][*:3]):[#6:4](-[#7+](-[#8-])=[#8]):[#6]:1.[CX3H1:6](=O)[*:5]>>[#6:0]1(-[*:1]):[#6]:[#6]:[#6:2]2-[#7]([*:3])-[#6:6](-[*:5])=[#7]-[#6:4]:2:[#6]:1",
                "reactants": ["product_4", "DEL006_BBC"]
            }
        ],
        "dna_barcode_on": "DEL006_BBA"
        "scaffold": "C1CCCCC1"
    }

Notice that there is no "ID" or "Name" element. That is because DELi assumes/asserts that
the base name of the file is the ID of the library. Thus, if this file was named
``path\to\file\DEL006.json``, the ID of the library would be ``DEL006``.

These are the named elements supported by the library file:
- ``bb_sets``
- ``barcode_schema``
- ``reactions``
- ``dna_barcode_on``
- ``scaffold``

The only section that is always required is the ``bb_sets`` section.
Every other section is optional, though its best practice to add it
if you have the information.

.. _bb-set-sec:

``bb_sets``
^^^^^^^^^^^
This element defines which building blocks sets are used for each
cycle of the library.

It is a list of dictionaries, with each dictionary defining a single
cycle. That dictionary should include the following key-value pairs:
- ``cycle``: the cycle number of the building block set
- ``bb_set_name``: the name of the building block set.
  This can be either the full path to the building block set file
  or just its name if you have configured the :ref:`DELi data directory <deli-data-dir-ref>`

While you can misorder the cycles such that the order of the list
does not following an ascending cycle order (e.g. 1, 3, 2 instead of 1, 2, 3),
DELi will raise an exception if you do this, as it expects the sets to be
listed in order. The cycle number acts more as a check to make sure
this order is correct, because if it is not DELi can behave in unexpected
ways that are hard to debug.

.. _barcode-sec-ref:

``barcode_schema``
^^^^^^^^^^^^^^^^^^
This section is optional and used to define the DNA barcode design used
to tag compounds in the DEL. This is only needed if you are trying to use
DELi to decode; All other functions of DELi can operate with out it.
See :ref:`lib-vs-del` for more details on this.

This element is a dictionary of dictionaries. Each element of the outer dictionary
maps a barcode section "name" to a dictionary of the DNA info of that section.
There are two types of barcode sections: variable and static. Variable sections mean
we expect each compound in the DEL to vary in sequence for that section. Static
means all compounds in the DEL will have the same sequence for this section.
Variable DNA in DELi is represented by the letter 'N' (e.g. NNNNNNNN) and static is the
four nucleic acids: 'AGTC' (e.g. AGTTCGTA).

The DNA info dictionary can have the following elements:
- ``tag``: the DNA sequence of the section. This is a string of nucleic acids (A, T, C, G) or
  'N's if the section is variable (with the number of N's equal to the expect length)
- ``overhang``: the overhang sequence of the section. This section is optional, and only
  should be included if the section actually has an overhang. Overhangs are the regions used to promote
  DNA ligation and must be static, even in a variable section. That means it should be a string
  of nucleic acids (A, T, C, G)
- ``error_correction``: the error correction mode used for this section. This is an optional region and should only
  be included for a variable section that was built to be error corrected.
  This is often done for the building block sections.
  DELi support several modes of error correction, for details on how to specify them for sections,
  see the :ref:`error correction docs <error-correction-docs>`.

When it comes to barcode sections, some section names are reserved, and possibly required, by DELi.
These are:
- library: this section name is for the part of the DNA barcode used to decode the
  library the compound originates from. It is *required* and must be static.
- bb##: sections that are for building block cycles are BB<cycle number>. The number of building block
  sections should match the number of building block sets specified in
  the :ref:`bb_sets section <bb-set-sec>`. They are mapped by the building block cycle number, so
  the building block set for cycle 1 is represented by the barcode section name BB1.
  If there is a mismatch (e.g. a cycle without a tag section, or tag without a BB set)
  DELi will raise an exception. DELi requires these sections to be variable, and there must be
  at least two of them (since there must be at least two building block sets in any library).
- umi: this section is for the unique molecular identifier (UMI). It must be variable, though it
  it is optional. It should only be included if the library was designed to have a UMI in each barcode.
  It this section is include, DELi will automatically collect umi corrected counts, otherwise it will
  only provide raw count during decoding.

.. _reaction-sec-ref:

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
  for more info on how reaction SMARTS are defined,
  see `the Daylight docs <https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html>`_

* ``reactants``: this value should be a list, and contain the BB cycle ids or
  SMILES of the reactants. Only use the SMILES if all compounds in the library
  are reacting with the same compound. You can also include 'scaffold' if the
  reaction includes using the :ref:`scaffold <scaffold-sec-ref>` in the DEL.
  You can reference the product of any other steps by
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
    steps in the order they occur is preferred.


``dna_barcode_on``
^^^^^^^^^^^^^^^^^^

This element should be a string.
It should match the BB Set ID of one of the ``bb_sets``.
This tells DELi which building block the DNA is attached to.
Some DELi methods (and possible future
DELi uses this during some data analysis methods.


.. _scaffold-sec-ref:

``scaffold``
^^^^^^^^^^^^
This element should be a string.
In the case that this library has common central
structure that all compounds in the DEL share it's SMILES can be listed
in the "scaffold" element. If there is not scaffold for the library, this section
should not be provided

.. note::
    Scaffold should not have DNA tags in the Barcode.
    If you're scaffold DOES have a section of the Barcode
    that encodes it, you should treat it like another cycle
    (since it technically is, otherwise it should not be tagged)


.. _lib-vs-del-ref:

``Library`` vs ``DELibrary``
----------------------------
At a high level, DELi defines two "types" of DELs, one called ``Library``
and another called ``DELibrary``. This main distinction here is that a ``Library``
lack information about the DNA tag, where a ``DELibrary`` contains that information.

DELi only needs the DNA tag information if you are trying to do any type of DEL
decoding (converting the raw sequence output to compound counts). Otherwise, this
information is not needed, and is not required to be defined. DELi will let you
know if you are trying to do something that requires the DNA tag information but
your libraries lack the correct information. A notable case for this is if you are
just trying to do analysis on a DEL dataset you were given, with no knowledge of the
tags (or a need to use them)


Saving in the DELi Data Directory
---------------------------------
If you have configured the :ref:`DELi Data Directory <deli-data-dir-ref>`,
you should save the library files in the ``libraries`` sub-directory.
This way, when using DELi, you can reference libraries by their name, rather
than having to know the exact path to the file.
