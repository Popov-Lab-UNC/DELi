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
        "doped": [
            {
                "compound_id": "CONTROL_001",
                "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
                "bb1": "AGCTAGGTC",
                "bb2": "TTCGAACTG",
                "bb3": "GGCATCGTA"
            },
            {
                "compound_id": "CONTROL_002",
                "smiles": "C1=CC=C(C=C1)C=O",
                "bb1": "GATATTTCC",
                "bb2": "AAATTTGGG",
                "bb3": "CCGATGATG"
            }
        ],
        "reactions": {
            "step1": {
                "step_id": 1,
                "rxn_smarts": "[#6:1]-C(=O)-O.[15N:2]>>[15N:2]-C(=O)-[#6:1]",
                "reactants": ["DEL006_BBA", "[15N]"]
            },
            "step2": {
                "step_id": 2,
                "rxn_smarts": "[#8]=[#6](-[#8]-[#6]-[#6]1-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-[#6]2:[#6]-1:[#6]:[#6]:[#6]:[#6]:2)-[#7:1][*:2]>>[#7:1]([H])[*:2]",
                "reactants": ["product_1"]
            },
            "step3":
            {
                "step_id": 3,
                "rxn_smarts": "[Nh:2].[#6:3](=O)-O >>[#6:3](=O)-[N:2]",
                "reactants": ["product_2", "OC(C1=CC=C(F)C(N(=O)=O)=C1)=O"]
            },
            "step4":
            {
                "step_id": 4,
                "rxn_smarts": "[#6:0]1(-[*:1]):[#6]:[#6]:[#6:2](-[#9]):[#6:4](-[#7+](-[#8-])=[#8]):[#6]:1.[NH2:6][*:3]>>[#6:0]1(-[*:1]):[#6]:[#6]:[#6:2](-[NH:6][*:3]):[#6:4](-[#7+](-[#8-])=[#8]):[#6]:1",
                "reactants": ["product_3", "DEL006_BBB"]
            },
            "step5":
            {
                "step_id": 5,
                "rxn_smarts": "[#6:0]1(-[*:1]):[#6]:[#6]:[#6:2](-[#7][*:3]):[#6:4](-[#7+](-[#8-])=[#8]):[#6]:1.[CX3H1:6](=O)[*:5]>>[#6:0]1(-[*:1]):[#6]:[#6]:[#6:2]2-[#7]([*:3])-[#6:6](-[*:5])=[#7]-[#6:4]:2:[#6]:1",
                "reactants": ["product_4", "DEL006_BBC"]
            },
        "dna_barcode_on": "DEL006_BBA",
        "scaffold": "C1CCCCC1"
    }

Notice that there is no "ID" or "Name" element. That is because DELi assumes/asserts that
the base name of the file is the ID of the library. Thus, if this file was named
``path\to\file\DEL006.json``, the ID of the library would be ``DEL006``.

These are the named elements required by the library file:

- ``bb_sets``

These are the named elements that are optional but recognized by DELi:

- ``barcode_schema``
- ``doped``
- ``reactions``
- ``dna_barcode_on``
- ``scaffold``

You can include any other sections you like, but DELi will ignore them.
Below is a deeper dive into each of these sections.

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
maps a barcode section "name" to a dictionary of the section's info.

.. note::
    Variable DNA in DELi is represented by the letter 'N'; DELi only accepts 5 type of nucleotides:
    'A', 'G', 'C', 'T' and 'N'. If any other letters are used in the tag sequences, DELi will raise an exception.

The DNA info dictionary can have the following elements:
- ``tag``: the DNA sequence of the section. This is a string of nucleic acids (A, T, C, G) or
  'N's if the section is variable (with the number of N's equal to the expected length)
- ``overhang``: the overhang sequence of the section. This section is optional, and only
  should be included if the section actually has an overhang. Overhangs are the regions used to promote
  DNA ligation and must be static, even in a variable section. That means it should be a string
  of nucleic acids (A, T, C, G). This is most useful in cases where variable sections, like a building block section
  all share the same overhang sequence. Rather than adding it to each building block tag, you can specify it here
  and only list the variable DNA in the building block files. This can improve the quality of decoding, since error
  correction will not be applied to the overhang region (which is wasteful since it is static).
- ``error_correction``: the error correction mode used for this section. This is an optional region and should only
  be included for a variable section that was built to be error corrected.
  This is often done for the building block sections.
  DELi support several modes of error correction, for details on how to specify them for sections,
  see the :ref:`error correction docs <error-correction-docs>`.

There are three types of barcode sections: variable, static, and mixed. Variable sections mean
we expect that part of the barcode to vary *completely* read by read. This means the ``tag``
element is all 'N' (for an unknown nucleotide). Static sections means all barcodes in the DEL
will have the same sequence for this section (no unknown nucleotides). Mixed sections are
a mixture of both, they can have static regions mixed with variable regions. For example,
a mixed section could have a tag like "AGTNNNCGA", where the first three and last three nucleotides
are static, but the middle three are variable.

DELi will automatically determine the type of section based on the ``tag`` element.
Variable only have 'N's, static only have A, T, C, G, and mixed have a combination of both.

Some barcode sections are reserved by DELi. They are used to handle core processes for decoding.
They also get loaded as unique classes in DELi to further separate them from unreserved sections.
Unlike unreserved sections, reserved sections must follow specific naming conventions and rules.

There are two required reserved sections in the barcode schema for a DEL:

- ``library``: this section name is for the part of the DNA barcode used to decode the
  library the compound originates from. It is **required** and must be static.
- ``bb##``: sections that are for building block cycles are BB<cycle number>. The number of building block
  sections should match the number of building block sets specified in
  the :ref:`bb_sets section <bb-set-sec>`. They are mapped by the building block cycle number, so
  the building block set for cycle 1 is represented by the barcode section name BB1.
  If there is a mismatch (e.g. a cycle without a tag section, or tag without a BB set)
  DELi will raise an exception. DELi requires these sections to be variable, and there must be
  at least two of them (since there must be at least two building block sets in any library).

There are also some optional reserved sections:

- ``umi``: this section is for the unique molecular identifier (UMI). It must be variable, though it
  it is optional. It should only be included if the library was designed to have a UMI in each barcode.
  It this section is include, DELi will automatically collect umi corrected counts, otherwise it will
  only provide raw count during decoding.

.. Note::
    While the UMI section is optional to DELi, it really isn't optional in DEL selection (if you want
    your data to mean anything). DELi will warn you by default if your DEL lacks one.

- ``primer_*``: any section that starts with "primer\_" will be loaded
  as a primer section. Primer sections must be static. If your real primer section has variable
  or random bases in it, you should name it something else (for example, "random_primer") to avoid DELi
  trying to load it as a primer section.

.. _doped-compounds-sec::

``doped``
^^^^^^^^^
This section is optional and used to define any :ref:`doped compounds <doped-compounds>` into your library.
If you are not doping any compounds into your library, this section should not be provided.
This element is a list of dictionaries, with each dictionary defining a single doped compound.

Each dictionary should include the following key-value pairs:
- ``compound_id``: the unique identifier for the doped compound
- ``bb1``, ``bb2``, ``bb3``, ... : the building block tags for each cycle that
  correspond to this doped compound. The number of these elements should match
  the number of building block sets defined in the ``bb_sets`` section.
- ``smiles``: *Optional* the SMILES string for the doped compound

DELi will assert that the tags given in bb1, bb2, etc are unique from the building block tags
used in the building block set for the given cycle. This is to avoid miscalls during decoding.
It will also assert that the compound_id is unique from other doped compounds.

.. warning::
    DELi **will not** check that doped compound ids are unique from any given compound
    id in the library. DEL compound ids are generated on the fly, and running an exhaustive check is
    a waste of time. Make sure your ids will not conflict. You this is easy to do if you make sure the
    DEL library id does not show up in the doped compound id.

.. _reaction-sec-ref:

``reactions``
^^^^^^^^^^^^^
This section describes the reaction scheme used to synthesize
the compounds from a given library.
While not required, if the reaction schema is provided, it will enable
DELi to :ref:`enumerate <enumeration-docs>` the full SMILES for any compound in the library.

.. note::
    Enumeration is expensive, so it is recommend to run a full library
    enumeration once and the results saved in a database mapping DEL IDs
    to the full SMILES. See Enumerating Libraries for more info

The reaction scheme is defined as a dictionary of named of "steps" that
contain information about that reaction step, namely the reaction SMARTS
and the reactants used in that step. DELi uses the names of the
reactions steps and the reactants to build possible reaction paths.
This is discussed in more detail below.
Not all compounds in a given DEL will follow the same reaction path.
For example, some DELs might leverage two different reaction routes to
attach cycle 1 and cycle 2 building blocks. DELi is capable of handling
this scenario, but the user must provide detailed reaction information
regarding which reaction steps are used for which building blocks

The name for each reaction step in the reaction dictionary
should be unique with respect to the library. Each step can
contain the following key-value pairs:

**Required keys:**

* ``rxn_smarts``: the SMARTS that defines the reaction that will occur.
  for more info on how reaction SMARTS are defined,
  see `the Daylight docs <https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html>`_.
  You can also use the name of a predefined reaction in your :ref:`DELi data directory <deli-data-dir-ref>`
  in the reaction subfolder. DELi will auto detect this for you.

  This element can also be a list, rather than a single string. In this case,
  DELi will attempt to run each reaction SMARTS in order until one works. This is useful
  if you have multiple reaction SMARTS that could be used for a given reaction step that
  have a defined priority order. For example, reacting with primary then secondary amines.

  .. note::
      In cases where there are two or more entirely different reactions (that are not related)
      that could be taken, you should instead define a separate reaction step for each. You
      can read more about this below.

* ``reactants``: this value should be a list, and contain strings that reference the follow:
  * A SMILES string for a chemical reactant
  * A building block set id (should match the one of the ones specified in :ref:`the bb set section <bb-set-sec>`
  * A product from a previous reaction step (see below)

.. _reaction-bb-subset-docs:

  In the case that a building block set is used as a reactant, you can further specify if a specific
  subset is used by appending ":::" and a subset name to the building block set name. For example,
  if building block set "BB1" is used, but only the subset "subsetA" is reacting, you would
  specify the reactant as "BB1:::subsetA". See the building :ref:`block set docs <defining_building_blocks>` for more info on
  defining subsets.

  To specify a product from a different reaction step, you can use the syntax "product_<step_name>",
  where <step_name> is the name of the reaction step that produced the product you want to use.
  You can also use the ``step_id`` to keep things cleaner if your reaction step names are long.

  .. note::
      This means that reaction step names cannot match the pattern "product_*"
      (where * is anything). DELi will complain if you try to use such a name.

  You can also use a list of strings to specify the pooling of multiple reactant source.
  For example, if I wanted to react the products of step 1 and step 2 with with building block set "BB1",
  I could specify the reactants as ``reactants: [["product_1", "product_2"], "BB1"]``.
  This would tell DELi to pool the products of step 1 and step 2 together as a single reactant source, so
  in this case there are only two reactants. This is very different
  from ``reactants: ["product_1", "product_2", "BB1"]``, which would tell DELi you are doing a reaction
  that takes in three reactants, the products of step 1, step 2, and building block set "BB1".
  The reactants can be any number of items, as long as they match the number reactant in the reaction SMARTS for
  the given step.

  .. warning::
      DELi uses RDKit to perform the reactions, which requires that the order
      of reactants matches the order in the reaction SMARTS. Be sure that the order
      of the reactants you provide matches the order in the SMARTS.

  Lastly, there are two reserved reactant keywords beyond the three listed above. These are "linker" and
  "scaffold". These keywords tell DELi to use the library's linker or scaffold as a reactant
  in the reaction. The "linker" keyword is also used to set the root of the reaction tree. If you
  do not include one, DELi will create a dummy linker for you to use as the root (this will not impact the reaction products)

  .. warning::
      This means "scaffold" and "linker" (in any capitalization) are reserved names. DELi will complain
      if you try to use them as :ref:`building block set names <defining_building_block_set_names>`.

  This section is probably the most complex part of defining reactions in DELi, as it is what DELi uses to make
  the reaction trees, and if you have configured it wrong DELi will not be able to properly enumerate. There are
  some checks that exist to catch some mistakes, but not all. So be sure to double check that you are specifying
  your reacts correctly give your reaction schema.

**Optional keys**

* ``step_id``: an addition unique identifier for the reaction step.
  This is useful if your reaction step names are long and cumbersome as it can be referenced
  in place of the full reaction step name in files

  .. note::
      `step_id` is always used to reference previous products and reactions.
      If you do not set a `step_id`, DELi will assign it the same value as the step name.

* ``description``: a string description of the reaction step. This is only useful for documentation,
  DELi does not use this information.

* ``pick_fragment``: If your reaction SMART defines multiple product fragments,
  (contains a '.' in the product side)
  use this to specify which one to pick as the product.
  The indexing start at 1, so if you have two possible fragments, you can set this to 1 or 2.
  If this is not provided, DELi will pick the first fragment by default and warn you that you
  did not specify which fragment to pick.

* ``ignore_multiple_products``: A boolean (true/false) value that tells DELi to supress warnings
  in the case where multiple possible products are generated after running the reaction.
  By default, DELi will warn you everytime a reaction results in ambiguous products, then
  select the first to proceed with. If you set this value to true, DELi will not warn you
  and will just pick the first product to proceed with.

  .. note::
      Multiple products are not the same as multiple fragments. A reaction SMARTS can define multiple fragments
      (e.g. with a '.' in the product side), but still only produce one product per fragment.
      Multiple products means that the reaction SMARTS can produce multiple different products
      for the same set of reactants. This can happen if the reaction SMARTS is not specific enough
      or the reactants have multiple reactive sites.

* ``fail_on_multiple_products``: Similar to ``ignore_multiple_products``, but instead of suppressing
  multiple products warnings, DELi will raise an exception and fail the enumeration if multiple products are found.
  This is useful if you want to ensure that your reaction SMARTS are specific enough to only produce
  one product per reaction.

An example of a reaction step dictionary for an amide coupling between
a BB set 'BB1' and the scaffold would be
::
    "cycle_1_reaction": {
        "step_id": "1",
        "rxn_smarts": "[NH2:1].[C:2](=[O:3])(O)>>[C:2](=[O:3])[NH:1]"
        "reactants": ["BB1", "scaffold"]
    }

Reactions are not limited to 2 reactants, it can be any number that matches
the reaction SMARTS.
An example of a reaction step dictionary for a three step reaction is
::
    "cycle 1-2 reaction": {
        "step_id": 1
        "rxn_smarts": "[NH2].[C(=O)O].[OH]>>[C(=O)N].[C(=O)O]"
        "reactants": ["BB1", "BB2", "c1ccccc1[OH]"]
    }

.. note::
    Since DELi will load steps first then build a reaction tree, the order of steps
    in the dictionary does not matter.
    DELi will sort them out based on the reactants specified for each step.


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
