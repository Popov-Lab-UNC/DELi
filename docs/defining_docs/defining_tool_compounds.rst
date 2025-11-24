.. _defining-tool-compounds::

=======================
Defining Tool Compounds
=======================

:ref:`Tool compounds <tool-compounds>` are defined as JSON files.
They have two main sections: the ``barcode_schema`` section and the ``tool_compounds`` section.
Only the tool_compounds section is required.

Tool compound files can contain multiple tool compounds in a single file.
The only constraint is that all tool compounds in a single file must share the same barcode schema,
and the same library tag (much like all compounds in a given DEL).

Barcode Schema
--------------
Like the :ref:`DEL library definition files <defining_libraries>`, they can
implement a barcode schema. You can :ref:`read the docs <barcode-sec-ref>` on how these are defined.
Barcode schemas are not required if you are just loading a tool compound
not for decoding, but is required if you want to decode tool compounds (just like for
libraries). However, the required reserved sections are a bit different for tool compounds.

There are two required sections in the barcode schema for tool compounds:

- ``library``: This section defines the library tag for the tool compound.
  This is required for all tool compounds, and is used to identify that a read
  is going to map to a tool compound.
  .. note::
        The library tag for a tool compound must be unique from all other
        DEL libraries in the selection. If it is the same, you can try using
        :ref:`doped compounds <doped-compounds-sec>` instead.

- ``compound_tag``: This section defines the compound tag for the tool compound.
  This is required for all tool compounds and must be variable. It is used to identify the specific
  tool compound that a read maps to (if there are multiple tool compounds with the same
  library tag).

The optional reserved sections are the same between tool compounds and DEL libraries.

An example schema looks something like this:
.. code-block:: json

    "barcode_schema": {
        "library": {
            "tag":  "CCTTGGCACCCGAGAATTCCAATCGCTGA",
            "overhang": "CTA"
        },
        "compound_tag": {
            "tag": "NNNNNNNNNNNNNNNN",
            "overhang": "ACG",
            "error_correction": "hamming_dist_1"
        "extra_tags": {
            "tag": "GGGNNNAGGGNNNATC"
        },
        "preumi": {
            "tag": "AATGCCAGTACG"
        },
        "umi": {
            "tag": "NNNNNNNNNNN"
        }
    }


Tool Compounds
--------------
The tool compounds section is where the actual tool compounds are defined.
This section is a list of tool compound definitions.
Each tool compound definition has only one required field:

- ``compound_id``: A unique identifier for the tool compound.

The other sections are optional (but some are required for decoding):

- ``smiles``: The SMILES string for the tool compound.
- ``tag``: The DNA tag sequence for the tool compound.
  This is required for decoding tool compounds.
  It should be the DNA sequence of the "compound_tag" section of the barcode schema.
  This is what DELi will use to identify the specific tool compound during decoding.


Example Tool Compound Definition
--------------------------------
Here is an example tool compound definition file with two tool compounds defined:
.. code-block:: json
{
    "barcode_schema": {
        "library": {
            "tag":  "CCTTGGCACCCGAGAATTCCAATCGCTGA",
            "overhang": "CTA"
        },
        "compound_tag": {
            "tag": "NNNNNNNNNNNNNNNN",
            "overhang": "ACG",
            "error_correction": "hamming_dist_1"
        }
        "extra_stuff": {
            "tag": "GGGNNNAGGGNNNATC"
        },
        "umi": {
            "tag": "NNNNNNNNNNN"
        }
    },
    "tool_compounds": [
        {
            "compound_id": "ToolCompound_1",
            "smiles": "CCO",
            "tag": "ACGTACGTACGTACGT"
        },
        {
            "compound_id": "ToolCompound_2",
            "smiles": "CCN",
            "tag": "TGCATGCATGCATGCA"
        }
    ]
}

Tool Compound File FAQ
----------------------
- **Q:** My tool compounds have the same library tag but different barcode schemas. Does DELi support this?
    **A:** Not directly. Since DELi first separates reads based on library tags then aligns and decodes them, it
    requires that everything in a given library to have a single barcode schema (for the alignment).
    However, you can get around this by defining multiple tool compound files, each with a unique library tag.
    In this case your library tag will need to include enough of the barcodes from each to cover the part that makes
    the two schemas different.

    In the case that one of your tool compounds is a perfect subset of another, DELi offers no way to resolve this conflict.
    In general, this a bad idea for your DEL experiment as well, hiccups in ligation or synthesis could lead to miscalls
    between the two compounds.

Deli Data Directory
-------------------
Tool compound definition files are stored in the :ref:`DELi Data Directory <deli-data-dir-ref>`,
under the ``tool_compounds/`` sub-directory.
