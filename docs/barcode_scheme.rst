===============
Barcode Schemas
===============
Barcode Schemas define how the DNA tag (the barcode) is designed. Defining this schema is what allows DELi to take the sequences reads and converts them into information about the library selection results. If you configure the schema you provide DELi incorrectly, it can result in DELi producing an incorrect analysis; DELi can only determine if a schema is impossible, not if it has a typo in the primer sequence. Thus, it is very important to both understand how to write the schema properly but also to double check it for accuracy.

It should be noted that schema are separate from a full DNA barcode; the schema will include variable regions that we expect to be different among reads. The schema is meant to define where those variable regions are with respect to each other AND define any static regions expected to be the same among all reads (for example, a primer region or overhang)

Schema Format
=============

Barcode Sections
----------------
DELi breaks down a barcode schema into several smaller sections. These sections are
named and contain information about the DNA encoded in that section.
Not all sections are equal in the way we need to treat them.
Some are required, like the ``primer`` or ``umi`` regions.
Some allow for variable nucleotides (N) like ``library``
while others require the DNA is static, like ``primer``.
This is controlled for in DELi using the ``BarcodeSectionTraits``.
This defines the traits that correspond to each region.

Further, each region can have what we call an "overhang".
This is discussed in more detail a bit later,
but are never required and are always static.

.. note::
    DELi supports de-multiplexing during analysis. Most NGS software support this
    process during sequencing, so it is not required during analysis,
    hence the index regions are optional

Below are a descriptions each section that DELi supports

.. csv-table::
   :header: "Section Name", "Static/Variable", "Optional/Required", "Description"
   :widths: 10, 7, 10, 30

   "pre-index","Static","Optional","region before the index tag of the experiment"
   "index","Static","Optional","the index tag of the experiment (if not de-multiplexed)"
   "primer","Static","Required","the primer region of DNA tag"
   "library","Variable","Optional","the region with the library DNA tag"
   "bb<N>","Variable","Required(ish)","the region with the <N>th cycle building block DNA tag. DELi requires at least two of these are passed"
   "pre-umi","Static","Optional","region before the umi tag of the experiment"
   "umi","Variable","Required","the randomized UMI region of the tag"
   "closing","Static","Optional","the closing (primer) region after the UMI"

.. note::
    DELi requires at least 2 building blocks (``bb1_tag`` and ``bb2_tag``) but supports (optional) additional building blocks up to 99 (``bb3_tag`` - ``bb_99_tag``).
    All these tags can have overhangs (optionally) as well

Schema File Format
------------------
A JSON file is used to encode barcode schemas.
The format requires three named elements:
 * ``id``; the name of the schema
 * ``sections``; a dictionary of all the sections and their DNA tag definition"

Barcode Schema files should be save as a json with the same name as the id (i.e. ``MyExampleSchema.json``)
and saved in the ``$BARCODE_SCHEMA_FOLDER`` folder

Below is an example file

.. code-block::

    {
        "id": "MyExampleSchema",
        "sections": {
            "pre-index": null,
            "index": "NNNNNNNNNN",
            "primer": "AGCTAGCTAGCTAGCTAGCT",
            "library_tag": "NNNNNNNN",
            "library_overhang": "AGT",
            "bb1_tag": "NNNNN",
            "bb1_overhang": "GCT",
            "bb2_tag": "NNNNN",
            "bb2_overhang": "ACT",
            "bb3_tag": "NNNNN",
            "bb3_overhang": "TGC",
            "pre-umi": "GCTGCT",
            "umi": "NNNNNNNNNN",
            "closing_primer": "GCTGCACTTGCT",
        }
    }


While not required, sections are are optional and not included in your schema should be encoded as ``null`` to express explicitly in the config file that they are not used

Schema Ordering
---------------
DELi assumes the order of the sections is the same as the order they are listed.
Thus, in the ``sections`` part, list them in the correct order as they would
appear reading the DNA tag left to right

Loading a Schema
================
A schema can be loaded using the following code::

    from deli.barcode import BarcodeSchema
    barcode_schema = BarcodeSchema.load_from_json("MYFILE.json")

This is a quick and easy way to check if you schema is properly written: if this returns an error it means that something is wrong. The error message should give some guidance on what to fix

Overhangs
=========
Each valid barcode section can have an overhang section associated with it.
These *MUST* be named ``<section_name>_overhang``,
so for example ``library_overhang`` for the ``library`` section.
They also *MUST* come directly after the region they are associated with.
If they do not, DELi will complain and through ``BarcodeSchemaErrors``
