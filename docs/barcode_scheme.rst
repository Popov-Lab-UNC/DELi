===============
Barcode Schemas
===============
Barcode Schemas define how the DNA tag (the barcode) is designed. Defining this schema is what allows DELi to take the sequences reads and converts them into information about the library selection results. If you configure the schema you provide DELi incorrectly, it can result in DELi producing an incorrect analysis; DELi can only determine if a schema is impossible, not if it has a typo in the primer sequence. Thus, it is very important to both understand how to write the schema properly but also to double check it for accuracy.

It should be noted that schema are separate from a full DNA barcode; the schema will include variable regions that we expect to be different among reads. The schema is meant to define where those variable regions are with respect to each other AND define any static regions expected to be the same among all reads (for example, a primer region or overhang)

Schema Format
=============

Barcode Sections
----------------
DELi breaks down a barcode schema into several smaller sections. These sections are named and contain infomation about the DNA encoded in that section. Some, like ``library_tag`` are variable regions and are encoded with ``N`` rather than a normal base pair. Variable regions imply that we expect this region to vary during analysis. Others regions, like the ``primer`` region are static and contain explict base pairs that make up that regions. Static regions can contain variable base pairs, for example ``GCGNATNCA`` is a valid `static` base pair. However ``NNNNN`` is not a valid static base pair as it contains only vairable regions.
All regions are either required to be static **OR** variable but never both. Some sections are Optional, and can be set to ``null`` or and empty string if not utilized.


.. note::
    DELi supports de-multiplexing during analysis. Most NGS software support this process during sequencing, so it is not required during analysis, hence the index regions are optional

Below are a descriptions each section that DELi supports

.. csv-table::
   :header: "Section Name", "Static/Variable", "Optional/Required", "Description"
   :widths: 10, 7, 10, 30

   "pre-index","Static","Optional","region before the index tag of the experiment"
   "index","Static","Optional","the index tag of the experiment (if not de-multiplexed)"
   "primer","Static","Optional","the primer region of DNA tag"
   "library_tag","Variable","Required","the region with the library DNA tag"
   "library_tag_overhang","Static","Optional","a static overhang region after the library_tag"
   "bb1_tag","Variable","Required","the region with the first cycle building block DNA tag"
   "bb1_tag_overhang","Static","Optional","a static overhang region after the bb1_tag"
   "bb2_tag","Variable","Required","the region with the second cycle building block DNA tag"
   "bb2_tag_overhang","Static","Optional","a static overhang region after the bb2_tag"
   "bbN_tag","Variable","Optional","the region with the Nth cycle building block DNA tag"
   "bbN_tag_overhang","Static","Optional","a static overhang region after the Nth tag (if Nth tag is included)"
   "pre-umi","Static","Optional","region before the umi tag of the experiment"
   "umi","Variable","Required","the randomized UMI region of the tag"
   "closing_primer","Static","Optional","the closing primer region after the UMI"

.. note::
    DELi requires at least 2 building blocks (``bb1_tag`` and ``bb2_tag``) but supports (optional) additional building blocks up to 99 (``bb3_tag`` - ``bb_99_tag``).
    All these tags can have overhangs (optionally) as well

Schema File Format
------------------
A JSON file is used to encode barcode schemas.
The format requires three named elements:
 * ``id``; the name of the schema
 * ``sections``; a dictionary of all the sections and their DNA tag definition"

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
DELi assumes the order of the sections is the same as the order they are listed. If this is not the case and don't want to rearrange the sections, you can also include an ``order`` dictionary in the configuration file. If a order dictionary is included, all sections **must** have a defined order, orders cannot be assumed.

If a section is ``null``, the order must allso be ``null``. Ordering starts from 0. Order does not need to be continuous; just organized such that sections that come first must have a lower value than section that come later.

.. code-block::

    {
        "id": "MyExampleSchema",
        "sections": {
            "pre-index": "AGCGCAT",
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
            "closing_primer": "GCTGCACTTGCT"
        }
        "order": {
            "pre-index": null,
            "index": 0,
            "primer": 1,
            "library_tag": 8,
            "library_overhang": 9,
            "bb1_tag": 2,
            "bb1_overhang": 3,
            "bb2_tag": 4,
            "bb2_overhang": 5,
            "bb3_tag": 6,
            "bb3_overhang": 7,
            "pre-umi": 10,
            "umi": 11,
            "closing_primer": 12
        }
    }

Barcode Schema files should be save as a json with the same name as the id (i.e. ``MyExampleSchema.json``)
and saved in the ``$BARCODE_SCHEMA_FOLDER`` folder

.. note::
    if a section is ``null`` and the order is not, the schema is still parseable; that order will just ignored since that section will be ignored. This is not recommend since it makes the config file hard to interpret.

Loading a Schema
================
A schema can be loaded using the following code::

    from deli.barcode import BarcodeSchema
    barcode_schema = BarcodeSchema.load_from_json("MYFILE.json")

This is a quick and easy way to check if you schema is properly written: if this returns an error it means that something is wrong. The error message should give some guidance on what to fix
