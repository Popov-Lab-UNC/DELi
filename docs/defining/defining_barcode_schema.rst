========================
Defining Barcode Schemas
========================
Barcode schemas are used to define the format of the DNA barcode that is attached
to the chemicals of a DEL.
DELi requires this schema into order to figure out which part of the DNA tag
to need to read for calling to occur. DELi expects one schema per file and they look like this:
::
    {
        "id": "MyExampleSchema",
        "sections": {
            "index": "NNNNNNNNNN",
            "primer": "CCTTGGCACCCGAGAATTCCA",
            "library": "NNNNNNNN",
            "library_overhang": "AG",
            "bb1": "NNNNN",
            "bb1_overhang": "GT",
            "bb2": "NNNNN",
            "bb2_overhang": "GA",
            "bb3": "NNNNN",
            "bb3_overhang": "TT",
            "pre-umi": "AATGCCAGTACG",
            "umi": "NNNNNNNNNNN"
        }
    }

Schema Sections
===============
Barcodes are broken up into several section (called ``BarcodeSections`` in DELi).
Sections can either be "static", where all the nucleotides are known and the same
for every compound in the library, or "variable", where the nucleotides change for
different compounds. All sections are either static OR variable, they cannot be both
(though overhang regions can give a variable region some static portions).
Variable regions are represented as 'N' rather than one of the 4 base pairs (ATGC).

DELi support eight (8) specific types of barcode sections.
Only a handful are required to always be in the schema, the rest are optional.
Unless otherwise stated, all sections should only appear once (if at all) in the schema:

- pre-index
This section is static, and always comes before the index region.
It is not required, and can be left out if you don't have an index or
anything before it.

- index
If you require DELi to demultiplex your reads for you, you can include the index
region. This region is variable, and is used to determine which index/sample this
read came from.

- primer
The primer region is one of the only required regions that must be in your schema.
It static, and is used to locate a valid read to extract for alignment and calling.

- library
This region is also required, but is variable.
It is the region used to determine which library this read comes from.
This is crucial to allow for correct calling of the building blocks.

- bb##
This is the only section where multiple copies can be present,
as long as the numbers are different.
DELi only supports up to 99 bb (which would be a 99 cycle library and far bigger
than practical anyway).
All bb sections are variable.
DELi will require that a bb1 and bb2 section are present,
but all other from bb3-bb99 are optional, and only necessary if you have those cycles.
The only other limitation is that bb sections must appear in numerical order.
So you cannot start with bb2 and then next be bb5 and bb16. It must be bb1, bb2 and bb3. If you do not follow this DELi can have weird behavior

- pre-umi
This is an optional static region that comes before the umi region.
Generally it is techincally part of the closing primer, but it is split up to
help isolate the umi region

- umi
This is the last required section and is variable. It is used to remove the
noise added during random sampling of PCR.

- closing
This is a final static region behind the umi region.

Most DEL barcode design follow is exact order, but if yours does not,
sections can be placed in any order.
The only exception is the Index, Library and Primer sections, which must come before all bb sections.
This is more a limitation of DEL than of DELi, as these section should be present before synthesis or attached after selection occurs, thus it is not possible to
get the bb sections before them.

DELi iterpurts the order of the sections listed in the file as the order they
appear in the barcode, so for the example defined above, the full barcode would look like:
"NNNNNNNNNNCCTTGGCACCCGAGAATTCCANNNNNNNNAGNNNNNGTNNNNNGANNNNNTTAATGCCAGTACGNNNNNNNNNNN".

Overhangs
=========
Often times, after variable regions there are short 3-5 bp static regions
(these are to help with DNA ligation). To accommodate for this, every section
has a corispoding "_overhang" section. These must directly follow the section they
are named after.
For example
::
    {
        "sections": {
            "library": "NNNNNNNN",
            "library_overhang": "AG",
            "bb1": "NNNNN",
            "bb1_overhang": "GT",
            "umi": "NNNNNNNNNNN"
        }
    }

Is a valid section, but
::
    {
        "sections": {
            "library": "NNNNNNNN",
            "library_overhang": "AG",
            "bb1": "NNNNN",
            "umi": "NNNNNNNNNNN"
            "bb1_overhang": "GT",
        }
    }

is not. This is because DELi expect the section to be listed in the order they
appear. It makes no sense for a overhang to be after a umi region if it was for bb1

Hamming
=======
DELi can handle hamming encoded tags. The barcode schema defines which sections are hamming
encoded AND which hamming code they are using (see Defining Hamming Codes). This specification is
done in a seperate section in the json file keyed as 'hamming'
::
    {
        "id": "MyExampleSchema",
        "sections": {
            "index": "NNNNNNNNNN",
            "primer": "CCTTGGCACCCGAGAATTCCA",
            "library": "NNNNNNNN",
            "library_overhang": "AG",
            "bb1": "NNNNN",
            "bb1_overhang": "GT",
            "bb2": "NNNNN",
            "bb2_overhang": "GA",
            "bb3": "NNNNN",
            "bb3_overhang": "TT",
            "pre-umi": "AATGCCAGTACG",
            "umi": "NNNNNNNNNNN"
        },
        "hamming": {
            "bb2": "hamming7"
            "bb3": "hamming9"
        }
    }

In the hamming section, just specify the name of the section and attach it the correct hamming code used.
Just like with any DELi Data Dir compatable files, the hamming code could be the path to the hamming code file OR the name of it if it is in the hamming sub dir in the deli data dir.

DELi Data Dir
=============
Barcode schemas can be saved in the DELi Data Dir sub-dir named "barcodes"
