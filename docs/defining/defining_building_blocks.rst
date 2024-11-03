=============================
Defining Building Blocks Sets
=============================

Since (most) DELs are made up of thousands of building blocks,
defining each one as a separate file (like for other components)
would be too tedious. Instead, they are defined as single rows
if a large "set" of building blocks. This makes more sense
intuition wise as well, as DELs are composed of "cycles" where
an entire set of building blocks is used for that cycle.

As a results, you will only ever have to define a full building
block set (a ``BuildingBlockSet`` object in DELi). DELi can
only load building blocks as sets. It cannot load single
building blocks.

File Format
===========
The main file format is that of a comma separated value (CSV) file.
A header line is expected, naming the columns
There are two required columns: ``id`` and ``tag``.
The order of these columns doesn not matter.
As expected the ``id`` column holds the name/ID of a given building block.
The ``tag`` column hold the DNA tag associated with that building block

.. note::
    If your DEL building block regions have overhangs
    the overhang nucleotides **should not** be part of
    the tag included in the building block set file

An example building block set file might look like:
::
    id,tag
    BB1,AGCTG
    BB2,GCTAG
    BB3,GGCTT
    ...

Including SMILES
----------------
An optional SMILES column can also be included.
The name of this column should be ``smiles`` and
would hold the SMILES string for the given tag.

.. note::
    DELi assumes that these SMILES for each building block
    are the raw chemical used in the actual reactions to
    create the library. This includes any protecting groups
    or other reactive handles that might not make it into
    the final compound structure

Saving in the DELi Data Directory
=================================
For easy organiation, building block sets can be saved
as files in the ``building_blocks`` sub-directory in a
DELi Data Directory. This allows for the loading of the
set by just using the name. The name of the building
block set should be the same as the name of the file.
E.g. set ``BBA_DEL1`` should be saved as ``BBA_DEL1.csv``
