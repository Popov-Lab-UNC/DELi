.. _defining_building_blocks:

=============================
Defining Building Blocks Sets
=============================

Since (most) DELs are made up of thousands of building blocks,
defining each one as a separate file (like for other components)
would be too tedious. Instead, they are defined as single rows
of a large "set" of building blocks. This makes more sense
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
There is only one required column: ``id``
Additionally, you must include either a ``tag`` column or a ``smiles`` column.
It can be both, but one must be present (otherwise it is just of list of ids).
The order of these columns does not matter.

.. note::
    DELi is case sensitive when it comes to column names.

Below are the possible columns:

- ``id`` column holds the name/ID of a given building block,
- ``tag`` column hold the DNA tag associated with that building block (must be a valid DNA sequence),
- ``smiles`` column holds the SMILES string for that building block.
- ``subset_id`` column holds the subset ID for that building block.
  This is used when defining building block :ref:`subsets for reactions <reaction-bb-subset-docs>`.

Any other columns are ignored by DELi; it will not cause any issues to include them.

.. note::
    If your DEL building block regions have overhangs
    the overhang nucleotides **should not** be part of
    the tag included in the building block set file

An example building block set file might look like:
::
    id,tag,smiles,subset_id
    BB1,AGCTG,CCC,1
    BB2,GCTAG,CCCCC,1
    BB3,GGCTT,CCCCOC,2
    ...

.. _defining_building_block_set_names:

Naming building block sets
^^^^^^^^^^^^^^^^^^^^^^^^^^
The base name of the file is used as the name of the building block set by DELi.
If my building block file is ``path/to/my_building_blocks.csv``, then the ID/name
of this set is ``my_building_blocks``. DELi requires that all building block
set *within a given Library* have unique names.
Some keywords are also reserved by DELi and cannot be used as building block set names,
such as ``linker`` or ``scaffold``. DELi will raise an error if you try to load
a building block set with a reserved name.

Uniqueness of Building Blocks
-----------------------------
DELi will assert that the ``tag`` column (if provided) is unique for every building
of a given building block set. While tags between sets can be the same,
they cannot be the same within a given set.

To support the possible multi-tagging or building blocks, DELi does not require
the ``id`` column to be unique. However, if you include a ``smiles`` column it
will assert that the SMILES is the same for all building blocks with the same ``id``.
If there is a ``subset_id`` column it will also assert that all building blocks
with the same ``id`` have the same ``subset_id``.

DELi will warn you if it detects two identical building blocks in a file
that (all elements are the same). It will then just keep one of them and
continue on.

It is worth noting here that this is not the exact same as how DELi handles
``BuildingBlockSet`` objects in python. It will always assert that
all ``BuildingBlock`` objects in a ``BuildingBlockSet`` have unique. This is
because ``TaggedBuildingBlock`` objects can handle an arbitrary number of tags
for the same building block, something that you cannot really do in a CSV file.

.. warning::
    You must provide the **exact same** SMILES. Even if the SMILES are
    chemically equivalent but written differently, DELi will raise an error.

Saving in the DELi Data Directory
---------------------------------
If you have configured the :ref:`DELi Data Directory <deli-data-dir-ref>`,
you should save the building block set files in the ``building_blocks`` sub-directory.
This way, when using DELi, you can reference the building block sets by their name, rather
than having to know the exact path to the file.
