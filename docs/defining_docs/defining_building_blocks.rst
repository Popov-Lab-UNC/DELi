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
It can be both, but one must be present (otherwise this list of building blocks
is not useful to DELi).
The order of these columns does not matter.
As expected:

- the ``id`` column holds the name/ID of a given building block,
- the ``tag`` column hold the DNA tag associated with that building block (must be a valid DNA sequence),
- the ``smiles`` column holds the SMILES string for that building block.

Any other columns are ignored by DELi and it will not cause any issues to include them.

.. note::
    If your DEL building block regions have overhangs
    the overhang nucleotides **should not** be part of
    the tag included in the building block set file

An example building block set file might look like:
::
    id,tag,smiles
    BB1,AGCTG,CCC
    BB2,GCTAG,CCCCC
    BB3,GGCTT,CCCCOC
    ...

Naming building block sets
^^^^^^^^^^^^^^^^^^^^^^^^^^
The base name of the file is used as the name of the building block set by DELi.
If my building block file is ``path/to/my_building_blocks.csv``, then the ID/name
of this set is my_building_blocks. DELi does not require that building block sets
have unique IDs/names during runtime (as you might want to load copies of the same
set). To avoid confusion, it is best to use unique names for each building block set
when saving them. The :ref:`DELi Data Directory <deli-data-dir-ref>` is a good place to save these files

Uniqueness of Building Blocks
-----------------------------
DELi will assert that the ``id`` column is unique for every building
of a given building block set. While IDs between sets can be the same,
they cannot be the same within a given set.
This is also true for the ``tag`` column, if provided. All tags must be unique
within a given building block set, if not DELi will raise an error
when loading the building block set.
Unlike the other two columns, the ``smiles`` column does not have to be unique.
While it would make sense (in most cases) to have a unique SMILES, there are cases
where multiple building blocks can have the same SMILES string to help with reducing
DEL noise. Beware, this means DELi will not check for uniqueness of the SMILES column,
so best to double check that SMILES are correct for the set.

Saving in the DELi Data Directory
---------------------------------
If you have configured the :ref:`DELi Data Directory <deli-data-dir-ref>`,
you should save the building block set files in the ``building_blocks`` sub-directory.
This way, when using DELi, you can reference the building block sets by their name, rather
than having to know the exact path to the file.
