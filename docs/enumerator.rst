===========
Enumeration
===========

DELi includes support for enumerating DNA-encoded libraries (DELs).
All it requires is that your libraries files have the necessary reaction information
and that your building block sets have a SMILES column. See :ref:`the defining DELs docs <define-dels>`
for more info on how to define these files.

.. note::
    Enumeration does not require any DNA tag information, thus this data is not
    required in any of the files.

The enumerator performs the following:

1. Extracts reaction information from the library files
2. Builds a reaction workflow (using RDKit)
3. Runs the reaction for every combination of building blocks for all cycles

As you can imagine, this can be extremely memory/compute intensive for large libraries.
Careful consideration should be taken about when and how to enumerate a library.

Enumeration with DELi
---------------------
DELi provides enumeration functionality through the ``Library`` class.
Objects of this class can be enumerated using the ``enumerate`` method:

.. code-block:: python

    from deli.dels.library import Library
    library = Library.load("path/to/library.json")
    for compound in library.enumerate():
        type(compound)  # An `EnumeratedDELCompound` object
        smi = compound.smi  # Access SMILES of enumerated compound
        mol = compound.mol  # Access RDKit mol object

By default, the ``enumerate`` method will *not* fail if a reaction cannot be carried out for a set of blocks.
Instead it will return a ``DELCompound`` object rather than an ``EnumeratedDELCompound`` object. This means
the DEL will lack a `smi` and `mol` attribute. If you want the method to fail instead, you can just set
the ``fail_on_error`` argument to ``True``:

.. code-block:: python

    for compound in library.enumerate(fail_on_error=True):
        # If a reaction fails, an EnumerationRunError will be raised
        pass

Note that DELi will tell you which library and which set of building blocks caused the reaction to fail
if you want to try and figure out the issue. If you just want a list of pure ``EnumeratedDELCompound``s
You can also just drop failed compounds by setting ``drop_failed=True``:

.. code-block:: python

    for compound in library.enumerate(drop_failed=True):
        # Only `EnumeratedDELCompound` objects will be returned
        assert type(compound) is EnumeratedDELCompound
        smi = compound.smi
        mol = compound.mol

``fail_on_error`` and ``drop_failed`` are mutually exclusive, so only one of them can be True or DELi will
throw an exception.

Enumerating to a file
---------------------
More often than not, you will want to save the result of enumeration to a file.
While you can just generate all the enumerated compounds in memory and then write them to a file,
that can be very memory intensive. Instead you can use the ``enumerate_to_csv_file`` method:

.. code-block:: python

    library.enumerate_to_csv_file(
        out_path="enumerated_library.csv",
        separator=","
        dropped_failed: bool = False,
        fail_on_error: bool = False,
        use_tqdm: bool = True,
    )

This will write the results of enumeration directly to a CSV file *on the fly* to keep the
memory usage low. The CSV file will have the following columns:

- `DEL_ID`: the DEL id of the compound
- `SMILES`: the SMILES of the compound
- `LIB_ID`: the id of the library the compound belongs to
- *for each building block set in the library* :
    - `BB<cycle#>_ID`: the ids of the building blocks used to build the compound
    - `BB<cycle#>_SMILES`: the smiles of the building block sets used to build the compound

Single Compound Enumeration
---------------------------
You can also enumerate a single compound by providing the building block object *or*
the IDs of the building blocks to use:

.. code-block:: python

    bb_cycle1 = library.bb_sets[0].get_bb_by_id("BB_1")
    bb_cycle2 = library.bb_sets[1].get_bb_by_id("BB_234")
    bb_cycle3 = library.bb_sets[2].get_bb_by_id("BB_624")

    compound_a = library.enumerate_by_bbs([bb_cycle1, bb_cycle2, bb_cycle3])
    compound_b = library.enumerate_by_bb_ids(["BB_1", "BB_234", "BB_624"])

    assert compound_a.smi == compound_b.smi  # will be True

This can be useful if you only want a handful or specific subset of the DEL to be
enumerated.

Command-Line Interface
----------------------
DELi's CLI include a command for enumerating DELs from a library file: ``deli enumerate``.
All you need to do is provide the path to a library JSON file. See the :ref:`CLI docs <deli-enumeration-cli-docs>`
for more info.
