.. _enumeration-docs:

===========
Enumeration
===========

DELi includes support for enumerating DNA-encoded libraries (DELs).
All it requires is that your libraries files have the necessary :ref:`reaction information <reaction-sec-ref>`
and that your :ref:`building block sets have a SMILES column <defining_building_blocks>`.
See :ref:`the defining DELs docs <define-dels>` for more info on how to define these files.

.. note::
    Enumeration does not require any DNA tag information, thus information about
    the DNA tags and schema are not required.

The enumerator performs the following:

1. Extracts reaction information from the library files
2. Builds a reaction workflow (using RDKit)
3. Runs the reaction for every combination of building blocks for all cycles

How Enumeration Works
----------------------
DELi treats enumeration as a directed acyclic graph (DAG) of reactions.
Each node in the graph is a reaction step,
and each edge is a possible next reaction. This allows for complex
reaction schemes to be represented, including branching paths and multiple
routes through the synthesis process. For example, when adding cycle 2, some
building blocks might undergo a esterification reaction, while others might undergo
an amide coupling. DELi allows you to specify these different routes by provide two
different steps, one for the esterification and one for the amide coupling, and then
specifying which building blocks (using :ref:`building block subsets <TODO>`)
go through which step.

There is also support for reaction priority, where you can specify a preferred
priority of reactions to try for any given reaction step. Say you want to try and
find the primary amine first, then a secondary amine if there is no primary. DELi's
enumeration module can handle this we well. See the :ref:`reaction docs <reaction-sec-ref>`
for more info on how to configure this in a file.

.. warning::
    As you can imagine, enumeration can be extremely memory/compute intensive for large libraries.
    Careful consideration should be taken about when and how to enumerate a library.

Enumeration Python API
----------------------
The top-level object handling enumeration is the ``enumerator`` object.
It is possible to build these objects manually, but the easiest way is to
load a ``Library`` object from a library JSON file and then call the
``enumerate`` method on it *or* extract the enumerator object using the
``get_enumerator`` method.
However, below is an example of how to build on yourself.

First you need to create/load the building block sets that will be enumerated.

.. code-block:: python

    from deli.dels.building_block_set import BuildingBlockSet
    bb_set1 = BuildingBlockSet.load("path/to/bb_set1.csv")
    bb_set2 = BuildingBlockSet.load("path/to/bb_set2.csv")
    bb_set3 = BuildingBlockSet.load("path/to/bb_set3.csv")

    for bb_set in [bb_set1, bb_set2, bb_set3]:
        assert bb_set.has_smiles()  # Ensure SMILES are present in the building block sets

Then you need to create the ``ReactionTree`` object that defines the reaction scheme.
This is done by first listing out all the reaction steps.

.. code-block:: python

    from deli.enumeration.reaction import ReactionStep, Reaction, ReactionTree, BBSetReactant, ProductReactant, StaticReactant
    step1 = ReactionStep(
        step_name="Cycle 1",
        step_id="1",
        reaction=Reaction("[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[N:2]")
        reactants= [BBSetReactant(bb_set_reactant_id="bb_set1"), StaticReactant(smiles="[13C]N")],
    )
    step2 = ReactionStep(
        step_name="Cycle 2",
        step_id="2",
        reaction=Reaction("[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[N:2]")
        reactants= [ProductReactant("1"), BBSetReactant(bb_set_reactant_id="bb_set3")],
    )
    step3 = ReactionStep(
        step_name="Cycle 3",
        step_id="3",
        reaction=Reaction("[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[N:2]")
        reactants= [ProductReactant("2"), BBSetReactant(bb_set_reactant_id="bb_set3")],
    )

    reaction_tree = ReactionTree(steps=[step1, step2, step3])

Now we can define the enumerator object itself.

.. code-block:: python

    from deli.enumeration.enumerator import DELEnumerator
    enumerator = DELEnumerator(
        reaction_tree=reaction_tree,
        building_block_sets=[bb_set1, bb_set2, bb_set3]
    )

Enumeration with DELi
---------------------
DELi provides enumeration functionality through the ``Library`` class.
Objects of this class can be enumerated using the ``enumerate`` method:

.. code-block:: python

    from deli.dels.library import Library
    library = Library.load("path/to/library.json")
    for compound in library.enumerate():
        type(compound)  # An `EnumeratedDELCompound` object
        del_id = compound.compound_id  # Access DEL ID of enumerated compound
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

You can extract the ``Enumerator`` object from the ``Library`` object, but doing so detaches it from the
info about the library. This means it cannot create ``DELCompound`` objects with the proper library
information and DEL IDs. Instead it just tells you which building blocks were used and generated product

.. code-block:: python

    enumerator = library.get_enumerator()
    for result in enumerator.enumerate():
        bbs = result[0] # List of building block used to make this compound
        mol = result[1]  # RDKit mol object of the enumerated compound


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
DELi's CLI includes a command for enumerating DELs from a library file: ``deli enumerate``.
All you need to do is provide the path to a library JSON file. See the :ref:`CLI docs <deli-enumeration-cli-docs>`
for more info.
