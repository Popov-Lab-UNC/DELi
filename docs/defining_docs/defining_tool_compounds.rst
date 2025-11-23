.. _defining-tool-compounds::

=======================
Defining Tool Compounds
=======================

:ref:`Tool compounds <tool-compounds>` are defined as JSON files They are
defined by their library tag, single compound ID (and possible SMILES structure) directly.
All they need if for you to specify is the following fields:

- ``library_tag``: The library tag of the tool compound. This must be unique from all other
  library tags in the selection.
- ``compound_id``: The unique ID of the tool compound.
- ``smiles``: (Optional) The SMILES structure of the tool compound.

Example Tool Compound Definition
--------------------------------
Here is an example tool compound definition file:
.. code-block:: json

    {
        "library_tag": "AGCTAGCTAGCT",
        "compound_id": "TOOL_COMP_001",
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"
    }

Deli Data Directory
-------------------
Tool compound definition files are stored in the :ref:`DELi Data Directory <deli-data-dir-ref>`,
under the ``tool_compounds/`` sub-directory.
