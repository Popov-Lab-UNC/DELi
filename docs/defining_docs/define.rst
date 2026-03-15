.. _define-dels:

Defining DELs for use in DELi
=============================

In order to use DELi to process your DEL data, DELi needs to know
some information about your DELs. These are called the components of a DEL,
and they include things like libraries, building blocks, building block sets, reactions, tool compounds, etc.
Nearly any DEL can be defined using these two file types (and if
yours cannot, please open an issue on the DELi GitHub repository to add support).

Component IDs
--------------
In DELi each unique component requires an ID. The rules on whether an ID can be shared amoung different components varries based on the component type.
For example building block sets must have unique IDs all other sets (a globally unique ID), but building blocks only need to have unique IDs within the same
building block set (a locally unique ID).

DELi has little restriction how what IDs can look like. The only restriction is that they cannot contain four reserved tokens:
- The ``comp_id_sep`` token (default is "-"), which is used to separate components when :ref:`generating compound IDs <compound_ids>`.
- The ``BB_MASK`` token (default is "###"), which is used to mask building blocks when creating Di/Monosynthon IDs.
- The "." character
- The "," character

Not all IDs enforce this restriction. For example a building block set or library ID can contain the ``BB_MASK``, but as a rule of thumb avoid
using these token in any of your IDs to avoid confusion. DELi will raise an error if, during the loading of all required components, if finds and ID
that violates its rules.

The detailed rules for each type of component are outlined in the documentation for that component.

.. toctree::
    :maxdepth: 2

    defining_libraries
    defining_building_blocks
    defining_enumerators
    defining_reactions
    defining_tool_compounds
    defining_hamming
