.. _tool-compounds::

==============
Tool Compounds
==============

In some cases, compound many be used that are not part of the DEL itself, but still have a DNA tag.
These are often referred to as "tool compounds" and can be used for a variety of purposes,
such as positive controls, negative controls, or internal standards. For example, a known binder
to the target with a reference Kd may be spiked into the DEL selection to help calibrate the results.
More often they are used as a quality control measure to ensure that selection, sequencing and decoding
is working as expected.

Decoding Tool Compounds
-----------------------
Since tool compounds have a DNA tag, it has a section that is treated as the "library" portion of the tag.
This allows us to detect it during decoding. This library tag needs to be unique to the tool compound,
that way it can be easily identified from other libraries in the selection. It does not need to have the
same barcode design as any other library in the selection, nor does it need to have any other barcode section
besides the library tag (unlike DELs which need to have building block barcode sections). It can have them, but
this will be ignored during decoding as each tool compound is a single compound, so once the library tag is called,
DELi know exactly which compound it is.

The reason for this strictness of unique library tags for tool compounds is because they exist outside of the
combinatorial DEL. Adding them to the DEL would brake the assumptions DELi makes about compounds in a given DEL
(that they are made up of building blocks linked by a specific reaction scheme).
By making them unique single compound libraries we can avoid violating any
assumptions DELi makes about the combinatorial nature of other DELs

.. _doped-compounds::

Doping compounds into selections
--------------------------------
Sometimes you tool compounds do have the same library tag as a DEL library. This is most common
when "doping" in a compound to a specific DEL. Overall, this "doping" is not a very common situation.
Thus, DELi is not well equipped to handle this cases.
However, DELi support a "doping" mode for decoding that will allow some flexability here.

In doping mode, you must specify in the DEL library definition the tag, id (and SMILES if available) of the compound
being doped in. In this case the doped compound **must have** the same barcode scheme as the DEL it is being doped into.
This means if must have the same sections (and the same base pair lengths of those sections) in the same order as every
other compound in that DEL. DELi cannot check this for you (only you know what its barcode looks like) so make sure this
is correct to avoid errors and miscalls. DELi will also need to know the building block tags for each of the building block
cycles in the doped tag. DELi will assert that these are unique from the other building blocks in the DEL to avoid miscalls.
DELi will only call this doped compound if all three of the building blocks are called correctly and correspond to the doped compound.

Tool Compound vs Doped Compound
-------------------------------
Under the hood, DELi treats tool and doped compounds as ``ToolCompound`` objects that result is a
``DecodedToolCompound``. The main difference is that during decoding, doped compound are decoded
as an extra part of a DEL library, while tool compounds are decoded as their own separate library.
Tool compounds are thus much faster to decode compared to doped compounds, which can slow down the
decoding process as DELi needs to now check that each DEL compound lacks the doped compound's building
block tags (as that is now an invalid call). It is recommended to use tool compounds whenever possible
over doped compounds to avoid these issues.

Workarounds
^^^^^^^^^^^
Below are some times to help workaround issues to make your doped compounds valid tool compounds.


- My tool compound has a identical barcode and library tag as a DEL in my selection

Since doped compounds must have unique building block tags from the DEL they are being doped into,
you can simply make the library tag for this this compound include the first (or all) of the building
blocks.

- My doped compound has an identical library tag to a DEL, but different tag design.

Much like above you can just make the library tag include the unique parts of the tag design
to create a valid tool compound.

- My doped compound has the same barcode, library tag and (some) building block tags as a DEL.

Make the library tag the full tag of the doped compound. This will make it unique from the DEL.
If you share all the same building block tags your experiment itself if broken. You have added
two different compounds with the same DNA barcode to your selection; Only god can tell them apart now.

Adding Tool/Doped Compounds to Decoding Runs
--------------------------------------------
In your :ref:`decoding selection file <decoding-run-file-docs>` there is a section called ``tool_compounds``.
Here you can list out all the tool and doped compounds that were added for this specific selection.
Each entry in this list points to a tool compound that was used, much like the libraries section points
to DEL libraries that were used. Tool compounds are part of the :ref:`DELi Data Directory <deli-data-dir-ref>`

Defining Tool/Doped Compounds
-----------------------------
See the :ref:`defining DELs docs <defining_libraries>` for more information on how to define a doped compounds.
See the :ref:`defining tool compounds docs <defining-tool-compounds>` for more information on how to define tool compounds.
