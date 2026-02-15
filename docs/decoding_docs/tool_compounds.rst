.. _tool-compounds::

==============
Tool Compounds
==============

In some cases, a DEL selection may include DNA tagged compounds that exist outside of the DELs used.
A common example is trying to recover a known binder when validating your selection protocol.
DELi refers to these "tool compounds" and implements custom logic to handle these compounds during decoding.

.. note::
    Tool compounds could also refer to compounds added to the selection condition (but not tagged).
    DELi does not use this definition. A tool compound in DELi **must** have a DNA tag otherwise it
    should be considered part of the selection condition.

Decoding Tool Compounds
-----------------------
Obviously the main reason to add tool compounds is to be able to decode them and track their counts during decoding.
DELi will attempt to decode tool compounds in the same way as DEL compounds. However you need to declare them to DELi
if a special way, since they are not traditional DELs. There are two ways to do this in DELi: "doping" or "adding".

.. _doped-compounds::

Doping Tool Compounds
---------------------
Doped tool compounds are tool compounds that are "mixed" into a DEL. This means it shares the same barcode schema
as that DEL. It does not have to share the same reaction (or even number of building blocks). As long as the
doped compounds DNA tag follows the barcode schema (including share the DELs library tag sequence!),
the compound can be doped into the DEL. As a result, doped compounds are defined in the library definition file
directly. You can read about how to define a doped tool compound in the :ref:`defining DELs docs <doped-compounds-sec>`.

In nearly 99.99% of cases you should never do this, and it might get removed some day. This is because it rarely happens
in practice (outside of my old academic lab at least). "added" compounds cover nearly all cases where tool compounds
are used. The only case that doping would be useful is if you tool compound's DNA tag has no regions that could
differentiate it from another DELs DNA schema.

.. _tool-compounds-added:

Adding Tool Compounds
---------------------
Added tools compounds is a simple as telling DELi that they were used providing the correct configuration file
to let DELi know how to decode them.

You can read more about how to specify added tool compounds in a :ref:`selection file <selection-file-docs>`
and how to define the :ref:`tool compound configuration <defining-tool-compounds>`.

How DELi handles tool compounds
-------------------------------
Under the hood DELi will load in doped tool compounds as part of the DEL library. They will be detached from the library
size and library enumerator (since they exist outside the combinatorial space of the DEL), but their section tags
will be added the tag map to allow for decoding.

Added tool compounds are a bit different. These will get loaded as a "library" of one compound. DELi needs to do this
since it operates by first determining the library a read belongs to. Since there is only one compound in a "tool compound library"
DELi doesn't need to decode it post library demultiplexing. However, it will still do a quick alignment check to make
sure that the ``tool_compound_ref`` section is present (with some error tolerance) to avoid miscall. The UMI will also be
extracted. See the :ref:`defining tool compounds docs <defining-tool-compounds>` for more details.
