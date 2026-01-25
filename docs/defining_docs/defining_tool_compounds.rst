.. _defining-tool-compounds::

=======================
Defining Tool Compounds
=======================

:ref:`Tool compounds <tool-compounds>` are defined as JSON files.
Like with :ref:`DEL library definition files <defining_libraries>`, these files
can optionally implement a barcode schema for decoding.

Other than the ``barcode_schema`` key, tool compound definitions also have

- ``compound_id``: A unique identifier for the tool compound.
- ``smiles``: (optional) The SMILES string for the tool compound.

.. _tool-compound-barcode::

Barcode Schema
--------------
Like the :ref:`DEL library definition files <defining_libraries>`, they can
implement a barcode schema. You can :ref:`read the docs <barcode-sec-ref>` on how these are defined.
Barcode schemas are not required if you are just loading a tool compound
not for decoding, but is required if you want to decode tool compounds (just like for
libraries). However, the required reserved sections are a bit different for tool compounds.

There are only three required sections in the barcode schema for tool compounds:

- ``library`` or ``tool``: This section defines the 'library' tag for the tool compound.
  The naming here is a bit confusing, as tool compounds don't actually have a library.
  However, DELi decoding operates by first determining the library a read belongs then decoding the sequence
  based on that library call. For tool compounds to operate inside this logic, they become their own 'library'
  (of just one compound) for the purposes of decoding. Thus they need a tag to separate them from other libraries.
  This is what the ``library`` tag does when you define it for DELs. To make this logic clearer, you can also use
  the keyword ``tool`` instead of ``library`` here to define the library tag for tool compounds. Under the hood DELi
  will treat these the same way.
  .. note::
        The library tag for a tool compound must be unique from all other
        DEL libraries in the selection. If it is the same, you can try using
        :ref:`doped compounds <doped-compounds-sec>` instead.

- ``tool_compound_ref``: This section defines the static DNA sequence that DELi will use to check that a tool compound
  was called correctly. This is because unlike DELs there are no required decodable sections for tool compounds. So
  after the library call is made, the tool compound is the only choice (since it is a library of one). However, this
  can result is weird behavior if we just stop here. For example, if your tool compounds library tag happens to be a
  subset of another library's building block tags, and the read is missing some parts, DELi might miscall the library
  and suggest it is a tool compound when it is not. This isn't normally an issue, as this would break the downstream
  decoding and thus result in a failed read. However, with tool compounds there is no downstream decoding, so it
  would just accept this read as correct. To prevent this, DELi requires that tool compounds have a static
  ``tool_compound_ref`` section that it can use to verify the tool compound call was correct. If the tag has too many
  errors compared to the expected tag, DELi will reject the read.

- ``umi``: This section defines the UMI for the tool compound.
  UMIs are required for tool compounds to help identify PCR duplicates during decoding.
  This is :ref:`the same for DELs <barcode-sec-ref>`

The optional reserved sections are the same between tool compounds and DEL libraries.

Example JSON
------------
An example schema looks something like this:
.. code-block:: json
    {
        "compound_id": "ToolCompound_1",
        "smiles": "CCO",
        "barcode_schema": {
            "tool": {
                "tag":  "CCTTGGCACCCGAGA",
                "overhang": "CTA"
            },
            "compound_tag": {
                "tag": "ATTCCAATCGCTGA",
                "overhang": "ACG",
            "extra_tags": {
                "tag": "GGGNNNAGGGNNNATC"
            },
            "preumi": {
                "tag": "AATGCCAGTACG"
            },
            "umi": {
                "tag": "NNNNNNNNNNN"
            }
        }
    }

Tool Compound File FAQ
----------------------
- **Q:** My tool compounds have the same library tag but different barcodes. Does DELi support this?
    **A:** Not directly. DELi will only allow for one tool compound with a given library tag.
    However, you can get around this by defining multiple tool compound files, each with a unique library tag
    by include a region of the tag that varies between the two tool compounds.

    In the case that one of your tool compounds is a perfect subset of another, DELi offers no way to resolve this conflict.
    In general, this a bad idea for your DEL experiment as well, hiccups in ligation or synthesis could lead to miscalls
    between the two compounds.

- **Q:** My tool compound has the same library tag as a library in my selection. Can I still use it?
    **A:** Not directly. DELi requires that all library tags be unique across both libraries and tool compounds.
    However, you can try using :ref:`doped compounds <doped-compounds-sec>` instead, which allows you to define
    a compound within a library that behaves like a tool compound. You could also extend out the library tag to include
    more sequence to make it unique.

Deli Data Directory
-------------------
Tool compound definition files are stored in the :ref:`DELi Data Directory <deli-data-dir-ref>`,
under the ``tool_compounds/`` sub-directory.
