==================
Defining Reactions
==================

This section covers how to define reaction files that DELi can load, not how to define a reaction scheme
for a DEL. For that information see the :ref:`defining DEL reaction scheme docs <reaction-sec-ref>`.

File Format
===========
The reaction file format is a simple single line plain text file, with a ``.rxn`` file extension containing
a single reaction SMARTS.
DELi will read in the entire file, strip the whitespace at the front and end and try to load the reaction.
If it is not a valid SMARTS, an error will be raised. DELi does not support multi line reaction SMARTS.
An example reaction file might look like:

::
    [c:1][OH:2].[N:3]>>[c:1][O:2][N:3]
