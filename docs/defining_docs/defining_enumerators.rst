.. _define-enum-docs::

====================
Defining Enumerators
====================

Enumerator definitions are nearly identical to DEL library definitions, as both require the
building blocks and reaction scheme to be defined. The main difference is that enumerators do not
require barcode information, as they are not tied to DNA tags. Those you can define a enumerator
using a JSON file much like a DEL library file, but without the barcode schema section.

You can also load a enumerator directly from a DEL library definition file, or by by accessing the
``enumerator`` attribute of the loaded ``Library`` object.

This is only useful in a context where you don't know the tag information but do know the SMILES and
reaction scheme. This could be helpful when designing a DEL library, where you want to see the enumerated
compounds before assigning barcodes.
