########
Development FAQs
########

This section addresses some of the frequently asked questions (FAQs) related to the development decisions made in DELi.
This section is more historical to provide context on why certain choices were made during the development of DELi that
may not be immediately obvious.
It is not intended to be a guide for users or developers working with DELi, nor does it cover anything an average user
of DELi might need to know.

Q: What is a ``Library`` object in DELi?

A: A ``Library`` object in DELi represents a collection of compounds that share a single library ID. In the case that
these compounds have DNA tags, it also means they share the same DNA barcode schema. Outside of this, DELi does not
impose any additional requirements. Note that a `Library` object **does not** make any constraints on how many
compounds can be in the library. This present some confusion between the concept of a "library" in DELi and
the concept of a "library" in chemistry, where a library is often thought of as more than one compound.
This can result is some oddities in implementation, especially during the demultiplexing stage of decoding.
In order to demultiplex reads, DELi first separates reads based on library tags, so ever read in the system must
have a unique tag that identifies the library it belongs to. Yet there times when single tagged compounds are used
outside of the normal DELs. These still need to be separated and and decoded, but they don't generally have a "library"
tag, since to the chemists they are not part of a library.

To address this, DELi masks a bit of behavior away from the users to try and avoid confusion. To users, single compounds
are implemented as ``Compound`` objects, even if they have tags, as you would expect. Yet when DELi looks at the DNA
barcode, it will map the part of the tag defining that compounds ID to a ``LibraryBarcodeSection``. In this way, DELi
is treating these compounds as if they are in a library, even if the user does not think of them that way. This get more
concrete when DELi is demultiplexing reads. Since every read must belong to a library, DELi converts any added single
compounds to ``Library`` objects behind the scenes. This makes it very easy to handle all reads in a uniform way,
even if some of them are from single compounds and are not ``Library`` objects to start.
