.. _defining_hamming-docs:

Defining Hamming Codes
======================

Most modern DEL libraries utilize a hamming encoded scheme to enable fixing of single base
pair misreads during sequencing (or during PCR) of the selection. This is most common in
the building block sections, helping to recover as many as 10% of reads to valid calls.

What is a Hamming Code?
-----------------------
Hamming codes allow for single "bit" errors to be detected and corrected.
This is done using a specific arrangement of bits, where some bits (parity) are
defined by the state of the other bits (data). The most common hamming code is
a hamming(7, 4) code, which uses 4 data bits and 3 parity bits. Using some
matrix multiplication, errors in the 7 bits can be detected and corrected, since
the code will assert that all possible combinations of 7 bits are a distance of 3
substitutions (swaps) apart. This means that if a single bit is flipped, the code will
be able to detect the original 7 bit code without any ambiguity.

In the case of DNA, binary bits are instead quaternary bits and the same
principle applies, allowing SNPs to be detected and corrected. However, INDELs
cannot be detected or corrected using hamming codes.

Hamming Codes VS Hamming Sets
-----------------------------
While its possible to use a single hamming matrix to build a code, thus requiring only
a single hamming matrix to decode, a valid hamming "set" is just a set of codes that
have a given minimum hamming distance from every other member of the set.
While you can generate the largest possible hamming set using a single matrix, it is
possible to randomly sample codes and check the hamming distance between them, resulting
in a valid hamming set, but one generated many possible hamming matrices. This second
"random" approach is far more common way of doing hamming encoding in DEL.

.. _deli-hamming-design-docs:

Building a Hamming Set of DNA Tags in DELi
------------------------------------------
DELi has a built in hamming encoding design module. It will build hamming codes
following a single hamming matrix approach. All you need to do it provide the
size of the barcode, and whether to use extra parity (enabling detection but not
correction of some double SNP errors):
::
    from deli import DELi.design.tag_factory import QuaternaryHammingTagFactory

    # create a tag factory for 12bp tags
    tag_factory = QuaternaryHammingTagFactory(desired_tag_size = 12)

    # extract all the tags
    tags = tag_factory.get_all_tags()

    # extract the Nth tag
    tag = tag_factory.get_tag(154)

If instead you desire a code that can handle a specific number of tags,
you can specify the number of tags you want instead and DELi will determine
how many base pairs your tag must be to support that:
::
    tag_factory = QuaternaryHammingTagFactory(desired_number_of_tags = 1000)
    tag_factory.total_size
    >>> 9
    tag_factory.max_tags
    >>> 1024

.. note::
    DELi does not support filtering a "bad" tags, so if your desired number of tags
    is close to the maximum number of tags, you may end up with a lot of tags that
    are not desirable, like "AAAAAAAAAAAA". It is best to create a tag factory much
    larger than you need so you can reject tags that are undesirable.

If extra parity is desired, you can set the ``extra_parity`` argument to True.
Note that extra parity will not increase the number of tags you can generate,
but will increase the size of the tags. For example, a 12bp tag with extra parity
become 13bp long. Extra parity is used only to assert additional checks to
detect more errors.
::
    tag_factory = QuaternaryHammingTagFactory(desired_tag_size = 12, extra_parity=True)
    tag_factory.total_size
    >>> 13

.. _deli-decoding-hamming-docs:

Decoding Hamming Codes with DELi
--------------------------------
If your DEL setting files include :ref:`barcode sections with hamming codes <barcode-sec-ref>`,
DELi will automatically determine how to decode the hamming codes for you.

If you are using a single matrix, all you need is the name of the hamming matrix to use.
If you used a DELi ``TagFactory`` to generate the tags, you can get the
matrix name using ``tag_factory.get_hamming_matrix_name()``. Otherwise you will need to know the
length of the code, and the number of parity bits used to generate the name, and make sure the matrix
is :ref:`defined correctly for that name in the DELi config file <deli-custom-hamming-docs>`.
Decoding this way is the most efficient from a speed and memory perspective.

If you do not know the matrix used, or if you used a random approach that used many matrices,
you can use the ``random`` method. This will generate a mapping of all possible Levenshtein
distance (Hamming but with INDEL support) 1 neighbors to the codes given in the building block
files. The number of 1-distance neighbors scale as the following: :math:`8x + 4` where x is the number of
nucleotides in the code. While :math:`O(n)`, they do all need to be loaded into memory (in a hash map).
If you have lots of libraries with large numbers of building blocks, it can start to take up a large
amount of memory. For example, using a 12bp codes with a decoding run using 1,000,000 across all
DELs in the collection used, that would be 100,000,000 tags to map: ~100 MB of memory (plus overhead).
However, since hash maps are fast in python, this is likly just as fast (maybe faster) than
the single matrix decoding approach.

Why not prefix trees
^^^^^^^^^^^^^^^^^^^^
Currently DELi only enable a hamming/levenshtein distance of 1 while decoding.
However, if a hamming set is not full, it is possible that we can determine the
correct code even if it is distance is greater than 1 from the true code.
This is done using a prefix tree to find the matches with the smallest levenshtein
to the read. If there is only one possibility with the smallest distance, that is
is the correct code (assuming the number of errors during sequencing is not exceedingly high).
In practice this doesn't have too much of an impact: the error rate of sequencers is about
1/100 bp, so the odds two errors occur within a 12bp code is pretty low. Right now, DELi
does not support this approach, in part because prefix trees in python are pretty slow.

.. _deli-custom-hamming-docs:

Designing Custom Hamming Codes/Matrix in DELi
---------------------------------------------
DELi support customization of the hamming matrices it uses for
specific size/parity combinations. DELi includes default
hamming matrices for code ranging from length 7-15 without extra parity, and 8-16
with it. This means the codes can handle a max hamming set size ranging from 256
to 4,194,304 members (before pruning undesirable DNA sequences).
You can also add custom codes if you used a different version for an existing
size parity pair, or if you need support for a larger code. You can read more
about hamming codes `here <https://en.wikipedia.org/wiki/Hamming_code>`_

.. note::
    It is not recommended to define custom codes.
    DELi defines default hamming codes for you, and creating
    custom codes is not always a straight forward process and can
    result in silent errors. Unless you really need it, it is best
    to avoid doing this.

Defining a hamming matrix in DELi
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Hamming matrices are defined in the :ref:`DELi config file <deli-config-docs>`
under the section ``deli.hamming``. In this section you must first delcare the
mapping of nucleotide (AGCT) to numbers (0123) that your code uses. This is
done using the ``nuc_2_int`` variable. The syntax is <nuc>:<int>,<nuc>:<int>,...
For example, if your code uses A=0, G=1, C=2, T=3, you would write:
::
    nuc_2_int: A:0,G:1,C:2,T:3

.. note::
    DELi assumes all specified hamming matrices use the same mapping of
    nucleotides to numbers. This is not a requirement if you are only
    using a single hamming matrix, but if you are using multiple hamming
    matrix and all use a different mapping, DELi does not support this (yet).

To then define custom hamming matrix you write a new section named ``deli.hamming.<NAME>``
where <NAME> is the name of your hamming matrix. DELi follows a syntax of <size>_<parity>,
for example ``deli.hamming.7_3`` means size of 7bp and 3bp are for parity.
In this section you need to define two variables:
``hamming_order`` and ``custom_order``. The hamming order is the "power of 2" bit order
associated with the code. This must follow a formate were parity bits are only present at
a power of 2 location (1,2,4,8,16...) and data bits are at all other locations.
parity bits are denoted with a ``p`` and data bits with a ``d`` and indexed from 1
for each type. For example the order for a 7_3 code is p1,p2,d1,p3,d2,d3,d4
If using extra parity, the order would be p0,p1,p2,d1,p3,d2,d3,d4, as p0 is the extra
parity bit. This follows the notation most common to hamming codes.
The custom order is how you might have remapped the bits in your code. For example,
if it sometime common place to put all the parity bits at the start of the code and data bits
at the end: p1,p2,p3,d1,d2,d3,d4. This will alter the hamming matrix needed to decode this code,
as data and parity bits are in different locations. DELi needs to know this order so you must
specify it. In the case your codes order of bits is the same as the hamming order, you just
like the same order again.

To tie it all together, here is an example of four different hamming codes defined in the config file:
::
    [deli.hamming]
    nuc_2_int = A:0,T:1,C:2,G:3

    [deli.hamming.8_4]
    hamming_order = p0,p1,p2,d3,p4,d5,d6,d7
    custom_order = p0,p1,p2,d3,p4,d5,d6,d7

    [deli.hamming.16_5]
    hamming_order = p0,p1,p2,d3,p4,d5,d6,d7,p8,d9,d10,d11,d12,d13,d14,d15
    custom_order = p0,p1,p2,d3,p4,d5,d6,d7,p8,d9,d10,d11,d12,d13,d14,d15

    [deli.hamming.7_3]
    hamming_order = p1,p2,d3,p4,d5,d6,d7
    custom_order = p1,p2,d3,p4,d5,d6,d7

    [deli.hamming.15_4]
    hamming_order = p1,p2,d3,p4,d5,d6,d7,p8,d9,d10,d11,d12,d13,d14,d15
    custom_order = p1,p2,d3,p4,d5,d6,d7,p8,d9,d10,d11,d12,d13,d14,d15
