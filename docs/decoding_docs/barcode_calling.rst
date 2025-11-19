======================
Barcode  Calling
======================
DELi uses hashmaps (python dictionaries) to map DNA barcode sequences to their corresponding objects.
There are generics in DELi, so it can map tags to any type of object, but in practice we only care
about DELs and building blocks.

There are several ways to build these hashmaps. The simplest is to just map the exact sequence
to the object. This works well when there are no errors in the reads. However, in practice there are often errors
in the reads, especially in the building block region of the barcode. To mitigate this, DELi supports
several error correction methods that can be used when building the hashmaps for building block calling.

.. _error-correction-docs:

Error Correction
================
Error correction works by building a "spell checker"
(see Peter Norvig's `article on spell checking <http://norvig.com/spell-correct.html>`_).
This when decoding is first initiated, a one time preprocessing step is done to build a mapping of
all possible erroneous sequences to their corresponding object of origin.
There are many ways to determine what erroneous sequences are possible from a valid sequence,
and how to handle things like ambiguous matches
(where more than one valid sequence is equally close to the erroneous sequence).
DELi implements several type of methods outlined below.

Types of Error Correction
-------------------------
Ideally, the users know how this was done for the library. For example, if all DNA tags of a given building block set
have a hamming distance of 3 from all others (called a hamming-3 set) then we can nearly guarantee that any invalid read with
a hamming distance of 1 from a valid read is most likely the result of an error during sequencing/PCR, thus should be
corrected to the original sequence.

This is just one example of possible error correction methods. DELi implements/supports several methods outlined below.

Hashmap Error Correction
^^^^^^^^^^^^^^^^^^^^^^^^
Hashmap error correction involved an initial brute force mapping of all possible error containing sequence to their
corresponding sequence of origin (the correct building block sequence). This works by using somme discrete distance metric
to determine all possible sequences that are within a certain distance of the valid sequence. As you can imagine, for large
distances this can lead to a large number of possible sequences, and thus a large amount of memory and run time to build.
However, the building only is required once per building block set, so if you are doing a large decoding run with many reads,
it can be worth the time to build the hashmap, since once built error correction has a complexity of O(1) for each read.

By default, the implementation of this in DELi assert there are no "collisions" in the hashmap. This means that for all pairs
of valid sequences, the set of sequences within their distance thresholds have no overlap. This is because if there is a collision,
it become unclear/unconfident which sequence the read should be corrected to.
Therefore, it is best to use the hashmap error correction method when you know the build block set tags are a given distance away
from each other (meaning there are no collisions).

Before decoding starts, DELi will check this and fail if it detects a collisions.
In case a user doesn't care about collisions, or always wants the best match,
this behavior can be modified by using the `asymmetrical mode <#asymmetrical-mode>`_ outlined below.

DELi currently supports two types of distance functions for building the hashmaps for error correction:
`Hamming distance <https://en.wikipedia.org/wiki/Hamming_distance>`_ and
`Levenshtein distance <https://en.wikipedia.org/wiki/Levenshtein_distance>`_.

Hamming Hashmap
~~~~~~~~~~~~~~~
The hamming hashmap will generate all sequences within a hamming distance of the valid sequence from the building block set.
Collision can be guaranteed to not happen if the building block set is designed such that all valid sequences are at least
a hamming distance of 3 apart for hashmap cutoff distance 1, or 5 for hashmap cutoff distance 2, or :math:`2d + 1` for a cutoff
distance of :math:`d`.

To give you an idea of scale, if my building block set 12 bp long sequences, the equation to determine how many sequences
have a hamming distance of :math:`d` from the original is :math:`\sum_{n=1}^{d}3^n\binom{12}{n}`, where :math:`\binom{12}{n}`
is the binomial coefficient (12 choose :math:`n`). For :math:`d=1` this is only 36, for :math:`d=2`
it is 630, and for :math:`d=3` it is 6,570. So if we have a BB set of 5000 tags, that is ~30 million sequences in the
hashmap for a cutoff distance of 3.

Levenshtein Hashmap
~~~~~~~~~~~~~~~~~~~
The Levenshtein hashmap will generate all sequences within a levenshtein distance of the valid sequence from the
building block set. Just like with the hamming hashmap, collisions can be guaranteed to not happen if the building block set
is designed such that all valid sequences from teh BB set have a Levenshtein distance of :math:`2d + 1` for a cutoff
distance of :math:`d`.

The Levenshtein hashmap is probably the best option to use for error correction, as it allows for insertions and deletions (INDELS).
However, the can take longer to build and use more memory than the hamming hashmap, especially for larger distances.
Unlike the hamming distance, the levenshtein distance allows for insertions and deletions, which makes a closed form equation
hard to write. However, taking the length to be 12 bp again, with a distance of 1 it is about 100 sequences, with :math:`d=2` it
is about 5,000 sequences, and with :math:`d=3` it is about 150,593 sequences. This is much large than the hamming maps,
so be careful when building these. On an average laptop, building a Levenshtein hashmap with a distance of 2 takes about 30-40 seconds.
for a set of 5,000 tags. This scales linearly with the number of tags, so for 10,000 tags it would take about 70 seconds.

.. note::
    Given the size of the barcodes, the error rate of sequencing/PCR, and the size of the building block tags, it is unlikely to
    observe more than 2 errors in a specific building block tag, so in practice a max distance of 2 is often sufficient.

.. _asymmetrical-mode:
Asymmetrical Mode
~~~~~~~~~~~~~~~~~
While the default behavior of the hashmap error correction is to fail if it detects a collision, this is not always desired
given the nature of a DEL. Designing a set of 5,000 DNA barcodes that have a Levenshtein distance of 5 or more is a
difficult task, so more often than not the cap is set to distance of 3, meaning out distance cutoff should be 1 for the hashmap.
However, the space of possible barcodes is 16,777,216, so there are bound to situations where even if the distance of a query into
set of valid sequences has a distance of 2 for the nearest valid sequence, it is the only valid sequence within a distance of 2. We
can be more confident that this is the correct sequence to go back to. But this will not be true all the time, only sometimes.

Normally, we might use a Trie or a BK-tree to handle this, but these are still somewhat slow (especially in pure python).
Instead DELi implements an "asymmetrical mode" for the hashmap error correction, where when building out the hashmap with an upper
distance cutoff, it will allow for collisions, mapping to the valid sequence that has the smallest distance. This will still fail for
ambiguous matches (meaning more than one valid sequence has the same minimum distance to the query), but does enable some more flexibility.

So, even if you built a BB set with tags of a Levenshtein distance of 3 from each other, you could using a hashmap with a distance cutoff
of 2 and the asymmetrical mode enabled. This could help correct even more errors, at the cost of a longer initialization time to built
the hashmap.

Quaternary Hamming Decoder
^^^^^^^^^^^^^^^^^^^^^^^^^^
Hamming distance was invented in tandem with the Hamming code, a way to enable single bit error correction in binary data.
Hamming codes require specific sequences of the bits, some holding data some holding parity bits (which are set based on the
data bits and used to locate errors).
When you generate a set of sequences at random and check if their hamming distance to all others in the set is
at least 3, you cannot always guarantee that the same hamming code (order of bits) is the same for all sequences.

However, you can just pick a hamming code and use it to generate a set of sequences that are guaranteed to have a hamming
distance of 3 from each other *while* having the same order of bits. If this was done, you can then use a hamming decoder
to decode the reads and correct single SNP errors in the building block region of the DEL DNA barcode.

DELi implements a hamming decoder generalized from binary to quaternary space (the space of DNA sequences). If you know your
library was designed this way, you can use the hamming decoder. You just need to make sure DELi is configured to use the
your Hamming code (see :ref:`custom hamming docs <deli-custom-hamming-docs>`).

.. note::
    In practice, this is slower than the hashmap error correction, and far less likely to be needed, as
    DELi has on of the few (if not the only) implementations of a quaternary hamming encoder. Odds are,
    your tags were created at random, not with a hamming code. Also it is limited to a single SNP error correction,
    which is not always sufficient for DELs with longer building block regions.

.. _error-correction-methods-config:

Specifying Error Correction Methods
-----------------------------------
When you :ref:`define a DELs barcode schema <barcode-sec-ref>`, you can specify the error correction method to use for
the building block DNA tag region. This is done by setting the ``error_correction`` for that section of the schema.
The format for the value of this takes the form "<method>:<arguments>", where <method> is a reserved keyword for the
type of error correction to use, and <arguments> are some comma seperated values to configure of it runs.

The currently supported method keys are:
- ``hamming_dist``: Use a hashmap error corrector with hamming distance
- ``levenshtein_dist``: Use a hashmap error corrector with levenshtein distance
- ``hamming_code``: Use a quaternary hamming decoder

.. _hamming-dist-err-correction-format:

``hamming_dist``
^^^^^^^^^^^^^^^^
hamming_dist takes a single argument, the distance cutoff to use for the hashmap.
You can also add ``asymmetrical`` as an optional argument to enable :ref:`asymmetrical building <asymmetrical-mode>`.
For example, to use a hamming distance of 2 with asymmetrical mode enabled, you would set the value to
``hamming_dist:2,asymmetrical``. If you wanted to use a hamming distance of 1 without asymmetrical mode,
you would set the value to ``hamming_dist:1``.

``levenshtein_dist``
^^^^^^^^^^^^^^^^^^^^
levenshtein_dist takes the same arguments as :ref:`hamming_dist <hamming-dist-err-correction-format>`.
For example, to use a levenshtein distance of 2 with asymmetrical mode enabled, you would set the value to
``levenshtein_dist:2,asymmetrical``. If you wanted to use a levenshtein distance of 1 without asymmetrical mode,
you would set the value to ``levenshtein_dist:1``.

``hamming_code``
^^^^^^^^^^^^^^^^
.. note::
    This mode has been deprecated and will be removed in a future release. It is recommended to use the
    :ref:`hamming distance hashmap error correction <hamming-dist-err-correction-format>` instead.

hamming_code takes a single argument, the name of the hamming code to use. This much match a name of a hamming code
listed in your :ref:`deli config <deli-config-hamming-section>`.

For example, if you have a hamming code called "my_hamming_code" in your config file, you would set the value to
``hamming_code:my_hamming_code``.

.. warning::
    This makes ":" and "," reserved characters for naming hamming codes. DELi will raise an error if you try to use
    these characters in the name of a hamming code.
