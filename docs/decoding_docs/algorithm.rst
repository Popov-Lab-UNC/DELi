.. _decoding-algorithm:

==================
Decoding Algorithm
==================

DELi uses a two step algorithm to decode raw DEL selection data into compound enrichment data.

Library Demultiplexing
----------------------
First it "demultiplexes" the DNA reads based on which DEL the read originated from.
Most DEL selection actual screen a collection of DELs, not just a single library.
Since DELi supports decoding of "heterogeneous DELs" (DELs that might have very different DNA barcodes designs),
without knowledge of the library the read originated from, it is not possible to know how the barcode is designed.
In order to decode the read, the correct barcode design for that read is needed (so we know how to map each section).
The library demultiplexing step solves this problem as it can "call" the library for each read.

How this is done is by using a modified semi global sequence alignment algorithm. Specifically, DELi uses
`cutadpat <https://cutadapt.readthedocs.io/en/stable/>`_ to handle the demultiplexing under the hood. If you
are interested in more details on how demultiplexing with semi global alignment works, you should read the
`cutadapt documentation <https://cutadapt.readthedocs.io/en/stable/guide.html#demultiplexing>`_.

Building Block Calling
----------------------
After the demultiplexing step, each read is assigned to a specific library, thus a specific barcode design.
This design is then converted into a DNA sequence that can serve as a reference for a sequence alignment.
Some nucleotides in reference are unknown, and represented as "N". These are what DELi will try to match
to a query to figure out what the unknown nucleotides are, and then figure out what chemical that corresponds to.
This is also done using a semi global alignment, but this time with padding on the ends to avoid trimming the front/back of the read
and with a custom scoring matrix that does penalize for mismatches with N. Assuming there are only a handful of alignments with low scores,
DELi will loop through them, each time trying to decode the building block barcodes based on the alignment. If it fails (because that
barcode doesn't exist) it will try the next alignment and continue until it finds a valid decoding or runs out of alignments.

To do the conversion of the mapped nucleotides to building blocks, DELi uses the hashmaps, where the keys are the nucleotides and the values are the
building block they correspond to. While the default is to just use a direct, zero error tolerance mapping, you can also enable several
:ref:`error correction <error-correction-docs>` methods to allow for some flexibility in matching.

UMI Degeneration
----------------
After the read has been called as a specific DEL compound (library called and all building blocks decoded), the last step is to handle
UMI Degeneration (for qPCR). This is where reads are group as the same if they share the same UMI (their count is degenerated so to speak).
In DELi we call these "UMI Corrected Counts" and they are the final output of the decoding process.
This is optional in theory (and DELi supports that) but in practice it is required if you want "good" DEL results, as UMI corrected counts
remove the large amount of noise that can be introduced by PCR amplification.
DELi supports a few methods of UMI Degeneration.

Exact UMI Degeneration
^^^^^^^^^^^^^^^^^^^^^^
This is the simplest method, where all reads with the same UMI are grouped together and their counts are summed.
There is no room for error in the UMI. Given that even high end sequencers have error rates of about 1/100 base pairs,
This is going to introduce some noise, as AACCTTGG and AACCTTG*A* will be counted as the **different** UMIs, even though,
statistically speaking, the odds are far more likely one is a misread of the other. [*]_

.. [*]_
    To demonstrate this, say we have 100bp DNA barcode with a 10 bp UMI. on average 1 out of every 10 reads (10%) will have a single error
    in the UMI, giving it a hamming distance on 1. Now, with 10bp UMIs, there are 1,048,576 possible UMIs, but for any given UMI, there
    are only 30 other UMIs with a hamming distance of 1, or a 0.0029% chance of a read being misread as one of those. If we assume that the
    chance there are far more possible UMIs that those actually in the DEL (Since DELs are synthesized at such low concentrations) it is 4 fold
    more likely that two UMIs with hamming distance of 1 are just one single UMI (one has a misread) than it is to be two unique UMIs.
    In practice these odds are even lower, as UMIs are often longer than 10bp.

Clustered UMI Degeneration
^^^^^^^^^^^^^^^^^^^^^^^^^^
To compensate for misreads and PCR errors during sequencing, DELi implements a "clustered" UMI degeneration method.
It is similar the to clustering method outlined in Smith et al. 2017 [1]_. The major difference is that DELi uses
a greedy approach, where the centroid of the cluster just the first UMI observed that didn't fit into an existing cluster.
In Smith et al, the centroids are picked as the most common UMIs not already in a cluster.

References
----------
.. [1] Smith, T., Heger, A., & Sudbery, I. (2017). UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Research, 27(3), 491-499. doi:10.1101/gr.209601.116
