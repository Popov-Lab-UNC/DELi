.. _decoding-algorithm:

==================
Decoding Algorithm
==================

DELi uses the following process to decode raw DEL selection data into compound enrichment data.

1. DELi will attempt to align the raw read to the various possible library barcode designs.
2. Those alignments to decode the single library the read is from.
3. A refined alignment for the called library is performed to map barcodes regions to the read.
4. Variable barcode regions are decoded to determine the final compound.
5. Compounds are grouped by unique UMIs to produce final enrichment data for each compound.

It any read fails at any of these steps (due to sequencing or PCR errors, or because there are multiple possible answers),
that read is discarded as not decodable. For a highly quality sequencing,
DELi can get anywhere from 80-90% of the reads decoded successfully.

Sometimes, different approaches blur these steps together. For example you could align each read to a all
possible libraries at the same time, which would both generated a refined align and decode the library in one step.
To break things down a bit better, DELi separates these steps into three distinct modules:

1. Library demultiplexing
2. Barcode decoding
3. UMI degeneration

Library Demultiplexing
----------------------
Library Demultiplexing is responsible to determining the library a read originated from as well as
generating a set of possible :ref:`alignments <alignments>` of that read to the
:ref:`library's barcode design <barcode-sec-ref>`.
This is a crucial part of decoding, as in a typical DEL selection experiment screens a collection of DELs,
not just a single library. Further, if there is not alignment between the libraries barcode design and the
read, it will not be possible to decode the variable barcode regions to determine that compounds ID.

At first glace this is a bit tricky, since to call the library, you need to know how the library barcode is,
but to know that you need an alignment, but to do that you need to know the library. DELi provides three solutions
to this problem.

.. _full-barcode-demulti::

Full barcode demultiplexing
^^^^^^^^^^^^^^^^^^^^^^^^^^^
A trivial solution is to just align the read to the full barcode design of each library, and pick the best alignment.
This solves both problems at once. DELi support this method with the `full` demultiplexing mode.
DELi uses Biopython's `PairwiseAligner <https://biopython.org/docs/dev/Tutorial/chapter_pairwise.html>`_ for a semi-global
alignment adjusted for `N` wildcards to generate these full barcode alignments.

There are drawbacks to this approach. Primarily, it is computationally expensive.
In this frame work, if you have 50 libraries, you need to perform 50 alignments per read. Early stopping on a
perfect alignment can help, but not by much. Things are even worse if you need to search the reverse complement.
This will double the number of alignments that are needed.
However, there are benefits as well. This method is the most flexible to read quality, and generally
decodes to the most reads. If you have nearly infinite "free" compute, this isn't a bad option if you want
to squeeze out every last read. It is also a good choice if your reads have lots of errors in them, either due
to poor PCR or poor sequencing quality.

Single primer demultiplexing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A more efficient method is to just align to a small, known (static) section of the barcode design that is shared among
all or many libraries. Then you can located and align to this static section and, determine the library tags,
and lookup which library matches that (see :ref:`barcode decoding <TODO>` for more info on how that lookup occurs). Then
a final alignment between the whole read and the full barcode design of the called library can be performed
to locate the other variable regions.

This method is called `single` multiplexing mode in DELi. The reason it is refered to as primer here is because
for most DELs, if you have a static region shared amoung many DELs, it is usually a primer (used for qPCR).
However, any static section will suffice.

Currently, DELi will not let you choose which static region to use. Instead it picks the static region that is
present in the most libraries, thus can eliminate the most work. Overhangs are not considered static for this
purpose.

.. note::
    For DELi to know that two libraries share the same static region, the static region must be defined
    exactly the same in both library designs. This means the DNA tag *and* the section name must be the same
    in both library designs.

This method is far more efficient than the full barcode demultiplexing method, and in many cases can demultiplex
in 1-2 alignments. However, if all your DELs share no static regions, this will not be any faster *or* better than
just using the full barcode demultiplexing approach.

Flanking primer demultiplexing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
As mentioned above, DELs often share static regions that are primers used for qPCR. This only works if there is a
primer on both sides "flanking" the whole barcode. Therefore you could align to both primers. This can help eliminate
incorrect alignment from a single primer. This is because the single primer lacks the ability to determine any info
about the distance between sections. However, with two primers flanking either side of the library section, you can
make sure the distance between the primers roughly matches the expected distance for that library, which lets you
rule out some libraries right away. It can also help with :ref:`static section alignment <static-section-alignment>`
as it will be able to better detect and address INDELs that could shift the tag.

Like with the single primer, DELi will pick the flanking static regions for you, and they don't have to be the
primers, any two static regions that appear in many libraries on either side of the variable regions will work.

Generally, this is not much more expensive than the single primer method, and can be more accurate in some cases.
It is hard to say when it will be better and might be context dependent, so if you have the compute, it is worth
trying both single and flanking to see which works best for your data.

Core demultiplexing algorithms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
DELi's goal is not to reinvent the wheel, so it takes advantage of existing, well tested tools to perform
alignment/search part of the demultiplexing. Currently, there is support for three different open source tools:

- `Biopython <https://biopython.org/docs/dev/Tutorial/chapter_pairwise.html>`_ PairwiseAligner
- `Cutadapt <https://cutadapt.readthedocs.io/en/stable/>`_ Adapter search
- `Regex <https://github.com/mrabarnett/mrab-regex>`_ fuzzy search

.. note::
    `regex` is not the same as `re` from the python standard library. It is a third party library that
    supports more complex regular expression syntax (but is still computable with python `re`).

All of these tools have different strengths and weaknesses. Biopython is very flexible, but slow.
Cutadapt is very fast, but less flexible. Regex is middle of the road, giving good speed and flexibility.
By default, DELi uses Regex for demultiplexing.

.. _biopython::

Biopython
~~~~~~~~~
This will use Biopython's PairwiseAligner to perform semi-global alignments between the read
and the static region(s). Right now it uses a custom substitution matrix that treats `N` as a wildcard
(0 score for match with anything). Gaps are penalized normally with opening as -2 and continuing as -1.
Matches are scored as +1.

Cutadapt
~~~~~~~~
This will use Cutadapt's alignment algorithm to find the static region(s) in the read.
Cutadapt is a fantastically optimized tool for locating adapters in sequencing reads. However,
adapters are meant to be at the ends of reads, while static regions could be anywhere in the read.
Usually, primers are at the ends of reads so things work out okay, but if you have static regions
in the middle of the read, Cutadapt might not be able to find them as well as other methods.
This is because the alignment algorithm prioritizes finding matches at the ends of reads, and will
pick them over matches in the middle of reads *even if* there match in the middel is better (less errors).
This is not a bug but a feature of cutadapt, as it is designed for adapter trimming, not general sequence
searching. However, it can cause strange behavior in DELi. Only use this approach if you are sure your static
regions are at the ends of your reads (like primers). DELi will warn you if it detects static regions
in the middle of reads when using Cutadapt.

.. warning::
    If you are considering using cutadapt over regex for demultiplexing, you should read the above
    paragraph, there are important caveats to be aware of.

Regex
~~~~~
This will conduct a fuzzy regular expression search for the static region(s) in the read. Unlike a normal
regular expression which cannot handle insertion, deletion or substitution errors appearing arbitrarily in a
specific section of the query, the `regex` library supports fuzzy searching. DELi will build out regex expressions
that look like ``(?:(?P<primer1>AGCTAGCT)){e<=1}.{{10,20}}(?:(?P<primer2>CGTACGTA)){e<=1}``. While far more confusing
looking than a sequence alignment, this approach is faster than a full semi-global alignment, and provide a solid
foundation for :ref:`static section alignment <static-section-alignment>`

.. _alignments::

Alignments
^^^^^^^^^^
After using the demultiplexing step to call the library, DELi needs to generate alignments between the read and the
full barcode design of the called library. In DELi, alignments aren't full blown global or local alignments with scores,
but rather they are simply a mapping between the barcode sections and the DNA sequence of the read that corresponds
to that section (based on the alignment algorithm).
This is needed to map the variable barcode regions to the read, so that
the building block barcodes can be decoded. DELi implements two approaches to generate these alignments:
semi-global alignments and static section alignments. There are some pretty big differences between these two approaches.

.. note::
    In both approaches, DELi will return several possible alignment depending on the settings.

Semi-global alignments
~~~~~~~~~~~~~~~~~~~~~~
This is the more traditional approach, where DELi uses :ref:`Biopython's PairwiseAligner <biopython>` to perform
a semi-global alignment between the read and the full barcode design of the called library.
This will generate an alignment that lines up the read and the barcode design as best as possible taking into
account the overhangs and other static regions to better locate possible INDELs.

This approach will completely ignore whatever initial alignment was generated to help call the library (unless
a :ref:`full alignment was done already to call the library <full-barcode-demulti>`) and build a new one from
scratch. It can be triggered by telling DELi `realign = True`. It is more expensive than the other approach, but
can recover an addition 2-5% of reads in some cases, especially if there are a lot of INDELs in your reads.

.. _static-section-alignment::

Static section alignments
~~~~~~~~~~~~~~~~~~~~~~~~~
This approach takes advantage of the initial alignment(s) generated during the demultiplexing step to call the
library. Since those alignments already map the static regions of the barcode design to the read, DELi can use
those as "anchor points" to build the full alignment. This is much faster than doing a full semi-global alignment,
as it only needs to figure out how to map the variable regions in between the static regions.

This is done by looking at the distance between the static regions and the other regions in the barcode design,
and then uses the known location of the static regions used to find the library to guess where the other regions should be.
For a single known static region this will just seek out the variable regions based on their expected distance from the
known static region. For example, if the first building block section is known to be 15 bp ahead of the static region
and is 11 bp long, it will grab the DNA sequence 15-25 bp ahead of the known static region alignment and assign that as
the barcode sequence. However, if there are INDELS between these two section this could result in a misalignment.

With two known static regions on either side of the section in question, DELi can do a bit better.
It can use the observed vs expected distance between the two known static regions to adjust the
alignment of the section that fall between to account for INDELs. For example I could have the read
and barcode schema of::

AGCTAGCTGGGGAGCTAGACTAGCTAGCT
AGCTAGCTJJJJNNNNNNNKKAGCTAGCT

Where the N is the region we want and J and K are some other regions we don't care about. If three was an INDEL such that::

AGCTAGCTGG-GAGCTAGACTAGCTAGCT
AGCTAGCTJJJJNNNNNNNKKAGCTAGCT

Then by looking at the distance between the two static regions (AGCTAGCT) DELi can see that there is a 1 bp deletion
in the read between them. Therefore it can adjust the variable regions accordingly, so it would extract
``NNNNK`` instead of ``NNNK``. When calling the library barcode, `NNNNK` is only 1 levenshtein distance away from
the expected barcode making it easy to fix, but ``NNNK`` is 2 distances away, making it much harder to correct.

In cases where there might be more than one INDEL between two known static regions, (say it the above example had 2 more
deletions and created ``JJNNNNK`` as the variable region) DELi will walk through all sequential sets of nucleotides
that are of length equal and minus one to the expected length of the variable region to try and find a match.
This is still far faster than a full semi-global alignment, but also still a bit more error prone.

Barcode Decoding
----------------
Once the read has been assigned to a specific library and an alignment(s) between the read and the library
are created, then the sections can be decoded. This is pretty straight forward. DELi will loop through each alignment
and try to decode the building block barcodes based on the mapped nucleotides for each variable region.
If it fails (because that barcode doesn't exist) it will try the next alignment and continue until it finds a valid
decoding or runs out of alignments. Hashmap based :ref:`barcode calling <>` is used to allow for possible errors in the
barcode regions when decoding while keep runtime as low as possible (O(1) in this case).

Once all the variable sections are called, DELi will then locate the UMI by both using the alignment to find
where the UMI is located in the read *and* looking taking the closest aligned variable section to the UMI and
searching *n* bp away to find the UMI sequence. Both these (if different) are save as possible UMIs for the read,
with the *n* bp away being the primary UMI call. Future updates will build upon this to do more advanced UMI
correction.

This then produced a full decoded read, mapping it to a unique compound ID and associating it with a UMI.

Wiggle
^^^^^^
DELi supports a "wiggle" mode when calling building blocks. This is a simple method to help recover
reads that have insertions or deletions (INDELs) in the building block region of the barcode.
When wiggle mode is enabled, after the initial alignment of the read to the barcode schema, DELi will
extend the aligned region by 1 bp on each side, and then scan all possible chunks of the expected
building block tag length, 1 smaller and 1 larger than expected length (in that order).
If a match is found in any of these chunks, the building block is called.

Turning on wiggle is recommended when not using a semi-global alignment during demultiplexing,
as it can help recover reads lost to INDELs in regions outside of the sections of interest.


UMI Degeneration
----------------
For any given compound ID, DELi will group all reads that share the same UMI together and count them as one.
This is to correct for PCR amplification bias, where some molecules get amplified more than others.
Right now DELi does not implement any advanced algorithm to account for errors in the UMI (like UMI clustering or
UMI graphs [1]_). This is planned for future releases.



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
It is similar the to clustering method outlined in . The major difference is that DELi uses
a greedy approach, where the centroid of the cluster just the first UMI observed that didn't fit into an existing cluster.
In Smith et al, the centroids are picked as the most common UMIs not already in a cluster.

References
----------
.. [1] Smith, T., Heger, A., & Sudbery, I. (2017). UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Research, 27(3), 491-499. doi:10.1101/gr.209601.116
