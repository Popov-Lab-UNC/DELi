.. _decoding-algorithm:

==================
Decoding Algorithm
==================

Broadly speaking, DELi decodes raw DNA read into compound IDs using a simple 2 step system:

- Barcode Alignment: Determine which spans of the read belong to which sections of a DELs DNA barcode
- Section Calling: Use the alignments to call the decodable sections of the barcode

It would be nice if these steps were completely separate and flowed linearly, but in DELi these steps are intertwined.
For example, to align the barcode you need to know which library to align the read too. To do that you need to know
which library the read is from (so you can get to correct alignment). To do that you need to align the read and
call the library section, which we already established cannot be done until we know the library. This doesn't make
decoding impossible, instead it means a partial alignment to locate the library section must be done first, then a
the library section called, then the alignment can be finished and the rest of the sections called.

Most of this confusing interweaving happens during the first stages of decoding, where the reads library is determined.
To isolated this, DELi refers to this stage as "library demultiplexing". It is during this stage the the reads library
is called and an alignment between that libraries barcode and the read is generated.

After reads are converted to compound IDs, DELi can also do UMI correction ("degeneration" as it is called in DELi)
to create a final count of unique compounds observed for a given compound ID.

.. _library-demultiplexing:

Library Demultiplexing
----------------------
DELi implements four different approaches for library demultiplexing:

- Full barcode demultiplexing
- Library section demultiplexing
- Single section demultiplexing
- Flanking section demultiplexing

The last three approaches are roughly grouped together as :ref:`"static alignment" <static-section-alignment>` approaches, as they all rely on
finding static (known) sections of the barcode to help align the read to the barcode design rather than doing a full
sequence alignment.

Each have benefits and drawbacks under different situations.
The goal of so many methods is to give flexibility between speed and robustness.

.. _full-barcode-demulti:

Full barcode demultiplexing
^^^^^^^^^^^^^^^^^^^^^^^^^^^
An obvious way to determine which library is the read is from is to just try a full alignment to each possible
library barcode schema. The one with the best alignment score is (likely) the correct library. In ``full`` mode
DELi uses a or a semi-global alignment with a custom scoring matrix that does not penalize for the 'N' wildcards in
the barcode schemas during alignment. It is implemented using Biopython's
`PairwiseAligner <https://biopython.org/docs/dev/Tutorial/chapter_pairwise.html>`_.

The benefits of this approach are that the semi-global alignment is very flexible to read quality and it can be "reused"
when calling the other decodable sections, since the alignment to the full library barcode schema was already done to
find the library. However, this approach is "brute". If your DEL selection used 70 libraries, then for every read you
will need to run 70 semi-global alignments, one per library. It also scales directly with the size of reads and number
of libraries. This makes the full barcode demultiplexing the slowest of the four approaches


.. _static-section-alignment::

Static section alignments
^^^^^^^^^^^^^^^^^^^^^^^^^
Rather than spending the resources to do full dynamic sequence alignments, it can be alot faster to just look
for smaller static sections of the barcode to help align the read. Once the location of a static section is known in the
read, DELi can use that as reference to build the full alignment. This is much faster than doing a full
semi-global alignment, as it only needs to figure out how to map the decodable sections in between the static sections.
Luckly this can be done once and stored for the full decode run, allowing for nearly instant alignment generation after
the static section(s) are found.

This is done by looking at the distance between the static sections and the other sections in the barcode design,
and then uses the known location of the static sections used to find the library to guess where the other sections should be.
For a single known static section this will just seek out the variable sections based on their expected distance from the
known static section. For example, if the first building block section is known to be 15 bp ahead of the static section
and is 11 bp long, it will grab the DNA sequence 15-25 bp ahead of the known static section alignment and assign that as
the barcode sequence. However, if there are INDELS between these two section this could result in a misalignment.

With two known static sections on either side of the section in question, DELi can do a bit better.
It can use the observed vs expected distance between the two known static sections to adjust the
alignment of the section that fall between to account for INDELs. For example I could have the read
and barcode schema of::

AGCTAGCTGGGGAGCTAGACTAGCTAGCT
AGCTAGCTJJJJNNNNNNNKKAGCTAGCT

Where the N is the section we want and J and K are some other sections we don't care about. If three was an INDEL such that::

AGCTAGCTGG-GAGCTAGACTAGCTAGCT
AGCTAGCTJJJJNNNNNNNKKAGCTAGCT

Then by looking at the distance between the two static sections (AGCTAGCT) DELi can see that there is a 1 bp deletion
in the read between them. Therefore it can adjust the variable sections accordingly, so it would extract
``NNNNK`` instead of ``NNNK``.

DELi also supports a "wiggle" mode for static alignments, where the alignment span is extended by 1 bp to the right.
This dramatically increases the ability to recover reads with INDELs during :ref:`barcode calling <barcode-calling>`.

Under the hood, DELi uses the ``StaticSectionAligner`` class to handle this logic.

.. _library-section-demulti:

Library section demultiplexing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Using static section alignment, you could instead just
look for the each unique library section rather than align to each DELs barcode.
Since all libraries must have unique library section barcodes this will still allow for determining the library.
This is called ``library`` mode and is about ten times faster that using
:ref:`full barcode demultiplexing <full-barcode-demulti>`.

There are still downsides with ``library`` mode that is not easy to overcome. The library barcodes being so short is
one. Searching for barcodes that are as small as 8-12bp has an increase risk of mismatches cause by SNPs or random chance.
DELi has some extra bells and whistles to try and compensate (using the known size of the barcode to help it pick
the best match for example), but it is still a problem. There is also
the issue of scaling with the number of libraries. While the searches are faster, it would still be better if we
could get away with doing less searches rather than one per library.

.. _single-section-demulti:

Single section demultiplexing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Single section demultiplexing, called ``single`` in DELi, addressed both of the issues presented by the library section
demultiplexing method. Instead of searching for the library barcode a static section is searched for instead to
generate the static alignment. The main benefit here is that unlike
library sections, which must be unique between each library (thus requiring the a single search per library) static
sections can be shared among many library barcode designs. This is often the case, as DELs made by a single
organization or vendor will often have the same primer sections that are static across all the DELs. If we used 50 DELs
in our selection, but they all have the same primer (static) section, then we only need to do one search, not 50. In
practice, this means we can group DELs based on having the same static section, and then only do a search per group
rather than per library.

Of course this causes a new issue: how to tell which library the read is from when the generated match could be for any
of the libraries.
All that is needed to fix this is to use the alignment to locate the library section (just like we would
for any other section). Then we can use the :ref:`barcode calling module <barcode-calling>` to call the library section
and determine which library this read is from. This requires a small tweak to the static alignment method, which instead
will only do a partial alignment for just the library section, rather than the full barcode design. After we have the
library a second and full alignment can be generated between the read and the full library barcode.

Yet this tweak also causes a new issue: what if the libraries share a static section, but the distance between it and the library
section is not the same. Then the partial library section only alignment will not work for all the DELs in this group.
Instead, we need to further sub group the
DELs. These sub groups are determined by the distance between the static section and the library tag. If the distance
is the same, between DELs, they can be in the same subgroup. Otherwise they are in separate subgroups. Just like how
one search must be made per group, one partial library alignment must be made per subgroup in a group.

This might seem a bit too much to handle, having to group DELs based on various barcode schema compatibilities.
Lucky for us DELi can do this on its own. Given a set of libraries it can find the set of static sections to use that
will minimize the number of groups and subgroups (thus minimizing the number of searches needed). In practice, you
can often speed up a decoding run about 20-40x by using ``single`` mode over ``library``. However, in cases where
none of the libraries share any static sections, most of the benefits of this method are lost, and ``library`` mode
might be a better choice, especially if the size of the static sections are small and library barcodes are larger.

.. _flanking-section-demulti:

Flanking primer demultiplexing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The last method extends the single section demultiplexing method by using two static sections instead of one.
The benifits of this are discussed in the :ref:`static section alignment <static-section-alignment>` section, but the
main takeaway is this further improves the your ability to over come INDELs, at the cost of a minor increase risk in
:ref:`miscalling of barcodes <barcode-calling>`. DELi can also automatically determine which static sections to use as
the flanking sections, so you don't have to worry about that either.

Core demultiplexing algorithms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To actually carry out the above demultiplexing approaches DELi uses existing tools to perform the alignment and searches
As mentioned already, Biopython's `PairwiseAligner <https://biopython.org/docs/dev/Tutorial/chapter_pairwise.html>`_
is used for :ref:`full barcode demultiplexing <full-barcode-demulti>` and can be used for the other approaches
as well if you set `realign = True` (which will reject the search based alignment and replace it with a
semi-global alignment).

For the search based approaches, DELi can either use Cutadapt's
`Adapter search <https://cutadapt.readthedocs.io/en/stable/>`_ or the Regex packages
`fuzzy search <https://github.com/mrabarnett/mrab-regex>`_

.. note::
    ``regex`` is not the same as ``re`` from the python standard library. It is a third party library that
    supports more complex regular expression syntax (but is still compatible with python ``re``).

Cutadapt
~~~~~~~~
This will use Cutadapt's alignment algorithm to find the static section(s) in the read.
Cutadapt is a fantastically optimized tool for locating adapters in sequencing reads. However,
adapters are meant to be at the ends of reads, while static sections could be anywhere in the read.
Usually, primers are at the ends of reads so things work out okay, but if you have static sections
in the middle of the read, Cutadapt might not be able to find them as well as other methods.
This is because the alignment algorithm prioritizes finding matches at the ends of reads, and will
pick them over matches in the middle of reads *even if* there match in the middel is better (less errors).
This is not a bug but a feature of cutadapt, as it is designed for adapter trimming, not general sequence
searching. However, it can cause strange behavior in DELi. Only use this approach if you are sure your static
sections are at the ends of your reads (like primers). DELi will warn you if it detects static sections
in the middle of reads when using Cutadapt.

.. warning::
    If you are considering using cutadapt over regex for demultiplexing, you should read the above
    paragraph, there are important caveats to be aware of.

Regex
~~~~~
This will conduct a fuzzy regular expression search for the static section(s) in the read. Unlike a normal
regular expression which cannot handle insertion, deletion or substitution errors appearing arbitrarily in a
specific section of the query, the `regex` library supports fuzzy searching. DELi will build out regex expressions
that look like ``(?:(?P<primer1>AGCTAGCT)){e<=1}.{{10,20}}(?:(?P<primer2>CGTACGTA)){e<=1}``. While confusing to read,
these queries are quick powerful and efficient. In nearly all cases, regex is faster at finding static sections in a
given read than any other method.

Regex is the default and recommend method to use.

.. _barcode-calling:

Barcode Calling
---------------
After the alignment is complete and we know which spans of the read contain which sections of the barcode, we can
call the decodable sections to determine the final ID associated with the read. On the surface, this is pretty straight
forward, as we just need to loop through each decodable section spans and look up its match from the available options
provided. Yet is isn't so simple.

First, looking for perfect matches leaves alot on the table. Depending on sequencing quality and
number of decodable section as many as 30-40% of reads could have an error that would cause them to be uncallable if
a perfect match was needed. For DELs this is especially wasteful, as the set of valid barcodes is often far smaller than
the set of possible barcodes. This means even with a single error in the barcode we can be pretty confident about
what the original barcode was and still call it correctly.

To address this DELi uses a hashmap based approach, where prior to decoding a list of each barcode and all its possible
1-Levenshtein distance (1 INDEL or SNP) neighbors is generated and stored. This allows for O(1) lookup of the correct
barcode even if there is an error in the read, as long as it is only 1 error away from a valid barcode. If there is
a case where two of the valid barcodes are similar enough that they share a 1-Levenshtein distance neighbor
(making it too ambiguous to know the correct barcode), DELi can detect this and correct all the read as non-decodable.
This is the ``ambiguous`` mode in the error correcting settings.

There is one more "problem" to address. If you read the above sections about
:ref:`library demultiplexing <library-demultiplexing>`, you will notice that the alignments generated are not
guaranteed to be the correct size, or even just 1 bp larger. This would automatically cause the barcode calling to fail,
as the hashmap does not cover that many Levenshtein edits. Worry not, DELi is designed to handle this, and in fact,
it is a core feature in making barcode calling more robust.

When DELi get asked to call a section. It is given the aligned span from the read. It will then sweep through all possible
windows of the expected size, plus one smaller and one larger (to account for possible INDELs) for the section to try and
find a match in the hashmap.
There are several modes for this, but the default is a "greedy first perfect" which will pick the first window with a perfect
(zero errors) match or pick the window that results in a match with the least number of errors (if no perfect match is found).
This is still really fast, as
the matching is an O(1) operation thanks to the hashmap callers. Using this DELi can recover reads with as many as 2-3
INDELs with high accuracy.

This also address another problem we have yet to address, alignment spans that overlap. If decodable barcode sections
are close, and there are INDELs, it is possible for the alignments to overlap. By sweeping the windows, DELi can find
an alignment that either avoids overlaps between matching windows (which would result in an impossible decode output
since one base pair can't be claimed by two different sections) or allows DELi to fail when an non-resolvable
overlap is detected. This is generally quite rare, but DELi is designed to handle it when it does happen.

Barcode Decoding
----------------
Once the read has been assigned to a specific library and an alignment(s) between the read and the library
are created, then the sections can be decoded. This is pretty straight forward. DELi will loop through each alignment
and try to decode the building block barcodes based on the mapped nucleotides for each variable section.
If it fails (because that barcode doesn't exist) it will try the next alignment and continue until it finds a valid
decoding or runs out of alignments. Hashmap based :ref:`barcode calling <>` is used to allow for possible errors in the
barcode sections when decoding while keep runtime as low as possible (O(1) in this case).

Once all the variable sections are called, DELi will then locate the UMI by both using the alignment to find
where the UMI is located in the read *and* looking taking the closest aligned variable section to the UMI and
searching *n* bp away to find the UMI sequence. Both these (if different) are save as possible UMIs for the read,
with the *n* bp away being the primary UMI call. Future updates will build upon this to do more advanced UMI
correction.

This then produced a full decoded read, mapping it to a unique compound ID and associating it with a UMI.


UMI extraction
--------------
After the compound ID is decoded, the alignment is also used to extract the UMI sequence from the read.
DELi *does not* use the alignment for this. Instead it uses the known distance between the UMI and the
closest decoded section to find the UMI sequence. This is to adjust for INDELs that may have occurred in the read,
since window scanning does not work on the UMI (there is nothing to call, the region is meant to be random).


Decode Parallelization and Collection
-------------------------------------
The process of decoding a single read into a compound ID and UMI is embarrassingly parallel,
as each read can be decoded independently of the others. This is why DELi's decoding algorithm does not
implement any native multithreading or multiprocessing; it is easy (and way faster) to just run multiple processes of
DELi on different chunks of the input FASTQ files.

However, before count and UMI degeneration can occur, the decoded compounds need to be aggregated; we need to know
how many time each compound showed up and with what UMIs. This process is not easily parallelizable, since all the
data must be looked at in a single process (at some point). This is also the reason this process is sepatated from
decoding. DELi can do this collection for you; it is just a
massive groupby and count operation. Polars' streaming engine is used for this to keep the memory usage low.

UMI Degeneration
----------------
Once the compounds are collected and counted, DELi can do UMI degeneration to correct for PCR amplification bias.
Technically this is optional, but non-UMI corrected counts are not very useful due to all the noise introduced by
PCR amplification bias. DELi implements two methods for UMI degeneration, :ref:`Exact <exact-umi-degeneration>` and
:ref:`Clustered <clustered-umi-degeneration>`.

.. _exact-umi-degeneration:

Exact UMI Degeneration
^^^^^^^^^^^^^^^^^^^^^^
This is the simplest method, where all reads with the same UMI are grouped together and their counts are summed.
There is no room for error in the UMI. Given that even high end sequencers have error rates of about 1/100 base pairs,
This is going to introduce some noise, as AACCTTGG and AACCTTG*A* will be counted as the **different** UMIs, even though,
statistically speaking, the odds are far more likely one is a misread of the other.

To demonstrate this, say we have 100bp DNA barcode with a 10 bp UMI. on average 1 out of every 10 reads (10%) will have a single error
in the UMI, giving it a hamming distance on 1. Now, with 10bp UMIs, there are 1,048,576 possible UMIs, but for any given UMI, there
are only 30 other UMIs with a hamming distance of 1, or a 0.0029% chance of a read being misread as one of those. If we assume that the
chance there are far more possible UMIs that those actually in the DEL (Since DELs are synthesized at such low concentrations) it is 4 fold
more likely that two UMIs with hamming distance of 1 are just one single UMI (one has a misread) than it is to be two unique UMIs.
In practice these odds are even lower, as UMIs are often longer than 10bp.

.. _clustered-umi-degeneration:

Clustered UMI Degeneration
^^^^^^^^^^^^^^^^^^^^^^^^^^
To compensate for misreads and PCR errors during sequencing that end up in the UMI section, DELi implements a
"clustered" UMI degeneration method. It is an implementation of the "directional" method outlined in Smith et al [1]_.

References
----------
.. [1] Smith, T., Heger, A., & Sudbery, I. (2017). UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Research, 27(3), 491-499. doi:10.1101/gr.209601.116
