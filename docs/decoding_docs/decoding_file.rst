.. _decoding-run-file-docs:

=================
Decoding Run File
=================

Decoding run files are a special form of a :ref:`selection file <selection-file-docs>` that DELi uses to configure a
decoding run. They are based on YAML and look like this:

.. code-block:: yaml

    selection_id: "TEST_SELECTION"
    target_id: "Target1"
    selection_condition: "200uL Library with 50pM Protein"
    date_ran: "2025-03-30"
    additional_info: "This is a test selection for the DELi decode."

    sequence_files:
      - "./test_data/example.fastq"
      - "./test_data/example.fastq"

    libraries:
        - "DEL004"
        - "DEL005"

    decode_settings:
        demultiplexer_mode: "single
        demultiplexer_algorithm: "regex
        realign: NO
        revcomp: YES


Basically this is a selection file with a few extra keys that DELi uses to configure the decoding run settings.
The most important of these is the ``sequence_files`` key, which is a list of paths to the raw sequence files that DELi
needs to decode. The paths to the sequence files can be relative or absolute (absolute is recommend).
If they are relative, they are relative to the current working directory.

.. note::
    DELi does not handle demultiplexing of sequence files to separate selection conditions.
    If your sequencer does **not** demultiplex your selection after sequencing (most do),
    you will need to do that yourself before decoding.
    All listed sequence files should be for the same selection, and DELi will treat them that way.

The ``libraries`` key is the same as in a selection file, and it lists the libraries that were used in the selection.
DELi will also look for these files in the :ref:`DELi Data Directory <deli-data-dir-ref>` if they are not found in the
current directory.

Decoding settings
-----------------
The ``decode_settings`` key is a dictionary that contains the settings for the decoding run.
Every decode setting is optional as DELi has default values for all of them
The following settings are available:

- ``demultiplexer_mode``: Literal["cutadapt", "regex", "full"], default = "regex"
    The demultiplexing algorithm to use

    - "cutadapt": use a cutadapt to locate sections
    - "regex": use a regular expression based demultiplexer
    - "full": use a full alignment based demultiplexer

- ``demultiplexer_sections``: Literal["library", "single", "flanking"], default = "flanking"
    The demultiplexing section strategy to use

    - "library": demultiplex by matching just the library tag
    - "single": demultiplex by matching a single static barcode section
    - "flanking": demultiplex by matching barcode sections that flank the library tag

    (flanking means one before and one after the tag)
- ``realign``: bool, default = False
    if true, will perform a local realignment of the read to the
    libraries barcode schema *after* demultiplexing determine the library.
    This could help recover reads that have complex alignments due multiple indels
- ``library_error_tolerance``: int, default = 1
    The number of errors you are willing to tolerate in any given barcode
    section during library demultiplexing. Will apply to each section
    independently. For example, a flanking demultiplexer will allow for
    1 error in *each* of the flanking sections.
- ``library_error_correction_mode_str``: str, default = "levenshtein_dist:2,asymmetrical"
    The error correction mode string to use for library barcode
    calling during demultiplexing.
- ``min_library_overlap``: int , default = 8
    if using a cutadapt style demultiplexer, this is the minimum number of bases
    that must align to the expected barcode section for a match to be called.
    See the cutadapt documentation for more details on this parameter.
- ``wiggle``: bool, default = False
    if true, will extend aligned sections by 1 bp on each side
    and then and scan all possible chunks of the expected barcode length,
    1 smaller and 1 larger than expected length (in that order). If not
    using a local realignment post demultiplexing this can help recover
    reads lost to indels in the barcode region.
- ``revcomp``: bool, default = False
    If true, search the reverse compliment as well.
    In most cases it is faster to use an external tools
    to align and reverse compliment reads before decoding
- ``max_read_length``: int or None, default = None
    maximum length of a read to be considered for decoding
    if above the max, decoding will fail
    if `None` will default to 5x the min_read_length
- ``min_read_length``: int or None, default = None
    minimum length of a read to be considered for decoding
    if below the min, decoding will fail
    if `None` will default to the smallest min match length of
    any library in the collection considered for decoding
    with 10bp of buffer
- ``default_error_correction_mode_str``: str, default = "levenshtein_dist:1,asymmetrical"
    The default error correction mode string to use for decoding.
    If a barcode section lacks a specified error correction mode,
    this mode will be used.
    See the documentation for more details on the format of this string.
