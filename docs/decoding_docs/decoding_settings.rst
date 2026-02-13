.. _decoding-settings-docs:

=================
Decoding Settings
=================

There are several settings that can be configured for a decoding run in DELi.
They are implemented in the ``DecodeSettings`` class in the ``deli.decode`` module.
These settings can be customized by adding a ``decode-settings`` key to you :ref:`selection file <selection-file-docs>` you submit for decoding,
or by passing a yaml file with the settings to the :ref:`decode CLI <decode-cli-docs>`

Every decode setting is optional as DELi has default values for all of them
The following settings are available:

- ``ignore_tool_compounds``: bool, default = ``False``
    if true, will ignore any tool compounds during decoding
- ``demultiplexer_algorithm``: Literal["cutadapt", "regex", "full"], default = "regex"
    The demultiplexing algorithm to use.
    - "cutadapt": use a cutadapt to locate sections
    - "regex": use a regular expression based demultiplexer
    - "full": use a full alignment based demultiplexer
- ``demultiplexer_mode``: Literal["library", "single", "flanking"], default = "flanking"
    The demultiplexing section strategy to use.
    - "library": demultiplex by matching just the library tag
    - "single": demultiplex by matching a single static barcode section
    - "flanking": demultiplex by matching barcode sections that flank the library tag
    (flanking means one before and one after the tag)
- ``realign``: bool, default = ``False``
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
- ``revcomp``: bool, default = True
    If true, search the reverse compliment as well.
    In most cases it is faster to use an external tools
    to align and reverse compliment reads before decoding
- ``library_wiggle``: bool, default = ``False``
    if true, will extend library aligned sections by 1 bp on each sidd.
    Similar to wiggle, but used during library demultiplexing.
- ``wiggle``: bool, default = ``False``
    if true, will extend aligned sections by 1 bp on each side during decoding.
    Can help recover reads with INDELs if not using realignment after demultiplexing
    (and is much faster).
- ``decode_matching_approach``: Literal["greedy", "first_best", "first_perfect", "search_all"] = "first_perfect",
    the approach to use when matching decoded sections to the library during decoding.
    - "greedy": find the best scoring match across all possible windows and return it, even if imperfect
    - "first_best": search all possible windows and return the first occurring best match, even if more
    than one window has the same best score
    - "first_perfect": search all possible windows and return the first perfect match found, and otherwise
    fail if more than one window (with non-perfect score) share the same best score
    - "search_all": search all possible windows and return the best scoring one. Unlike 'first-perfect', if
    multiple windows have a perfect score the read will be rejected as ambiguous rather than just taking
    the first perfect match found.
- ``max_read_length``: int or None, default = ``None``
    maximum length of a read to be considered for decoding
    if above the max, decoding will fail
    if `None` will default to 5x the ``min_read_length``
- ``min_read_length``: int or None, default = ``None``
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
