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
        library_error_tolerance: 0.1
        min_library_overlap: 8
        revcomp: YES
        disable_error_correction: NO


Basically this is a selection file with a few extra keys that DELi uses to configure the decoding run settings.
The most important of these is the ``sequence_files`` key, which is a list of paths to the raw sequence files that DELi
needs to decode. The paths to the sequence files can be relative or absolute.
If they are relative, they are relative to the current working directory.

.. note::
    DELi does not handle demultiplexing of sequence files to separate selections. If your sequencer does not delmultiplex the files
    for you, you will need to do that yourself before decoding. All listed sequence files should be for the same selection.

The ``libraries`` key is the same as in a selection file, and it lists the libraries that were used in the selection.
    DELi will also look for these files in the :ref:`DELi Data Directory <deli-data-dir-ref>` if they are not found in the current directory.

Decoding settings
-----------------
The ``decode_settings`` key is a dictionary that contains the settings for the decoding run.
Every decode setting is optional as DELi has default values for all of them
The following settings are available:
- ``library_error_tolerance``: default = 0.1; the maximum allowed error rate when calling a library. If between (0,1)
  DELi considers it a percentage, otherwise it considers it a as the absolute number of errors allowed when making a call.
- ``min_library_overlap``: default = 10; the amount of the overlap between the library tag and the read tag that is required to call a library.
  Should never be smaller than 8 unless you really know what you are doing.
- ``revcomp``: bool = False; search the reverse complement of each read for library tags. Only useful if your reads don't have a known direction.
  Illumina sequencing reads are always in the same direction, so this is *not* needed for those.
- ``max_read_length``: default = 5 times the longest barcode among the library collection; this is the longest read length DELi will attempt to
  decode. If longer, will fail to decode.
- ``min_read_length``: default = shortest barcode among the library collection; this is the shortest read length DELi will attempt to
  decode. If shorter, will fail to decode.
- ``disable_error_correction``: default = False; if True, DELi will not apply any error correction methods when decoding the reads
  (even if they are specified in the library file).
