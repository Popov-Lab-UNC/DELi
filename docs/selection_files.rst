.. _selection-file-docs:

===============
Selection Files
===============
DEL selections are experiments, and like all experiments, they have all kinds of metadata and conditions that are
associated with them. Some of this information is useful for DELi to know when processing the selection data.
To help streamline this DELi has a "selection file format" that allows you to easily enter information about your
DEL selection.

.. _selection-file-format:

Selection File Format
---------------------
The selection file format is a simple YAML file that contains at a bare minimum the following information:

- **selection_id**: An identifier for the selection. This should be a short, descriptive name that identifies the selection
  (preferably the same as any other experiment tracking you use).

- **libraries**: A list of the libraries that were used in the selection.
  They can be the name of libraries in the :ref:`DELi data directory <deli-data-dir-ref>` or the (full)
  path to the library definition file for that library

Any other sections are not (currently) used by DELi. But you can add whatever info you want. It can be a good way to
track your experiments and keep all the relevant information in one place. For example, you might want to include
information about the target protein, the selection conditions, or any other relevant details.

.. code-block:: yaml

    selection_id: "EXAMPLE_SELECTION"

    libraries:
        - "DEL004"
        - "DEL005"


.. _sequenced-selection-file-format:

Sequenced Selection File Format
---------------------
If your selection has been sequenced, you can store the location of the sequence data in the selection file as well.
DELi calls this a "SequencedSelection" and when it is load will automatically set up readers to parse the FASTQ files.

You just need to add a "sequence_files" section to the selection file with the paths to the FASTQ files. For example:

.. code-block:: yaml

    selection_id: "EXAMPLE_SELECTION"

    libraries:
        - "DEL004"
        - "DEL005"

    sequence_files:
        - "path/to/selection1.fastq"
        - "path/to/selection2.fastq"

.. note::
    If you are using the sequence file to run a decoding job, it **must** be a sequenced selection file.
    DELi needs to know where the sequence files are to run the decoding job.

.. _decode-settings-in-selection-file:

Decode Settings in the Selection File
---------------------
In the case that you are using your selection file to run a :ref:`decoding job <running-decode-docs>`, you can also
include the decoding settings in the selection file itself. You can add them under the "decode_settings" key. For example:

.. code-block:: yaml

    selection_id: "EXAMPLE_SELECTION"

    libraries:
        - "DEL004"
        - "DEL005"

    sequence_files:
        - "path/to/selection1.fastq"
        - "path/to/selection2.fastq"

    decode_settings:
        ignore_tool_compounds: true
        demultiplexer_algorithm: "cutadapt"
        demultiplexer_mode: "library"
        realign: true
        library_error_tolerance: 2
        library_error_correction_mode_str: "levenshtein_dist:2,asymmetrical"
        min_library_overlap: 8
        revcomp: true
        library_wiggle: true
        wiggle: true
        decode_matching_approach: "first_perfect"
