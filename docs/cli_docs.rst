.. _deli-cli-docs:

===========================
DELi Command Line Interface
===========================

While you can use DELi as a normal package in python, DELi provides a
handful of command line tools to make common tasks more straight forward to run.

Config commands
---------------
``deli config`` provides a set of commands for managing the DELi config file

``deli config init``:
^^^^^^^^^^^^^^^^^^^^^
Create a default DELi configuration file

PATH is the path to the deli config directory to initialize. If not
provided, defaults to ~/.deli
.. code-block:: text
    Usage: deli config init [OPTIONS] [PATH]

    Options:
    -o, --overwrite  Overwrite any existing config directory
    --help           Show this message and exit.

DELi Data Directory Commands
----------------------------
``deli data`` provides a set of commands for managing the DELi Data Directory

``deli data init``
^^^^^^^^^^^^^^^^^^
Initialize the configuration directory

PATH is the path to the deli data directory to initialize.

NOTE: fix-missing will not overwrite existing sub-directories, while
overwrite will. 'overwrite' will also add any missing sub-directories, but
also replace existing ones.
.. code-block:: text
    Usage: deli data init [OPTIONS] [PATH]

    Options:
    -f, --fix-missing  Fix a deli data directory missing sub-directories
    -o, --overwrite    Overwrite any existing data directories
    --help             Show this message and exit.

``deli data set``
^^^^^^^^^^^^^^^^^
Set the DELi data directory to use for decoding

PATH is the path to the deli data directory to set.

NOTE: if not using --update-config, you will need to set the DELI_DATA_DIR
environment variable manually; the command required will be printed after
running.
.. code-block:: text
    Usage: deli data set [OPTIONS] PATH

    Options:
    -u, --update-config  Update the DELi config to use this data directory as default
    --help               Show this message and exit.

``deli data which``
^^^^^^^^^^^^^^^^^^^
Print the current DELi data directory.
If the DELi data directory is not set, will print a message and exit with a
non-zero status code.
.. code-block:: text
    Usage: deli data which [OPTIONS]

    Options:
      --help  Show this message and exit.


Decoding
--------
``deli decode`` is used to take a decoding experiment and fastq file from a DEL
selection and run decoding.
.. code-block:: text
    Usage: deli decode [OPTIONS] DECODE_FILE [FASTQ_FILES]...

      Run decoding on a given fastq file of DEL sequences

      DECODE is the path to a YAML file describing the decoding run settings.
      FASTQ_FILES is a path, list of paths or glob of FASTQ files to decode.

      NOTE: if the DECODE file contains a `selection` field, it will be used to
      select the

    Options:
      -o, --out-dir PATH        Output directory, defaults to CWD
      -i, --ignore-decode-seqs  Ignore the fastq sequence files in the decode file
      -p, --prefix TEXT         Prefix for output files
      -t, --tqdm                Show tqdm progress
      --save-failed             Save failed decoding results to a separate file
      --save-counter            Save decoding counters as a JSON file
      --debug                   Enable debug mode
      --disable-logging         Turn off DELi logging
      --skip-report             Skip generating the decoding report at the end
      --deli-data-dir PATH      Path to DELi data directory to read libraries from
      --help                    Show this message and exit.

.. _deli-enumeration-cli-docs:

Enumeration
-----------
``deli enumerate`` is used to enumerate a DEL library to generate all possible compounds
and write them to a file.
.. code-block:: text
    Usage: deli enumerate [OPTIONS] LIBRARY_FILE

      Enumerates compounds from a given library

      If out_path is not provided, will save to the current working directory as a
      CSV file named <library_id>_enumerated.csv

      LIBRARY_FILE is the path to a DELi library file to enumerate.

    Options:
      -o, --out_path PATH  Output CSV file path
      -t, --tqdm           Enable TQDM progress bar
      --help               Show this message and exit.
