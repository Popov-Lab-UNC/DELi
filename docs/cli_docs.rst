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
    Usage: deli config init [OPTIONS] [PATH]
    Options:
    -o, --overwrite  Overwrite any existing config directory
    --help           Show this message and exit.

    -o, --overwrite  Overwrite any existing config directory
    --help           Show this message and exit.
``deli data`` provides a set of commands for managing the DELi Data Directory

``deli data init``
^^^^^^^^^^^^^^^^^^
``deli data`` provides a set of commands for managing the DELi Data Directory

PATH is the path to the deli data directory to initialize.

Initialize the configuration directory
overwrite will. 'overwrite' will also add any missing sub-directories, but
PATH is the path to the deli data directory to initialize.
.. code-block:: text
NOTE: fix-missing will not overwrite existing sub-directories, while
overwrite will. 'overwrite' will also add any missing sub-directories, but
also replace existing ones.
    -o, --overwrite    Overwrite any existing data directories

``deli data set``
^^^^^^^^^^^^^^^^^
Set the DELi data directory to use for decoding
    -f, --fix-missing  Fix a deli data directory missing sub-directories
    -o, --overwrite    Overwrite any existing data directories
    --help             Show this message and exit.
NOTE: if not using --update-config, you will need to set the DELI_DATA_DIR
environment variable manually; the command required will be printed after
running.
.. code-block:: text
Set the DELi data directory to use for decoding

PATH is the path to the deli data directory to set.
    -u, --update-config  Update the DELi config to use this data directory as default
NOTE: if not using --update-config, you will need to set the DELI_DATA_DIR

``deli data which``
Print the current DELi data directory.
non-zero status code.
.. code-block:: text
    Usage: deli data which [OPTIONS]

    -u, --update-config  Update the DELi config to use this data directory as default
    --help               Show this message and exit.


Decoding
--------
``deli decode`` is used to take a decoding experiment and fastq file from a DEL
selection and run decoding.
.. code-block:: text


      DECODE is the path to a :ref:`YAML file describing the decoding run settings <decoding-run-file-docs>`.
      FASTQ_FILES is a path, list of paths or glob of FASTQ files to decode.

      NOTE: if the DECODE file contains a `selection` field, it will be used to
      select the

    Options:
      -o, --out-dir PATH        Output directory, defaults to CWD
      -i, --ignore-decode-seqs  Ignore the fastq sequence files in the decode file
      -p, --prefix TEXT         Prefix for output files
      --save-failed             Save failed decoding results to a separate file
      --save-counter            Save decoding counters as a JSON file
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
