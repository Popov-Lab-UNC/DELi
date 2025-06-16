.. _deli-config-docs:

DELi Configuration
==================

DELi is built with Poetry, and subscribes to the notion of 'no arbitrary code during install'.
By default no DELi figuration directory or file will be created.
If you want to control these defaults, you can use ``deli init_config <DESIRED_LOCATION>``
to generate a default configuration file.

By default, DELi looks for a ``.deli`` config folder in your ``$HOME`` directory.
You can override where to look for this file by setting the ``DELI_CONFIG_DIR`` environment
variable.

This folder contains a directory called ``deli_data``, which is the default DELi Data Directory
(See :ref:`deli-data-dir-ref`) and default DELi configuration file named ``.deli``

.. _deli-config-file-ref:

The ``.deli`` config file
-------------------------
This file follow the config (.ini) format and controls some of the default setting DELi uses
when decoding. This controls most of the more advanced settings that you might need for
DELi applications (like ``deli decode``). Most users will not have to alter these, outside
of setting the ``DELI_DATA_DIR``.

For most command line tools you can pass a ``--deli-config`` argument to specify
a config file. If no file is specified DELi will look to load one from ``$HOME\.deli\.deli``.
If there is no such file, DELi will use the default settings

DELI_DATA_DIR
^^^^^^^^^^^^^
This should be the path to the DELI_DATA_DIR, which defaults to ``$USER\.deli\deli_data``.
You can override this by setting the environment variable ``DELI_DATA_DIR``

BB_MASK
^^^^^^^
This is the token used by DELi when masking building blocks.
Masking building blocks is done when creating IDs for non-full DEL compounds,
like a disynthon for a trisynthon. The default is "###"

MAX_INDEX_RISK_DIST_THRESHOLD
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is the max Levenshtein distance a called index can be from the observed sequence.
During calling, DELi will pick the index with the shortest Levenshtein distance to the observed
sequence as the correct index (unless there is a tie, in which case it fails).
This distance might be very big. This parameter puts a cap to that distance.
DELi will ignore this is there is only 1 index to call.
The default value is 3

MAX_LIBRARY_RISK_DIST_THRESHOLD
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is the max Levenshtein distance a called library can be from the observed sequence.
During calling, DELi will pick the library with the shortest Levenshtein distance to the
observed sequence as the correct library (unless there is a tie, in which case it fails).
This distance might be very big. This parameter puts a cap to that distance.
DELi will ignore this is there is only 1 library to call.
The default value is 4

NUC_2_INT
^^^^^^^^^
This is a mapping of nucleotide to integer for all 4 bases.
This solely used for hamming decoding (and encoding) as the bases are converted to numbers
to calculate the parity. You should make sure yours matches DELi's.

The default is ``A:0,T:1,C:2,G:3``

The format is simply a comma seperated list of the 4 bases: A,T,G,C seperated by a colon to
their numeric conversion: 0,1,2,3. They can appear in any order, but all 4 must be present.
If they are duplicated an exception will be raised. If you use numbers outside of 0,1,2,3 an
exception will be raised.
