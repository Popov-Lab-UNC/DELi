DELi Configuration
==================

Upon installing DELi with ``pip``, a config folder named ".deli" should appear in your
USER directory. This folder contains the default ``DELI_DATA_DIR`` and the ``.deli`` config file

``.deli`` config file
-----------------
This file follow the config format and controls some of the default setting DELi uses
when decoding.

- DELI_DATA_DIR
This should be the path to the DELI_DATA_DIR, which defaults to ``$USER\.deli\deli_data``.
This is the last place DELi looks for the DELI_DATA_DIR, the evironment variable $DELI_DATA_DIR (if set and existing) will override this location

- BB_MASK
This is the token used by DELi when masking buiding blocks.
Masking building blocks is done when creating IDs for non-full DEL compounds,
like a disynthon for a trisynthon. The default is "###"

- MAX_INDEX_RISK_DIST_THRESHOLD
This is the max Levenshtein distance a called index can be from the observed sequence.
During calling, DELi will pick the index with the shortest Levenshtein distance to the observed sequence as the correct index (unless there is a tie, in which case it fails).
This distance might be very big. This parameter puts a cap to that distance.
DELi will ignore this is there is only 1 index to call.
The default value is 3

- MAX_LIBRARY_RISK_DIST_THRESHOLD
This is the max Levenshtein distance a called library can be from the observed sequence.
During calling, DELi will pick the library with the shortest Levenshtein distance to the observed sequence as the correct library (unless there is a tie, in which case it fails).
This distance might be very big. This parameter puts a cap to that distance.
DELi will ignore this is there is only 1 library to call.
The default value is 4

- NUC_2_INT
This is a mapping of nucleotide to integer for all 4 bases.
This solely used for hamming decoding (and encoding) as the bases are converted to numbers
to calculate the parity. You should make sure yours matches DELi's.

The default is ``A:0,T:1,C:2,G:3``

The format is simply a comma seperated list of the 4 bases: A,T,G,C seperated by a colon to their numeric conversion: 0,1,2,3. They can appear in any order, but all 4 must be present.
If they are duplicated an exception will be raised. If you use numbers outside of 0,1,2,3 an
exception will be raised.
