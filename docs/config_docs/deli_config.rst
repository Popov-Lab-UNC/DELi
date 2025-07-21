.. _deli-config-docs:

DELi Configuration
==================

DELi uses a ``.deli`` config file to control some of its more advanced settings.
To avoid arbitrary code execution, DELi does not initialize the config file when
installed via ``pip``. Instead it will create one the first time you every use DELi
(and it will raise a warning to let you know it did so).
You can also create one yourself using ``deli config init``.
A valid DELi config file is required for DELi to run *or* be imported.
DELi will raise a ``DELiConfigError`` if it detects something wrong with the config file.

.. note::
    Any settings that might impact the outcome of a DELi command (like ``deli decode``)
    are not part of the config file. The config file if for settings related to how DELi
    outputs data, where it looks for files, and defining mappings.

File Format
-----------
The config file follows the config (.ini) format, with sections and key-value pairs.
DELi looks for a handful of sections:

- ``deli.data``: holds settings related to the DELi Data Directory
- ``deli.buildingblocks``: holds settings related to how to handle building blocks
- ``deli.hamming``: holds settings related to hamming encoding/decoding

.. warning::
    The config file is case sensitive.

``deli.data``
^^^^^^^^^^^^^
This section has only one key: ``deli_data_dir``. This is the path to the
:ref:`DELi Data Directory <deli-data-dir-ref>`. This is the default location DELi will look
for your DELi data files. All other methods of setting the DELi Data Directory (like with
the ENV variable or the ``set_deli_data_dir`` function) will override this setting.

``deli.buildingblocks``
^^^^^^^^^^^^^^^^^^^^^^^
This section has only one key: ``bb_mask``. This is the token DELi will use when masking
building blocks. This happens during things like Di/Monosynthon ID creation, where specific
cycles are aggregated (and then masked). The default is "###".

. _deli-config-hamming-section

.. _deli-config-hamming-section:

``deli.hamming``
^^^^^^^^^^^^^^^^
This section has one key: ``nuc_2_int``. This is the mapping of nucleotide to integer values
Hamming encoding/decoding happens using quaternary number space, thus must first covert the DNA
into the numbers 0, 1, 2, and 3. How this is done must be the same for encoding and decoding.
If you used a specific mapping when encoding, you must use the same mapping when decoding and
this is how you tell DELi to change that. The default is ``A:0,T:1,C:2,G:3``.

.. note::
    The ``nuc_2_int`` mapping is only used for Hamming encoding/decoding using a hamming matrix.
    If you are using the "random" decoding method, this setting is not used.

The hamming section also has an infinite number of sub sections, with each being the name of a
specific hamming matrix. These follow the naming convention ``deli.hamming.<NAME>``. They
each have to keys to fill out: ``hamming_order`` and ``custom_order``. You can read more about
how this works in the :ref:`custom hamming docs <deli-custom-hamming-docs>`.

.. warning::
    The ":" character is reserved for the :ref:`error correction parser <error-correction-docs>`,
    thus DELi will not allow the use of it in the names of hamming matrices. If you try to use it,
    DELi will raise a ``DELiConfigError``.
