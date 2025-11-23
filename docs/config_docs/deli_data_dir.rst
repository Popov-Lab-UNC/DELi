.. _deli-data-dir-ref:

DELi Data Directory
===================
Large DEL operations might have hundreds of libraries and even more building block sets.
Each of these will need a configuration file to define them for DELi.
To help organize all these files, DELi has something called the "DELi Data Directory".
This directory is the root source for all DELi files, and when configured, DELi will
look for files here by default. This also allows you to reference libraries and building block
sets by name instead of full path.

This feature is not required to use DELi, but it is highly recommended to help with managing files
related to your DELs. If you have 100 DELs that could be 300 build block sets. That is 400 files for needed
to configure DELi (and frankly to just manage your DELs in general). The DELi Data Directory help make this a
bit more organized and avoid issues with file paths and references.

.. note::
    On installation, DELi will not create a default DELi Data Directory anywhere in your filesystem.
    If you want to use the DELi Data Directory, you will need to create one yourself.

Creating a DELi Data Directory
------------------------------
A valid DELi Data Directory must have the following sub folders:

* ``libraries``
* ``building_blocks``
* ``reactions``
* ``tool_compounds``

Additional files and folders can be added, but these DELi will ignore these.
You can create an empty DELi Data Directory by just creating a directory with the
two named subfolders *or* by running the command
``deli data init PATH/TO/DATA/DIR``.

As their names suggest, the two sub folders should hold the definitions files
the libraries (in ``libraries``) and building block sets (in ``building_blocks``).
The both are treated as globs, meaning any file in any subdirectory that is a valid definitions
file can be found and used by DELi. For example, the following file tree would work
just fine with:

.. code-block:: text

    my-deli-data-dir/
    ├── SOME_INFO.md
    ├── other-stuff/
    ├── libraries/
    │   ├── old-libraries/
    │   │   ├── lib1.json
    │   │   └── lib2.json
    │   ├── new-libraries/
    │   │   ├── lib3.json
    │   │   └── lib4.json
    │   └── lib5.json
    ├── buildingblocks/
    │   ├── junk/
    │   │   ├── junk/
    │   │   │   ├── junk/
    │   │   │   │   ├── junk/
    │   │   │   │   │   └── my-bb-set1.csv
    │   │   │   │   └── my-bb-set2.csv
    │   │   │   └── my-bb-set3.csv
    │   │   └── my-bb-set4.csv
    │   └── my-bb-set5.csv
    └── more-stuff/

Fixing a DELi Data Directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
DELi will check that any path you provide as a DELi Data Directory is valid.
If it detects its not, a ``DELiDataDirError`` will be raised. Generally, this
will tell you what is wrong with the directory and suggest how to fix it.
You can also use the command ``deli data init --fix-missing PATH/TO/BROKEN`` to
try and fix a broken DELi Data Directory. See the command :ref:`deli-cli-docs` for more info.

How to use the Data Dir
-----------------------
The purpose of the Data Directory is to simplify referencing other DEL components when
defining other components. For example, libraries require references to Building Block Sets.
Instead of using the full paths you can just specify the name of the Building Block set.
This helps to both centralize file locations and make the definition files easier to write.
For example, instead of referencing the BB sets in the library json file as
``\home\my_user\del\data\building_blocks\BB1.csv``, you can just set the DELi
Data Directory to ``\home\my_user\del\data`` and list ``BB1``.
DELi will handle locating that file for you.
If a file is not found in the DELi Data Directory, DELi will raise a
``DELiFileNotFoundError`` to let you know which one is missing and why it
needed to look for it.

.. note::
    DELi can also handle file paths over names if preferred. It can also handle
    URIs to external resources (like S3 buckets or web URLs).

Setting the DELi Data Directory
-------------------------------
DELi will search for a path to the DELi Data Directory following this hierarchy:

* the value passed to the ``deli.configure.set_deli_data_dir`` function
* the value passed to the a CLI command using the ``--deli-data-dir`` option
* the value of the ``DELI_DATA_DIR`` environment variable
* the value of DELI_DATA_DIR in the .deli config file in your home directory
* the current working directory (if all other options are not set)

DELi will *not* continue to look for valid DELi Data Directory paths after
it finds any of the higher priority options to be non-null. For example, if you
set the DELI_DATA_DIR environment variable to ``not/a/valid/data/dir``, but you
are running a DELi command in a directory that *is* a valid DELi Data Directory,
DELi will still fail, as it will check the ``$DELI_DATA_DIR`` env first, see that
it is a non-null value, and then check if that is a valid DELi Data Directory. If
it is not, it will raise an exception. This is to avoid confusion about which DELi
Data Directory DELi is using.

How to Best Use the Data Directory
----------------------------------
Tracking DEL data files can be tricky, especially when working in a team.
The best practice for the DELi Data Directory is to have a single, version controlled
repository that holds your teams DEL info. This way, it can be pull from or copied
to local machines or containers or anything else in a more reproducible way.
At UNC, we use a private GitHub repository for this.
