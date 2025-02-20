.. _deli-data-dir-ref:

DELi Data Directory
===================
To help organize DEL information, DELi searches for DEL components in a centralized location.
This is called the DELi Data Directory. You can specify that path to this directory in the the
.deli config file (See :ref:`deli-config-file-ref` for more info) or by setting the
``DELI_DATA_DIR`` environment variable to the path


Setting the DELi Data Directory
-------------------------------
DELi will search for a path to the DELi Data Directory following this hierarchy:

* the value of DELI_DATA_DIR in the .deli config file passed to the command (if applicable)
* ``DELI_DATA_DIR`` environment variable
* the value of DELI_DATA_DIR in the .deli config file in your home directory
* defaults to ``$HOME/.deli/deli_data``

You can also set the DELi data directory uisng the ``set_deli_data_dir`` function
in your code

Directory Structure
-------------------
There are 5 subfolders in the DELi Data directory:

* libraries
* barcodes
* building_blocks
* hamming
* indexes

As their names suggest, each sub folder should hold the definitions files
(See :ref:`define-dels`) for the components they represent. So all library json definitions
should go into the the libraries sub folder in the DELi Data Directory folder.

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
