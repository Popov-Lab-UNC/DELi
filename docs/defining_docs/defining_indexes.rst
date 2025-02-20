Defining Indexes
================

Sometimes, you might need DELi to de-multiplex your sequences for you.
If that is the case, you will need to define the indexes to use.

This is very simple: indexes should be only one line of text in a ``.txt`` file.
That text should be the DNA tag linked to that index.
For example, my file index1.txt contains
::
    TGGAGCAGTG

The name of the file is the name of the index. So if I told a DELi experiment I was using
index1, it will look for the index file named "index1.txt"

Every unique index should have its own file

DELi Data Dir
=============
Index files can be saved in the DELi Data Dir sub-dir named "indexes" and end with ".txt" file extensions
