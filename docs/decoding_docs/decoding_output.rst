.. _decoding-output-docs:

===================
DEL Decoding Output
===================
After a decoding run, DELi can generate a couple of different outputs containing information about the decoding process and results.
The main two of these are the "cube" file and the "decoding statistics".

Cube File
---------
The cube file is the main output of the decoding process. It contains the decoded compound IDs and their associated enrichment values.
It can also include additional information such as the enumerated SMILES, Library IDs, Building Block IDs, and Building Block SMILES.

Decoding Statistics
-------------------
The decoding statistics file is in JSON format and contains information about the decoding process, including the number of reads decoded,
the number of reads that failed to decode, and the reasons for any failures. This is is mostly information to help with
downstream enrichment analysis, as it contains information about the sequencing and read depth (which some enrichment analysis tools require).

Decoding Report
---------------
The decoding report is a summary of the decoding process rendered in HTML. It includes all the information from the decoding statistics file,
but in a far more user friendly format with graphs and tables.

Counter File
------------
In the case that you want to save the raw counts (including the UMIs sets) from the decoding process, DELi can also generate a counter file.
This is a JSON file that contains all that information in a format that can be easily read by other tools.
The JSON is formatted as a nested dictionary, where the first layer seperates compounds by library ID, and the second layer contains the
compound ID as the key and the value is another dictionary with the relevant count information for that compound (like raw count and UMIs).

This file is a very useful for parallelization, as UMI degeneration cannot be done without having all the raw UMIs (which the cube does not include).
This way, you can write a script that will read the counter JSON files for all parallel runs and correct the UMI counts to merge the results. You
can see an example of this in the ``nextflow`` example folder in the DELi repository.
