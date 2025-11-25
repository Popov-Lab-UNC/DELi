.. _decoding-docs:

============
DEL Decoding
============

DELi support decoding of raw DEL selection data (DNA sequences, or "reads") into DEL compound enrichment data (what
scientist and the analysis module use to find "hits"). This what the ``decode`` module in DELi handles.

The Basics
----------
DELi assumes a given decoding effort is for a single experimental selection. It does not have the ability to
demultiplex sequence files into multiple selections, nor decode more than one selection at a time. You will need
to launch a separate decoding run for each selection you want to decode.

.. note::
    The DELi roadmap includes plans to support "selection_id" tags which will essentially enable demultiplexing
    selections, as it will route calls to separate files for each unique selection tag. Raise an issue if you
    would like to see this feature sooner rather than later.

Run Decoding
------------
DELi supplies both a command line interface (CLI) and a Python API for decoding DELs.
For most users, the CLI is the easiest way to decode DELs.
All you need to do is provide a :ref:`decoding run file <decoding-run-file-docs>` and DELi will
start decoding.

If you want to decode via the Python API you can use the ``DecodingRunner`` class. You can
get even more control over the decoding process by using if you want. Refer to the decoding
module API documentation for more details.


.. toctree::
    :maxdepth: 2

    algorithm
    barcode_calling
    decoding_outcomes
    tool_compounds
    decoding_output
    decoding_file
