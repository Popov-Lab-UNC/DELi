====================
DELi Documentation
====================

Welcome to DELi's documentation!
================================

DELi is a Python library for working with DNA-encoded libraries (DELs). It requires Python 3.13 or newer.

You can install DELi using pip:
::
    pip install deli-chem

Optional GNN analysis dependencies: ``pip install 'deli-chem[ml]'``

End-to-end example workflows are in the `examples <https://github.com/Popov-Lab-UNC/DELi/tree/main/examples>`_ directory.


.. toctree::
   :maxdepth: 2
   :caption: Introduction

   deli_glossary
   compound_ids

.. toctree::
   :maxdepth: 2
   :caption: Configuration

   config_docs/deli_config
   config_docs/deli_data_dir

.. toctree::
   :maxdepth: 2
   :caption: Defining DELs

   selection_files
   defining_docs/define

.. toctree::
   :maxdepth: 2
   :caption: Running DELi

   cli_docs
   parallelization

.. toctree::
   :maxdepth: 2
   :caption: Enumeration

   enumerator

.. toctree::
   :maxdepth: 2
   :caption: Decoding

   decoding_docs/decode

.. toctree::
   :maxdepth: 2
   :caption: Analysis

   analysis_docs/analysis_readme
   analysis_docs/analysis_config
   analysis_docs/analysis_contribute

.. toctree::
   :maxdepth: 2
   :caption: Development

   development_docs/testing
   development_docs/faqs

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/deli

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
