DELi Analysis Documentation
===========================

Overview
--------

DELi can be utilized in a variety of ways. In efforts to provide maximal flexibility to
end users, we provided an "automated" analysis mode and a "manual" analysis mode.

Automated Analysis Mode
------------------------

The automated method allows you to :ref:`configure a single YAML file <analysis-config-docs>`, choose the analyses you want to run, and will automatically generate a report, figures, and CSVs you can provide to your chemistry team.

The automated analysis mode is the default mode and can be run with the following command:

.. code-block:: python

    python analysis.py --config <path_to_config.yaml>

This will generate a report and all the figures and CSVs you need to provide to your chemistry team.

Manual Analysis Mode
---------------------

The manual method will be preferred for those that want more control and is likely a better entry point for those that opt out of DELi's matching modules but would still like access to the analysis modules.

A Jupyter notebook is provided to help guide you through the manual analysis process. We
intend to provide future cookbooks and tutorials for more advanced analyses and to showcase
some regular use cases we encounter with our DEL work at UNC.
