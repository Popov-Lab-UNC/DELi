.. _analysis-config-docs:

Analysis configuration (YAML)
==============================

The automated analysis driver reads a single YAML file. A full skeleton with multiple experiments, comparisons, and optional **raw_indexes** is in the repository at ``examples/analysis_examples/DEL005_cube_config_sample.yaml``. How to run the driver is described in :doc:`analysis_docs/analysis_readme`; this page only documents the config file itself.

General section
---------------

All of these keys live under **general**:

- **data**: Path to the cube CSV (path is resolved relative to the process working directory unless you pass an absolute path).
- **lib_size**: Encoded library size (integer). Required for sampling-depth metrics (``SD_min``, ``NSC_values``), PolyO, and the normalized z-score path when enabled.
- **output_dir**: Base directory for outputs; the driver also creates a dated subdirectory for each run.
- **ID_col**: Column used as the compound identifier (for example ``DEL_ID``). If omitted, the first column of the CSV is used.

Indexes and controls
--------------------

- **indexes**: Maps each **experiment name** (a label you choose) to a list of count columns for that experiment—typically UMI-corrected or otherwise decodable read counts per replicate or channel.
- **control_cols**: Maps each experiment name to the NTC (or other control) column(s) used by metrics that compare selection to control—for example MLE and log-based z-scores. Keys should match the experiment names you use in **indexes** where a control exists.
- **raw_indexes**: Optional. If **indexes** is empty, **raw_indexes** is copied to **indexes**. For current implementations, NSC, ``SD_min``, and PolyO use the columns listed under **indexes** (not a separate raw-only path). Keeping **raw_indexes** in the file is optional and mainly useful for bookkeeping or future use.

Flags block
-----------

All pipeline switches sit under **flags**. They are booleans unless noted. The driver runs steps in a fixed order; some flags depend on others (for example disynthon-based plots require **disynthon_data** first).

Understanding the flags (grouped)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use the following groupings when reading an example config (such as the multi-target sample linked above):

**Enrichment and QC metrics**

- **SD_min** / **NSC_values**: Sampling depth and normalized sequence counts from **indexes** and **lib_size**.
- **MLE**: Maximum-likelihood enrichment style ratios using **indexes** and **control_cols**.
- **Z_score**: Normalized z-scores under a uniform-library null (expected fraction ``1 / lib_size``), using **indexes**; requires **lib_size**.
- **z_score_log_data**: Log-space z-style scores using experiment and control columns from **control_cols**; independent of **Z_score**’s theoretical null.

**Disynthon / PolyO / overlap**

- **disynthon_data**: Builds disynthon-level tables; required before PolyO and disynthon overlap steps that consume them.
- **polyO**: PolyO scores on the cube (uses **indexes** and **lib_size**).
- **trisynthon_overlap** / **disynthon_overlap**: Reproducibility-style plots at tri- and disynthon resolution.
- **disynthon_threshold**: Numeric cutoff for overlap filtering, or the string ``"auto"`` where supported.

**Normalization**

- **normalized_data**: Subtracts control counts from selection columns (per experiment). Run after steps that still need raw control columns if you want normalized exports downstream.

**Downstream CSVs and ML (optional)**

- **simple_spotfire_version**: Writes a reduced column set for Spotfire or similar tools.
- **ml_fingerprints_to_RF_reg** / **ml_fingerprints_to_RF_clf**: Random forest models from fingerprints; **clf_thresh** sets the classification cutoff.
- **gnn_classifier**, **gnn_threshold**, **gnn_arch**: Optional graph neural network classifier (architecture string such as ``GAT`` or ``GIN``).

**Hits and report**

- **top_hits** / **top_hits_metric**: How many compounds to list and whether to rank by **sum** or **avg** over index columns.
- **report**: Writes the bundled HTML report and a dated cube CSV next to the run.

Nested flag structures
^^^^^^^^^^^^^^^^^^^^^^

**top_disynthons**

A **comparisons** list; each item is one chart batch. Fields include **comparison** (``"control"``, ``"exp2"``, or ``"none"``), **exp_name**, optional **exp2_name**, **control_name** where required by the comparison type, **top_count**, and **comparison_metric** (``sum`` or ``avg``). Multiple rows mean multiple comparison sets.

**top_delta_compounds**

A **comparisons** list of **exp_name** / **control_name** pairs (keys into **indexes** / **control_cols**), plus **metric** and **top_count**.

**monosynthon_chemical_space**

Either enabled as a boolean (all experiments) or as a mapping with **experiments**: a list of experiment names to include in the t-SNE-style view.

**report** is typically ``true`` when you want the HTML summary; it collects outputs from whatever earlier flags actually ran.
