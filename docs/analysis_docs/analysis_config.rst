Analysis Configuration
=======================

General Settings
----------------
- **data**: Path to the data file.
    - Example: ``analysis/example_data/DEL003_cube_noid.csv``
- **lib_size**: Library size.
    - Example: ``58000``
- **output_dir**: Output directory.
    - Example: ``DEL3``

ID Column
---------
- **id_col**: Identifier column.
    - Example: ``DEL_ID``

Indexes
-------
- **indexes**: Columns with enrichment values (ideally UMI corrected), for your DEL members.
    - **Cycle1**: ``["corrected_PTMODD3index4", "corrected_PTMODD3index5", "corrected_PTMODD3index6"]``
    - **Cycle2**: ``["corrected_PTMODD3index1", "corrected_PTMODD3index2", "corrected_PTMODD3index3"]``

Control Columns
---------------
- **control_cols**: NTC control enrichment columns for your given experiments.
    - **Cycle1**: ``["corrected_D3index7"]``
    - **Cycle2**: ``["corrected_D3index8"]``

Raw Indexes
-----------
- **raw_indexes**: Raw or non-UMI corrected enrichment columns for each member. Indexes will default to using these values if UMI corrected values are not provided.
    - **Cycle1**: ``["raw_PTMODD3index4", "raw_PTMODD3index5", "raw_PTMODD3index6"]``
    - **Cycle2**: ``["raw_PTMODD3index1", "raw_PTMODD3index2", "raw_PTMODD3index3"]``

Flags
-----
- **SD_min**: Enable or disable Sampling Depth/Sampling Depth minimum threshold calculation.
    - Example: ``True``
- **NSC_values**: Enable or disable nornmalized sequence count calculation.
    - Example: ``True``
- **MLE**: Enable or disable maximum liklihood estimate calculation.
    - Example: ``True``
- **Z_score**: Enable or disable normalized Z_score calculation.
    - Example: ``True``
- **z_score_log_data**: Enable or disable log transformed normalized Z_score calculation.
    - Example: ``True``
- **normalized_data**: Enable or disable background subtraction using NTC calculation, should generally be performed after other operations that require usage of NTC as it will remove the NTC columns.
    - Example: ``False``
- **disynthon_data**: Enable or disable disynthon enrichment calculation.
    - Example: ``true``
- **trisynthon_overlap**: Enable or disable trisynthon reproducibility calculation which displays venn diagram of library members across replicate runs.
    - Example: ``True``
- **disynthon_overlap**: Enable or disable disynthon reproducibility calculation which displays venn diagram of library members across replicate runs.
    - Example: ``True``
- **disynthon_threshold**: Set the disynthon threshold value for overlap visualization.
    - Example: ``100``
- **simple_spotfire_version**: Enable or disable function to output simplified csv that we prefer for handing off to chemistry team for Spotfire/Data Warrior usage.
    - Example: ``true``
- **ml_fingerprints_to_RF_reg**: Enable or disable function that creates RF regressor using chemical fingerprints.
    - Example: ``true``
- **ml_fingerprints_to_RF_clf**: Enable or disable function that creates RF classifier using chemical fingerprints.
    - Example: ``true``
- **clf_thresh**: Set the threshold value for RF classifier.
    - Example: ``15``
- **gnn_classifier**: Enable or disable construction of a graph neural network classifier, this adds about 3 - 5 min of runtime on my macbook.
    - Example: ``false``
- **gnn_threshold**: Set the gnn threshold value.
    - Example: ``10``
- **gnn_arch**: Set the architecture you want to use either GAT or GIN.
    - Example: ``GAT``
- **top_hits**: Set the number of structure hits you want to pull for final visualizaition.
    - Example: ``20``
- **top_hits_metric**: Set the metric you want to use to pull top structures.
    - Example: ``sum``
- **report**: Enable or disable self-contained html report generation.
    - Example: ``true``

Top Disynthons
--------------
- **top_disynthons**: Configuration for top disynthons.
    - **enabled**: Enable or disable the comparison of top disynthons, regularly used to compare to NTC or an inhibitor binding experiment.
        - Example: ``True``
    - **comparison**: Type of comparison.
        - Example: ``"exp2"`` (Options: "control", "exp2", "none")
    - **exp_name**: Experiment name, in this case would be Cycle2.
        - Example: ``"Cycle2"``
    - **exp2_name**: Second experiment name, in this case would be Cycle 1.
        - Example: ``"Cycle1"``
    - **control_name**: Control name, if i chose control I'd add in Cycle2 here to compare my Cycle2 to the Cycle2 NTC.
        - Example: ``"none"``
    - **top_count**: Number of top disynthons to pull, these will be labeled and plotted for each disynthon pair.
        - Example: ``10``
    - **comparison_metric**: Metric for comparison.
        - Example: ``"avg"``
