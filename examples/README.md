# DELi examples
Here are some example for how to use DELi for decoding and analysis.
NOTE: these are user provided, so they are not guaranteed to be up to date with current versions.

## Background

**UNCDEL003** is a methyl-lysine reader focused DNA encoded library of around 60k members. This particular data provided was from a selection performed in triplicate with no NTC or naive, conducted on an Illumina sequencer. For more information, see this publication:
10.1021/acs.jmedchem.3c01192

To ensure your analysis is running correctly, you can check for validated features: in UNCDEL003, Feature C006 (methyl piperazine) has been experimentally validated.

**UNCDEL006** is an open-sourced benzamidazole-core small diversity library (~1M) featured in our DELi paper with a complete synthesis scheme and building blocks provided. This library includes chemistry that mimics prior known chemical matter for BRD4: 10.1021/acs.jmedchem.9b01670

For benchmarking your results, note that in UNCDEL006 BRD4, the library member A035-B040-C030 is a nanomolar (nM) binder with experimental validation.

## Enumerate

See `UNCDEL006_Enumerate` for an example of running enumeration using DELi. The folder includes three building block sets (A, B, C) with defined IDs and SMILES. The `UNCDEL006_Enumerate/enumerate.py` script is an example using DELi code to create building block sets and a reaction workflow for the enumeration of the library.

Running the example is as simple as:

```bash
python UNCDEL006_Enumerate/enumerate.py
```

## Decode

## Analyze

The `UNCDEL006_BRD4.csv.gz` file in the `UNCDEL006_BRD4_Analyze/` directory is output from the DELi decode module, containing decoded DEL sequences with enrichment data from BRD4 binding experiments.

To run analysis on this data:

```bash
# First unzip the data file
gunzip UNCDEL006_BRD4_Analyze/UNCDEL006_BRD4.csv.gz

# Run analysis
deli analyze --config UNCDEL006_BRD4_Analyze/analysis_config_del6.yaml
```

This will generate a comprehensive analysis report with plots, statistics, and chemical space visualizations.