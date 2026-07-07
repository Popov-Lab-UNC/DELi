# DELi examples

Example workflows for enumeration, decoding, and analysis using DELi.

## Background

**UNCDEL003** is a methyl-lysine reader focused DNA encoded library of around 60k members. This particular data provided was from a selection performed in triplicate with no NTC or naive, conducted on an Illumina sequencer. For more information, see this publication:
10.1021/acs.jmedchem.3c01192

To ensure your analysis is running correctly, you can check for validated features: in UNCDEL003, Feature C006 (methyl piperazine) has been experimentally validated.

**UNCDEL006** is an open-sourced benzamidazole-core small diversity library (~1M) featured in our DELi paper with a complete synthesis scheme and building blocks provided. This library includes chemistry that mimics prior known chemical matter for BRD4: 10.1021/acs.jmedchem.9b01670

For benchmarking your results, note that in UNCDEL006 BRD4, the library member A035-B040-C030 is a nanomolar (nM) binder with experimental validation.

## Enumerate

See `UNCDEL006_Enumerate` for an example of enumerating DEL006. The folder includes three building block CSV files (A, B, C) with IDs, SMILES, and barcodes. The `enumerate.py` script loads the library from `example_deli_data_dir` and enumerates a single compound.

```bash
cd UNCDEL006_Enumerate
python enumerate.py
```

You can also enumerate from the CLI using the library definition in `example_deli_data_dir`:

```bash
deli --deli-data-dir examples/example_deli_data_dir enumerate \
  examples/example_deli_data_dir/libraries/DEL006.json \
  -o DEL006_enumerated.csv
```

## Decode

Data for an example decoding run can be found in `UNCDEL006_BRD4_Decode`. Point DELi at the example data directory, then run decode:

```shell
cd UNCDEL006_BRD4_Decode

deli --deli-data-dir ../example_deli_data_dir decode run \
  -t -p EXAMPLE -o ./output -f example_decode.yaml
```

Output files are written to `./output` (decoded TSV, statistics JSON, HTML report, and failed reads).

## Analyze

The `UNCDEL006_BRD4.csv.gz` file in `UNCDEL006_BRD4_Analyze/` is example cube data from a DEL006 BRD4 selection, with enrichment columns for analysis.

```bash
cd UNCDEL006_BRD4_Analyze

gunzip -k UNCDEL006_BRD4.csv.gz   # skip if already decompressed
deli analyze --config analysis_config_del6.yaml
```

A smaller DEL003 example is in `UNCDEL003_53bp1_Analyze/`:

```bash
cd UNCDEL003_53bp1_Analyze
gunzip -k UNCDEL003_53bp1.csv.gz
deli analyze --config analysis_config_del3.yaml
```

These commands generate HTML analysis reports with plots, statistics, and chemical space visualizations.
