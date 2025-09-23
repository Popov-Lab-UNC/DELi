# DELi
![DELi](./DELi_logo.png)

> [!WARNING]
> **Disclaimer:** DELi is currently under active development. Features and documentation WILL change. Use with caution.

DELi (DNA-Encoded-Library informatics) is the software suite used by the CICBDD to do their in-house DEL analysis.
It incorporates the whole pipeline post base-calling/sequencing including:
- Barcode/DEL ID calling and cube file generation
- Enumeration of chemical structures from building blocks
- Disython and Monosynthon analysis
- Binary classification of enriched DELs (TODO)
- Generation of machine learning datasets and baseline models from DEL data
- Various digestible reports to understand the DEL results

You can read the detailed documentation [here](https://dna-encoded-library-informatics-deli.readthedocs.io/en/latest/).

## Installing DELi
You can install DELi using pip for any OS/Machine that supports Python 3.11+:
```shell
pip install del-informatics
```

## Why not a compiled language
DELi is written in Python for two reasons:
1. We wrote the first versions of it in Python
2. Python is the language most scientists in our field know, so it makes contributions from other DEL experts easier

It is true that DELi would be faster as a compiled C++ or Rust program, but we have optimized the DELi enough that runtime isn't much of an issue.
We hope to someday write a Rust version of DELi (at least for decoding and enumeration) but those plans are not yet in motion.
