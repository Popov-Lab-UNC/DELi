# DELi
![DELi](./assets/DELi_logo.png)

DELi (DNA-Encoded Library informatics) is a Python library for working with DELs.
It incorporates the whole pipeline post base-calling/sequencing including:
1. Barcode/DEL ID calling and cube file generation
2. Enumeration of chemical structures from building blocks
3. Disython and Monosynthon analysis
4. Generation of machine learning datasets and baseline models from DEL data
5. Various digestible reports to understand the DEL results

You can read the detailed documentation [here](https://dna-encoded-library-informatics-deli.readthedocs.io/en/latest/).

## Installing DELi
You can install DELi using pip for any OS/Machine that supports Python 3.10+:

```shell
pip install deli-chem
```

## Getting Started

You can use DELi as a command line tool (see the [docs](https://dna-encoded-library-informatics-deli.readthedocs.io/en/latest/cli_docs.html) for more details) or as a python package
```python
import deli
print(deli.__version__)
```

For an end-to-end workflow of running DELi with open source libraries and selections (Enumerate, Decode, Analyze), see the [examples documentation](examples/README.md).

## Why not a compiled language
DELi is written in Python for two reasons:
1. We wrote the first versions of it in Python
2. Python is the language most scientists in our field know, so it makes contributions from other DEL experts easier

It is true that DELi would be faster as a compiled C++ or Rust program, but we have optimized the DELi enough that runtime isn't much of an issue.
We hope to someday write a Rust version of DELi (at least for decoding and enumeration) but those plans are not yet in motion.

**Note for developers:** DELi is built using poetry. You can use `poetry build` to build from source after
cloning the repo. Be on the lookup for contribution docs in the near future!
