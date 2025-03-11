# DELi
![DELi](./DELi_logo.png)

> [!WARNING]
> **Disclaimer:** DELi is currently under active development. Features and documentation WILL change. Use with caution.

DELi (DNA-Encoded-Library informatics) is the software suite used by the CICBDD to do their in-house DEL analysis.
It incorporates the whole pipeline post base-calling/sequencing including
- Barcode/DEL ID calling and cube file generation
- Mapping DEL IDs to enumerated chemical structures 
- Disython and Monosynthon analysis
- Binary classification of enriched DELs (TODO)
- Generation of machine learning datasets and baseline models from DEL data (TODO)
- various digestible reports to understand the DEL results (TODO)

## Why not a compiled language
DELi is written in python for 2 reason:
1. We wrote the first versions of it in python (we are grad students not software devs)
2. Python is the langauge most scientist in our feild know, so it makes contributions from other DEL experts easier

It is true that DELi would likly be faster as a compiled C++ or Rust program, but we utilize Numba to help 
keep things fast enough to be usable for most academic settings. DELi can handle calling around 10 million 
reads in 1-2 hours on a (built in 2022) 32 core machine. If you have billions of reads, you DELi will need
more cores to have good runtimes

We are interested in writing the calling code of DELi into a faster language, probably Rust cause it sounds
like a fun language to learn. But we (the orginal authors) are too bust trying to get PhDs to spend enough 
time to write something good. If your are interested in do it for us, please reach out.

## Installing DELi
To install DELi,
first clone the repo and then set up a virtual environment
(conda or venv, both will work) by installing the `requirements.txt` file:
```shell
git clone https://github.com/Popov-Lab-UNC/DELi.git
cd DELi
conda create -n deli
conda activate deli
conda install pip
pip install -r requirements.txt
```
NOTE: DELi is built to run on linux or windows,
it has not been tested on MacOS, but in theory should be work there as well

DELi is built will multiprocessing as the matching and barcode calling is a CPU hungry process.
It is recommended that you install DELi on a computer with at least 20 cores.
Currently, the Popov lab uses a machine with a 32 core Threadripper (`glassfrog.dhcp.unc.edu`)

## Running DELi
While DELi can be run in discrete steps,
it is recommended that you run the entire pipeline in a single call to `deli_run.py`.
This command will run through all step of converting a `.fastq` file
(post base-calling, see [dorado](#dorado-setup-for-base-calling-))
into the cube file that can be analyzed in DataWarrior or Spotfire.

Details about the settings can be obtained with `python3 deli_run.py -h`
```
positional arguments:
  input                 input fastq file with reads to be searched for matches

options:
  -h, --help            show this help message and exit
  --make                the 'make' ID defining the makeup of the library
  -i INDEX [INDEX ...], --index INDEX [INDEX ...]
                        index(s) included in the DEL selection being called; pass "all" for all indexes
  -l LIBRARY [LIBRARY ...], --library LIBRARY [LIBRARY ...]
                        library(s) included in the DEL selection being called; pass "all" for all library(s)
  -t CONTROL, --control CONTROL
                        index of selection that is the negative control
  -o OUTDIR, --outdir OUTDIR
                        directory to write the output files to
  -s, --save-match      whether to save the matches as a file
  -r, --recursive       recursively loop through all file in a directory
  -m {full,minimum,double,single}, --mode {full,minimum,double,single}
                        which barcode query to use when searching for matches, see readme.md for more details
  -e ERROR, --error ERROR
                        how many edit distance errors you are willing to allow when attempting a match
  -c, --compliment      also search for the reverse compliment of the DEL
  -p PREFIX, --prefix PREFIX
                        prefix to add to the filename of output called csv file(s)
  -u, --umi-cluster     conduct a UMI clustering to further denoise results
  --strict              use a more strict alignment and pattern match that uses the closing primer post UMI
  -n, --normalize       normalize the DEL_ID sequence counts by the control (requires -t/--control is set)
  --monosynthon         conduct a monosynthon analysis on all generated cubes
  --disynthon           conduct a monosynthon analysis on all generated cubes
  -j N_WORKERS, --n-workers N_WORKERS
                        number of worked to use for calling
  --index_name          path to json file containing a dictionary mapping passed index to selection names
  --debug               enable debug logging
  --print               print updating status to console
```

While the settings allow for a fair amount of customization,
it is recommended that you use the default settings from the example below to get the best results.
This example assumes that we are calling a run that includes libraries 'DEL004' and 'DEL005'
(see the [library json](https://github.com/Popov-Lab-UNC/CICBDD-DEL-Pipeline/blob/main/data/libs.json) for list of libraries),
and indexes 'index5', 'index6', and 'index10'
(see the [index json](https://github.com/Popov-Lab-UNC/CICBDD-DEL-Pipeline/blob/main/data/experiment_index.json) for list of indexes) with 'index10' being the control (bead only) selection
```shell
python3 deli_run.py <PATH/TO/MY/FASTQ> -i index5 index6 index10 -l DEL004 DEL005 -o <PATH/TO/OUTPUT/DIR> -s -m single -e 3 -c -p <MY_PREFIX> -n -t index10 --disynthon --debug --print --index_name <PATH/TO/INDEX/JSON>
```

This will produce the desired cube files in the `output-dir`.
There will be a single file for each library, each containing all indexes.
There will also be a file for each index, each containing all libraries.
