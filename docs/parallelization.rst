=======================
Parallelization in DELi
=======================

First, Parallelization in DELi is really only important for the decoding module (at this
point in time).
There, DELi support minimal parallelization in the package/code itself. That is in part because most
tasks DELi performs are embarrassingly parallel, and pythons native parallelization
is much less efficient that manually running separate DELi jobs. For example, if I have
10 Billion sequences to decode, attempting to use the `multiprocessing` module in python
to form a pool of 40 workers is far slower and prone to strange Python shenanigans than just
splitting the 10 Billion sequences into 40
files and running 40 DELi decode commands in parallel.

Instead, we recommend splitting you job up into smaller decoding tasks and merging the
results. DELi implements a handful of merging commands to help streamline this process.
This also makes DELi easy to use on distributed systems, like HPC clusters or cloud.

Nextflow and DELi
=================
While DELi does not support parallelization in the code itself, it can be used with
very easily with workflow managers. The Authors of DELi prefer Nextflow for this
since it can natively split FASTQ files. To help streamline this, an example Nextflow
configuration file and script is provided in the DELi repository. This file can be used to run
easily parallelized DELi decode jobs. It could also be expanded to add fun custom additions,
like emailing progress updates or pushing files to remote storage, features DELi will not
support internally.

Running Nextflow for parallel decoding
--------------------------------------
TODO
