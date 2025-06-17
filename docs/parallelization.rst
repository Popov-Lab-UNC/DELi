=======================
Parallelization in DELi
=======================

DELi as a package does not implement any type of parallelization.
It is design to run on a single CPU core.
This is in part because most analysis methods are quick enough that no
parallelization is needed.

The only place where parallelization could be useful is in enumeration and decoding.
Yet, the ``mp`` module from python is not very great, and both these
processes are `embarrassingly parallel <https://en.wikipedia.org/wiki/Embarrassingly_parallel>`_.
This means it's better to just run separate jobs on your system than to let python do it for you.

There are a few hiccups to overcome. For example in decoding you need to break up you fastq file
into chunk, and at the end the UMI corrected counters need to all be loaded at once to get the
right counts (others the same UMIs could be double counted). For enumeration you need to
figure out how to split the library building block up into subsets so each job enumerates only
part of the library.

These are not always trivial issues, but also not very hard to write a script to overcome.
To help show how this works, DELi has some example scripts using `Nextflow <https://www.nextflow.io/>`_
to parallelize decoding. The DELi team will continue to provide example Nextflow scripts for
future tasks as well, as Nextflow works well with HPC system (which most academic labs are
limited too) *and* a cloud environment since it is file based.

Running parallel decoding with Nextflow
---------------------------------------
You can run the nextflow script very similar the the decoding CLI command in DELi.
The largest difference is you can provide a chunk size which will determine how many jobs
are needed.

.. code-block:: shell
    nextflow run decode.nf ./example_decode.yaml --chunk_size 100000

One thing to keep in mind is that since Nextflow is file based, copies of the FASTQ files
will be made. If these are large, you can need alot of disk space.
This workflow also requires the raw counts for each DEL compound and their UMIs are saved
and then merged at the end (all loaded into memory). This step would be required even if
DELi handled parallelization internally, but with Nextflow we can limit the large memory
usage to the end, limiting how long we need the big memory VM for (and saving cost).
