#!/usr/bin/env nextflow

params.fastq_files
params.experiment_setup

process Match {
    input:
    path fastq_file
    path experiment

    output:
    path '*.match'

    script:
    """
    deli match --input $fastq_file --experiment $experiment
    """
}

process Call {
    input:
    path match_file
    path experiment

    output:
    path '*.call.csv'

    script:
    """
    deli call --input $match_file --experiment $experiment
    """
}
