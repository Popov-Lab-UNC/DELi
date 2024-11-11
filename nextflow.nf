#!/usr/bin/env nextflow

params.fastq_file
params.experiment
params.out_dir = $PWD
params.prefix = "deli_decode"
params.debug = false

process Decode {
    publishDir '$params.out_dir/logs/', mode: 'move', pattern: "*.log"

    input:
    path fastq
    path exp

    output:
    path '*_calls.csv', emit: calls
    path '*.log', emit: logs
    path '*_seq_lengths.json', emit: seq_lengths
    path '*_report_stats.txt', emit: decode_stats

    script:
    """
    deli decode $fastq $exp --save_report_data --skip_report
    """
}

process MergeCalls {
    publishDir '$params.out_dir/', mode: 'move'

    input:
    path "*_calls.csv"

    output:
    path "${params.prefix}_all_calls.csv"

    script:
    """
    awk 'FNR==1 && NR!=1{next;}{print}' *_calls.csv > ${params.prefix}_all_calls.csv
    """
}

process MergeReport {

}

workflow {

}
