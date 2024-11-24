#!/usr/bin/env nextflow

params.fastq_file
params.experiment
params.out_dir = "${launchDir}"
params.prefix = "deli_test"
params.debug = false
params.chuck_size = 50

process Decode {
    publishDir "$params.out_dir/logs/", mode: 'move', pattern: "*.log"

    input:
    path fastq
    path exp

    output:
    path '*_calls.csv', emit: calls
    path '*.log', emit: log
    path '*_report_stats.json', emit: decode_stats

    script:
    """
    deli decode $fastq $exp --save_report_data --skip_report
    export sub_job_id=`echo $fastq | awk -F'.' '{print \$2}'`
    mv deli.log "deli.\$sub_job_id.log"
    """
}

process MergeCalls {
    publishDir "$params.out_dir/", mode: 'move'

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
    publishDir "$params.out_dir/", mode: 'move'

    input:
    path '*_report_stats.json'

    output:
    path "${params.prefix}_decode_report.html"

    script:
    """
    deli report merge *_report_stats.json --render_report --name ${params.prefix}_decode_report
    """
}

workflow {
    fastq_files = Channel.fromPath(params.fastq_file).splitFastq(by: params.chuck_size, file: true)
    experiment = Channel.fromPath(params.experiment).first()
    Decode(fastq_files, experiment)
    MergeCalls(Decode.out.calls.collect())
    MergeReport(Decode.out.decode_stats.collect())
}
