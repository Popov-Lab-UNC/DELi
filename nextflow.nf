#!/usr/bin/env nextflow

params.decode_run
params.out_dir = "${launchDir}"
params.prefix = ""
params.debug = false
params.chunk_size = 50
params.save_failed = false

// this determines the selection id from the decode file
process GetSelectionIDPrefix {
    input:
    path decode_file

    output:
    val prefix into prefix_channel

    script:
    """
    grep '^selection_id:' $decode_file | sed 's/^selection_id: "\(.*\)"/\1/' | sed 's/ /_/g' | awk -v prefix="${params.prefix}" 'prefix != "" {print prefix "_" $0} prefix == "" {print $0}'
    """
}

process Decode {
    publishDir "$params.out_dir/logs/", mode: 'move'

    input:
    path fastq_file

    output:
    path '*_cube.csv', emit: cubes
    path '*_decode_statistics.json', emit: decode_stats
    path '*.log', emit: log

    script:
    """
    export sub_job_id=`basename ${params.decode_file} | sed 's/\.[^.]*$//'`
    deli decode run $decode_file $fastq_files --ignore-decode-seqs --skip-report ${params.debug ? '--debug' : ''} \
    ${params.save_failed ? '--save-failed' : ''} --prefix \$sub_job_id
    mv deli.log "deli.\$sub_job_id.log"
    """
}

process MergeCubes {
    publishDir "$params.out_dir/", mode: 'move'

    input:
    path "*_cube.csv"
    val prefix from prefix_channel

    output:
    path "${prefix}_cube.csv"

    script:
    """
    awk 'FNR==1 && NR!=1{next;}{print}' *_calls.csv > ${prefix}_cube.csv
    """
}

process MergeStats {
    publishDir "$params.out_dir/", mode: 'move'

    input:
    path '*_decode_statistics.json'
    val prefix from prefix_channel

    output:
    path "${prefix}_decode_statistics.json"
    path "${prefix}_decode_report.html"

    script:
    """
    deli decode statistics merge *_report_stats.json --out-path ${prefix}_decode_statistics.json
    deli decode report generate ${params.decode_run} ${prefix}_decode_statistics.json --out-path ${prefix}_decode_report.html
    """
}

process ExtractSequenceFiles {
    input:
    path decode_run

    output:
    path "sequence_files.txt"

    script:
    """
    #!/usr/bin/env python

    import yaml
    with open('${decode_run}', 'r') as f:
        data = yaml.safe_load(f)
    sequence_files = data.get('sequence_files', [])
    with open('sequence_files.txt', 'w') as out_f:
        out_f.write('\\n'.join(sequence_files))
    "
    """
}

workflow {
    decode_run = Channel.fromPath( params.decode_run ).first()
    GetSelectionIDPrefix( decode_run )

    sequence_files = ExtractSequenceFiles( decode_run )
        .splitText()
        .map { file(it.trim() }

    chunk_sequence_files = sequence_files.splitFastq( by: params.chunk_size, file: true )

    Decode(chunk_sequence_files)
    MergeCalls(Decode.out.cubes.collect())
    MergeReport(Decode.out.decode_stats.collect())
}
