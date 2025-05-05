#!/usr/bin/env nextflow

params.decode_run
params.fastq_file = null
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
    stdout

    script:
    """
    #!/usr/bin/env python

    import yaml
    with open('${params.decode_run}', 'r') as f:
        data = yaml.safe_load(f)
    selection_id = data.get('selection_id', None)
    prefix = "${params.prefix}" + "_" + selection_id if selection_id else "Unknown"
    if prefix.startswith("_"):
        prefix = prefix[1:]
    print(prefix)
    """
}

process Decode {
    publishDir "$params.out_dir/logs/", mode: 'move', pattern: "*.log"

    input:
    path fastq_file

    output:
    path '*_counters.pkl', emit: counters
    path '*_decode_statistics.json', emit: decode_stats
    path '*.log', emit: log

    script:
    """
    export sub_job_id=`basename ${fastq_file}`
    deli decode run ${params.decode_run} $fastq_file --ignore-decode-seqs --skip-report ${params.debug ? '--debug ' : ''}${params.save_failed ? '--save-failed ' : ''}--prefix \$sub_job_id --save-degen
    mv deli.log "deli.\$sub_job_id.log"
    """
}


process MergeCounters {
    input:
    path "*_counters.pkl"

    output:
    path "merged_counters.pkl"

    script:
    """
    deli decode counter merge *_counters.pkl --out-path merged_counters.pkl
    """
}


process MakeCube {
    publishDir "$params.out_dir/", mode: 'move'

    input:
    path "merged_counters.pkl"
    val prefix

    output:
    path "${prefix}_cube.csv"

    script:
    """
    deli cube from-counter merged_counters.pkl ${params.decode_run} --out-path ${prefix}_cube.csv
    """
}

process MergeStats {
    publishDir "$params.out_dir/", mode: 'move'

    input:
    path '*_decode_statistics.json'
    path "merged_counters.pkl"
    val prefix

    output:
    path "${prefix}_decode_statistics.json"
    path "${prefix}_decode_report.html"

    script:
    """
    deli decode statistics merge *_decode_statistics.json --counter-file merged_counters.pkl --out-path ${prefix}_decode_statistics.json
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
    """
}

workflow {
    DecodeRunChannel = Channel.fromPath( params.decode_run ).first()
    PrefixChannel = GetSelectionIDPrefix( DecodeRunChannel ).map( {it.trim()} )

    if ( params.fastq_file ) {
        chunk_sequence_files = Channel.fromPath( params.fastq_file ).splitFastq( by: params.chunk_size, file: true )
    }
    else {
        sequence_files = ExtractSequenceFiles( DecodeRunChannel ).splitText().map{ file(it.trim()) }
        chunk_sequence_files = sequence_files.splitFastq( by: params.chunk_size, file: true )
    }
    Decode(chunk_sequence_files)

    CounterChannel = MergeCounters(Decode.out.counters.collect())
    MakeCube(CounterChannel, PrefixChannel)
    MergeStats(Decode.out.decode_stats.collect(), CounterChannel, PrefixChannel)
}
