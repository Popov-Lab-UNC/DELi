#!/usr/bin/env nextflow

params.decode_run
params.out_dir = "${launchDir}"
params.prefix = ""
params.debug = false
params.chunk_size = 1000000
params.save_failed = false
params.enumerate_smiles = false
params.include_bb_smi = false
params.save_counter = false


// this determines the selection id from the decode file
process GetSelectionIDPrefix {
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


process ValidateInputParams {
    output:
    val true, emit: has_umi // Placeholder, replace with actual validation logic
    val "pass", emit: exit_condition // Placeholder, replace with actual exit condition logic

    script:
    """
    echo "TODO"

    """
}


process CheckValidationExitCondition {
    input:
    val exit_condition

    when:
    exit_condition != "pass"

    script:
    """
    # Custom error logic here
    if [[ ${exit_condition} == "missing_reactions" ]]; then
        echo "Critical error: Requested SMILES enumeration but some libraries lack reaction information"
        exit 1
    elif [[ ${exit_condition} == "missing_bb_smiles" ]]; then
        echo "Critical error: Building block SMILES required but some libraries lack building block SMILES"
        exit 1
    fi
    """
}


process Decode {
    publishDir "$params.out_dir/decode_logs/", mode: 'move', pattern: "*.log"

    input:
    path fastq_file

    output:
    path '*_counter_subjob.json.gz', emit: counters
    path '*_decode_statistics_subjob.json', emit: decode_stats
    path '*.log', emit: log
    path '*_decode_failed.tsv', emit: failed, optional: true

    script:
    """
    python ${workflow.projectDir}/decode.py --decode_file ${params.decode_run} --fastq_file ${fastq_file} ${params.debug ? '--debug ' : ''}${params.save_failed ? '--save-failed ' : ''}
    """
}

process MergeResultsUMI {
    publishDir "$params.out_dir/", mode: 'move'

    input:
    path "*_counter_subjob.json.gz"
    path "*_decode_statistics_subjob.json"

    output:
    path "*_cube.csv"
    path "*_decode_report.html"
    path "*_decode_statistics.json"
    path "*_counter.json", optional: true

    script:
    """
    python ${workflow.projectDir}/merge_umi.py --decode_file ${params.decode_run} ${params.enumerate_smiles ? '--enumerate_smiles ' : ''}${params.include_bb_smi ? '--include_bb_smi ' : ''}${params.save_counter ? '--save_counter ' : ''}--counters *_counter_subjob.json.gz
    """
}


process MergeFailedDecodes {
    publishDir "$params.out_dir/", mode: 'move'

    input:
    path "*_decode_failed.tsv"
    val prefix

    output:
    path "*_failed_decodes.tsv"

    script:
    """
    awk 'FNR==1 && NR!=1{next;}{print}' *_decode_failed.tsv > ${prefix}_failed_decodes.tsv
    """
}


process ExtractSequenceFiles {
    output:
    path "sequence_files.txt"

    script:
    """
    #!/usr/bin/env python

    import yaml
    with open('${params.decode_run}', 'r') as f:
        data = yaml.safe_load(f)
    sequence_files = data.get('sequence_files', [])
    with open('sequence_files.txt', 'w') as out_f:
        out_f.write('\\n'.join(sequence_files))
    """
}

workflow {
    PrefixChannel = GetSelectionIDPrefix().map( {it.trim()} )

    ValidateInputParams()
    CheckValidationExitCondition(ValidateInputParams.out.exit_condition)

    sequence_files = ExtractSequenceFiles().splitText().map{ file(it.trim()) }
    chunk_sequence_files = sequence_files.splitFastq( by: params.chunk_size, file: true )

    Decode(chunk_sequence_files)

    if( ValidateInputParams.out.has_umi )
        MergeResultsUMI(Decode.out.counters.collect(), Decode.out.decode_stats.collect())
    else
        error "UMI information is required for this workflow. Please check your input parameters."

    if ( params.save_failed )
        MergeFailedDecodes(Decode.out.failed.collect(), PrefixChannel)

}
