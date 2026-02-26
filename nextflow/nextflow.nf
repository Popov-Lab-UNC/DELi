#!/usr/bin/env nextflow

/*
 * DELi Decoding Workflow
 *
 * This workflow:
 * 1. Chunks FASTQ files from a selection configuration
 * 2. Runs decoding in parallel on all chunks with --split-by-lib
 * 3. Collects decoded sequences per library into NDJSON format
 * 4. Counts compounds per library asynchronously
 * 5. Optionally merges final results
 */

nextflow.enable.dsl = 2

// ============================================================================
// PARAMETERS
// ============================================================================

params.selection_file = null // Required: path to selection YAML file
params.out_dir = "${launchDir}/decode_results"
params.prefix = "" // Empty string will use selection_id
params.chunk_size = 1_000_000 // Split FASTQ at this many reads per chunk
params.deli_data_dir = null
params.config_file = null
params.debug = false

// ============================================================================
// VALIDATION & SETUP
// ============================================================================

if (!params.selection_file) {
    error("--selection_file is required")
}

selection_file_path = file(params.selection_file)
if (!selection_file_path.exists()) {
    error("Selection file not found: ${params.selection_file}")
}

// ============================================================================
// PROCESSES
// ============================================================================

process ExtractSequenceFiles {
    /*
     * Extract sequence files, selection_id, and libraries from selection YAML
     * Outputs:
     * - selection_id: The selection ID from the selection config
     * - libraries: Newline-separated list of library IDs
     * - files: Newline-separated list of sequence file paths
     */

    input:
    path selection_file

    output:
    path "selection_id.txt", emit: selection_id
    path "files.txt", emit: files

    script:
    """
    #!/usr/bin/env python
    import yaml

    with open("${selection_file}") as f:
        config = yaml.safe_load(f)

    # Write selection_id
    selection_id = config.get('selection_id', 'unknown')
    with open('selection_id.txt', 'w') as f:
        f.write(selection_id + '\\n')

    # Write sequence files
    sequence_files = config.get('sequence_files', [])
    with open('files.txt', 'w') as f:
        if isinstance(sequence_files, list):
            f.write('\\n'.join(sequence_files) + '\\n')
        else:
            f.write(sequence_files + '\\n')
    """
}

process DecodeChunk {
    /*
     * Run deli decode run on a chunk of FASTQ
     * Outputs one TSV file per chunk (no split-by-lib)
     */
    tag "${fastq_chunk.simpleName}"

    input:
    path fastq_chunk
    path selection_file
    val prefix
    val deli_args

    output:
    path "${prefix}_${fastq_chunk.simpleName}_decoded.tsv", emit: decoded_tsv
    path "${prefix}_${fastq_chunk.simpleName}_decode_statistics.json", emit: decode_stats
    path "deli.log", emit: deli_log
    """
    mkdir -p decoded_output

    deli ${deli_args} decode run \
        "${selection_file}" \
        "${fastq_chunk}" \
        --out-dir ./ \
        --prefix "${prefix}_${fastq_chunk.simpleName}" \
        --skip-report
    """
}

process MergeDecodeStatistics {
    /*
     * Merge decode statistics JSON files from all chunks into a single JSON file with aggregated stats
     */
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path("*_decode_statistics.json", arity: '1..*')
    path selection_file
    val prefix
    val deli_args

    output:
    path "${prefix}_decode_stats.json", emit: merged_stats

    """
    deli ${deli_args} decode merge-stats \
        *_decode_statistics.json \
        --selection-file "${selection_file}" \
        --out-loc "${prefix}_decode_stats.json"
    """
}

process CollectDecodeChunks {
    /*
     * Collect all chunk decode outputs files into a single collected NDJSON file
     */

    input:
    path("*_decoded.tsv", arity: '1..*')
    val prefix
    val deli_args

    output:
    path "${prefix}_collected.ndjson", emit: ndjson

    """
    # Run deli decode collect on all TSV files
    deli ${deli_args} decode collect \
        *_decoded.tsv \
        --out-loc "${prefix}_collected.ndjson"
    """
}

process CountChunk {
    /*
     * Count compounds from a chunk of the collected NDJSON file
     * Input file contains up to 500,000 NDJSON lines
     */
    tag "${ndjson_chunk.simplename}"

    input:
    path ("*.ndjson", arity: '1..*')
    val prefix
    val deli_args

    output:
    path "${ndjson_chunk.simplename}_counted.parquet", emit: counted

    script:
    """
    deli ${deli_args} decode count \
        "${ndjson_chunk}" \
        --out-loc "${ndjson_chunk.simplename}_counted.parquet" \
        --output-format parquet \
        --cluster-umis \
        --keep-raw-count \
        --keep-dedup-count
    """
}

process CollectCountChunks {
    /*
     * Merge counted files from all chunks
     */
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path counted_files
    val prefix

    output:
    path "${prefix}_counts.parquet", emit: merged_counts

    script:
    """
    #!/usr/bin/env python

    import polars as pl
    files = sorted([f for f in "${counted_files}".split() if f.strip()])
    pl.scan_parquet(files).sink_parquet("${prefix}_counts.parquet")
    """
}

process SummarizeDecodeRun {
    /*
     * Summarize the decode run by merging the decode statistics with the final counts to produce a final decode summary JSON file
     */
    publishDir "${params.out_dir}", mode: 'move'

    input:
    path merged_counts
    path decode_stats
    val prefix
    val deli_args

    output:
    path "${prefix}_final_stats.json", emit: final_stats

    script:
    """
    deli ${deli_args} decode summarize \
        "${merged_counts}" \
        "${decode_stats}" \
        --out-loc "${prefix}_decode_summary.json"
    """
}

process WriteDecodeReport {
    /*
     * Generate the decoding HTML report
     */

    input:
    path final_stats
    val prefix
    val deli_args

    output:
    path "${prefix}_decode_report.html", emit: report

    script:
    """
    deli ${deli_args} decode report \
        "${final_stats}" \
        --out-loc "${prefix}_decode_report.html"
    """
}

// ============================================================================
// WORKFLOW
// ============================================================================

workflow {
    // Validate potentially user-controlled path parameters to prevent shell injection
    def safePathPattern = ~/^[\w.\-\/]+$/
    if (params.deli_data_dir && !(params.deli_data_dir ==~ safePathPattern)) {
        error("Invalid characters in --deli_data_dir parameter")
    }
    if (params.config_file && !(params.config_file ==~ safePathPattern)) {
        error("Invalid characters in --config_file parameter")
    }

    // Build deli CLI arguments
    def deli_args = ""
    if (params.debug) {
        deli_args += " --debug"
    }
    if (params.deli_data_dir) {
        deli_args += " --deli-data-dir '${params.deli_data_dir}'"
    }
    if (params.config_file) {
        deli_args += " --config-file '${params.config_file}'"
    }

    extract = ExtractSequenceFiles(selection_file_path)
    def final_prefix = params.prefix ?: null

    prefix_ch = extract.selection_id
        .splitText()
        .map { it.trim() }
        .map { final_prefix ?: it }

    fastq_chunks = extract.files
        .splitText()
        .map { it.trim() }
        .map { file(it) }
        .splitFastq(by: params.chunk_size, file: true)

    decoded = DecodeChunk(fastq_chunks, selection_file_path, prefix_ch, Channel.value(deli_args))

    collected_decodes = CollectDecodeChunks(
        decoded.decoded_tsv.collect(),
        prefix_ch,
        Channel.value(deli_args)
    )

    merged_stats = MergeDecodeStatistics(
        decoded.decode_stats.collect(),
        selection_file_path,
        prefix_ch,
        Channel.value(deli_args)
    )

    WriteDecodeReport(
        merged_stats.merged_stats,
        prefix_ch,
        Channel.value(deli_args)
    )

    count_chunks = collected_decodes.ndjson
        .splitText(by: 500_000)

    counts = CountChunk(
        count_chunks,
        prefix_ch,
        Channel.value(deli_args)
    )

    collected_counts = CollectCountChunks(
        counts.counted.collect(),
        prefix_ch
    )

    SummarizeDecodeRun(
        collected_counts.merged_counts,
        merged_stats.merged_stats,
        prefix_ch,
        Channel.value(deli_args)
    )
}
