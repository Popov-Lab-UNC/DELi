# DELi Decoding Nextflow Workflow

A scalable parallel workflow for decoding Display-Encoded Library (DEL) selection experiments. This workflow processes FASTQ sequencing data, decodes molecular barcodes, counts compounds, and generates comprehensive analysis reports.

## Workflow Overview

The workflow implements a multi-stage pipeline optimized for parallelization:

1. **Extract Configuration** - Parses selection YAML files to extract sequence files and metadata
2. **Chunk FASTQ Files** - Splits large FASTQ files into manageable chunks (default: 1 million reads per chunk)
3. **Parallel Decoding** - Decodes molecular barcodes for each FASTQ chunk independently
4. **Collect Decode Results** - Aggregates decoded sequences from all chunks into a single NDJSON file
5. **Merge Decode Statistics** - Consolidates decoding statistics from all chunks
6. **Chunk and Count** - Further chunks the collected NDJSON data and counts compounds in parallel
7. **Merge Counts** - Combines parquet count files from all chunks into a single count matrix
8. **Summarization** - Merges decode statistics with final counts to produce comprehensive analysis
9. **Report Generation** - Creates an HTML report with decoding statistics and analysis results

## Processes

### ExtractSequenceFiles
Extracts sequence files, library information, and selection IDs from the selection configuration YAML file.

**Outputs:**
- Selection ID
- List of FASTQ file paths

### DecodeChunk
Runs `deli decode run` on individual FASTQ chunks to extract and decode molecular barcodes. Executes in parallel across all chunks.

**Resources:** 4 GB memory, 1 CPU
**Outputs:**
- Decoded TSV file
- Decode statistics JSON
- Deli log file

### CollectDecodeChunks
Collects all decoded TSV chunks into a single NDJSON file using `deli decode collect`.

**Resources:** 16 GB memory, 1 CPU (high memory due to collection operation)
**Outputs:**
- Collected NDJSON file

### MergeDecodeStatistics
Merges decode statistics JSON files from all chunks into a single aggregated statistics file.

**Resources:** 4 GB memory, 1 CPU
**Outputs:**
- Merged decode statistics JSON

### CountChunk
Counts compounds from NDJSON chunks (500,000 lines per chunk) and outputs parquet format.

**Resources:** 4 GB memory, 1 CPU
**Outputs:**
- Counted parquet file per chunk

### CollectCountChunks
Merges all parquet count files from chunked counting into a single final count matrix.

**Resources:** 16 GB memory, 1 CPU (high memory due to merging large parquet files)
**Outputs:**
- Final counts parquet file

### SummerizeDecodeRun
Combines final counts with decode statistics to produce comprehensive summary statistics.

**Resources:** 4 GB memory, 1 CPU
**Outputs:**
- Final statistics JSON

### WriteDecodeReport
Generates an HTML report from final statistics.

**Resources:** 4 GB memory, 1 CPU
**Outputs:**
- HTML report file

## Configuration

The workflow supports three execution environments through profiles defined in `nextflow.config`:

### Local Execution (Development)

```bash
nextflow run nextflow.nf -profile local
```

Runs on your local machine using Nextflow's local executor. Ideal for testing and development.

**Configuration:**
- Executor: local
- All processes execute sequentially or with limited parallelization based on available CPU cores

### Gridengine HPC Cluster

```bash
nextflow run nextflow.nf -profile gridengine
```

Submits jobs to a Gridengine (SGE) HPC cluster for large-scale processing.

**Configuration:**
- Executor: sge (Sun/Oracle Grid Engine)
- Parallel environment: smp (shared memory parallel)
- Queue: default (adjust if your cluster uses different queue names)

To customize queue assignment, modify the `nextflow.config` file:
```
queue = 'your_queue_name'  // Change from 'default' to your queue
```

### AWS Batch

```bash
nextflow run nextflow.nf -profile aws
```

Executes the workflow on AWS using AWS Batch for fully managed job scheduling and scaling.

**Configuration:**
- Executor: awsbatch
- Region: us-east-1 (adjust to your AWS region)
- Queue: default (set to your AWS Batch compute environment queue)

To customize for your AWS environment, modify `nextflow.config`:
```
region = 'your-aws-region'  // Change from 'us-east-1'
queue = 'your-queue-name'   // Set to your compute environment
```

## Resource Allocation

All processes follow a standardized resource model:

- **Standard processes:** 4 GB memory, 1 CPU core
- **Collection processes (CollectDecodeChunks, CollectCountChunks):** 16 GB memory, 1 CPU core

Collection processes require additional memory due to the overhead of merging large datasets.

## Parameters

Key parameters that can be overridden at runtime with `--parameter value`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `selection_file` | null (required) | Path to selection YAML configuration file |
| `out_dir` | `$launchDir/decode_results` | Output directory for results |
| `prefix` | "" (auto-generate) | Prefix for output files (defaults to selection_id) |
| `chunk_size` | 1,000,000 | Number of reads per FASTQ chunk |
| `deli_data_dir` | null | Path to DELi data directory (if not in default location) |
| `config_file` | null | Path to DELi configuration file (if not using default) |
| `debug` | false | Enable debug mode for deli commands |

## Usage Examples

### Basic execution on local machine
```bash
nextflow run nextflow.nf --selection_file /path/to/selection.yaml -profile local
```

### Large-scale HPC execution with custom parameters
```bash
nextflow run nextflow.nf \
  --selection_file /path/to/selection.yaml \
  --out_dir /scratch/decode_results \
  --chunk_size 2000000 \
  -profile gridengine
```

### AWS Batch execution with custom region
```bash
nextflow run nextflow.nf \
  --selection_file s3://my-bucket/selection.yaml \
  --out_dir s3://my-bucket/results \
  -profile aws
```

### Enable debug logging
```bash
nextflow run nextflow.nf \
  --selection_file /path/to/selection.yaml \
  --debug \
  -profile local
```

## Output Files

The workflow generates the following output files in the specified output directory:

- `{prefix}_decoded.tsv` - Decoded barcode sequences (per chunk, then collected)
- `{prefix}_collected.ndjson` - Aggregated decoded sequences in NDJSON format
- `{prefix}_decode_stats.json` - Aggregated decoding statistics (per-channel, error rates, etc.)
- `{prefix}_counts.parquet` - Final compound count matrix
- `{prefix}_decode_summary.json` - Comprehensive summary statistics combining counts and decode stats
- `{prefix}_decode_report.html` - HTML report with visualizations and statistics

## Requirements

- Nextflow >= 21.04.0
- Python 3.7+ with pyyaml, polars, and pandasfor data processing
- DELi package installed and available in PATH
- For HPC/AWS: appropriate cluster/cloud credentials and access

## Configuration Customization

To customize the workflow for your specific environment, edit `nextflow.config`:

1. **Adjust queue names** for your HPC cluster or AWS Batch compute environment
2. **Modify memory values** if your processes require different resources
3. **Enable Nextflow Tower** for cloud-based monitoring and reporting by setting `tower.enabled = true`
4. **Add additional profiles** for other execution environments

## Troubleshooting

### Job submission failures on HPC
- Verify queue names match your cluster configuration
- Check that the parallel environment (penv) matches your SGE setup
- Confirm your user has permissions to submit jobs to the queue

### AWS Batch errors
- Verify AWS credentials are properly configured
- Check that the compute environment queue exists and is accessible
- Ensure Docker image for job containers is properly configured

### Memory or timeout issues
- Increase `chunk_size` to reduce total chunks and collection overhead
- Increase memory allocation in `nextflow.config` for specific processes
- Check system logs for out-of-memory errors

## Performance Considerations

- Larger `chunk_size` values reduce parallelism but decrease collection overhead
- Default 1 million read chunks balance parallelization with efficient resource use
- Collection processes bottleneck on I/O; ensure adequate disk space and speed
- Consider network bandwidth for AWS S3 access when using cloud profiles

## Support

For issues related to the workflow logic, contact the DELi development team.
For Nextflow-specific questions, see [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html).
