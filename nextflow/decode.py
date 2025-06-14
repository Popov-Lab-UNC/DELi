"""nextflow process script for decoding in parallel"""

import argparse
import os

from deli.decode import DecodingRunner


parser = argparse.ArgumentParser()
parser.add_argument(
    "--fastq_file",
    type=str,
    required=True,
)
parser.add_argument(
    "--decode_file",
    type=str,
    required=True,
)
parser.add_argument(
    "--debug",
    action="store_true",
)
parser.add_argument(
    "--save_failed",
    action="store_true",
)
args = parser.parse_args()

sub_job_id = os.path.basename(args.fastq_file).split(".")[0]
out_dir = os.path.abspath("./")

runner = DecodingRunner.from_file(
    args.decode_file, fastq_files=[args.fastq_file], ignore_decode_seqs=True, debug=args.debug
)

save_failed_to = out_dir if args.save_failed else None

results = runner.run(save_failed_to=save_failed_to, use_tqdm=False)

runner.logger.info(f"Saving outputs to {out_dir}")

statistics_out_path = os.path.join(out_dir, f"{sub_job_id}_decode_statistics_subjob.json")
runner.logger.debug(f"Saving decode statistics to {statistics_out_path}")
results.write_decode_statistics(statistics_out_path, include_read_lengths=True)

counter_out_path = os.path.join(out_dir, f"{sub_job_id}_counter_subjob.json.gz")
runner.logger.debug(f"Saving counters to {counter_out_path}")
runner.degen.to_json(counter_out_path, compress=True)

os.rename("deli.log", f"subjob_{sub_job_id}_deli_decode.log")
