"""command line functions for deli"""

import datetime
import os
from csv import DictWriter
from pathlib import Path

import click

from deli.decode import BarcodeCaller, BarcodeMatch, BarcodeMatcher, DELExperiment
from deli.logging.logger import setup_logger
from deli.sequence import read_fastq


def _timestamp() -> str:
    """Get the time as a timestamp string"""
    return datetime.datetime.now().strftime("%m_%d_%Y_%H%M%S%f")


def _setup_outdir(out_dir):
    """Makes sure the output directory exists, if not makes it"""
    if out_dir == "":
        out_dir = os.getcwd()
    out_dir = Path(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    return out_dir


@click.group()
def deli():
    """Main command group entry"""
    pass


@deli.command()
@click.argument("fastq_file", type=click.Path(exists=True), required=True)
@click.option("--experiment", "-e", type=click.Path(exists=True), required=True)
@click.option("--out_dir", "-o", type=click.Path(), required=False, default="")
@click.option("--prefix", "-p", type=click.STRING, required=False, default="")
@click.option("--debug", is_flag=True, help="Enable debug mode")
def decode(fastq_file, experiment, out_dir, prefix, debug):
    """
    run decoding on a given fastq file of DEL sequences

    Parameters
    ----------
    fastq_file: Path
        the path to a fastq file
    experiment: Path
        the path to a DEL experiment outline
    out_dir: Path
        the path to a directory to write the output files
    prefix: str, default=""
        a prefix to append to the output files
    debug: bool, default=False
        enable debug logging mode
    """
    logger = setup_logger("deli-decode", debug=debug)
    out_dir = _setup_outdir(out_dir)
    _start = datetime.datetime.now()

    # load in the experiment data
    logger.debug(f"reading experiment settings from {experiment}")
    _experiment = DELExperiment.load_experiment(experiment)
    _experiment_name = os.path.basename(experiment).split(".")[0]
    if _experiment.indexes:
        logger.info(f"loaded experiment {_experiment_name}")
        logger.debug("detected indexes; turning on de-multiplexing")
        logger.debug(f"detected indexes {[idx.index_id for idx in _experiment.indexes]}")
    logger.debug(f"detected libraries {[lib.library_id for lib in _experiment.libraries]}")

    primer_experiments = _experiment.break_into_matching_experiments()
    logger.debug(f"detected {len(primer_experiments)} primer experiments")

    # write the sub experiment files
    _sub_experiment_filenames = []
    for i, primer_experiment in enumerate(primer_experiments):
        _sub_experiment_filename = out_dir / f"{_experiment_name}_sub_{i}.txt"
        _sub_experiment_filenames.append(_sub_experiment_filename)
        primer_experiment.to_experiment_file(_sub_experiment_filename)
        logger.debug(f"saving experiment {i} file at {_sub_experiment_filename}")

    logger.info(f"reading sequences from {fastq_file}")
    sequences = read_fastq(fastq_file)
    logger.info(f"read in {len(sequences):,} reads")

    # prep call file
    call_file_name = f"{prefix}_{_experiment_name}_{_timestamp()}_calls.csv"
    call_file_path = out_dir / call_file_name
    logger.debug(f"writing calls to {call_file_path}")

    _headers_fieldnames = list()
    for lib in experiment.libraries:
        _headers_fieldnames.extend(lib.barcode_schema.to_csv_header())
    _headers_fieldnames = list(set(_headers_fieldnames))
    _headers_fieldnames.extend(["from_read", "from_match", "match_seq"])

    with open(call_file_path, "w") as csv_file:
        csv_writer = DictWriter(csv_file, fieldnames=_headers_fieldnames)
        csv_writer.writeheader()

    _total_match_count = 0
    _total_valid_match_count = 0
    _reads_with_matches = set()
    _total_call_count = 0
    _total_valid_call_count = 0
    _reads_with_calls = set()

    # Experiment matching loop
    for i, primer_experiment in enumerate(primer_experiments):
        _match_start = datetime.datetime.now()
        logger.info(f"running matching experiment {i + 1} of {len(primer_experiments)}")

        matcher = BarcodeMatcher(primer_experiment, **primer_experiment.matching_settings)
        logger.debug(f"detected primer sequence {primer_experiment.primer}")
        logger.debug(f"matching to pattern {matcher.pattern}")

        matches = matcher.match(sequences)

        _num_valid_matches = sum([isinstance(_match, BarcodeMatch) for _match in matches])
        _seqs_matched = set([_match.sequence.read_id for _match in matches])
        _num_seqs_matched = len(_seqs_matched)
        _num_matches = len(matches)
        logger.debug(
            f"matching experiment {i} made {_num_matches:,} match "
            f"attempts in {datetime.datetime.now() - _match_start}"
        )
        logger.info(f"matching experiment {i} found {_num_valid_matches:,} valid matches")

        _total_match_count += _num_matches
        _total_valid_match_count += _num_valid_matches
        _reads_with_matches = _reads_with_matches.union(_seqs_matched)

        # calling loop
        _sub_total_call_count = 0
        _sub_total_valid_call_count = 0
        _sub_reads_with_calls = set()
        _uncalled_matches = matches
        _call_start = datetime.datetime.now()

        logger.debug(
            f"detected {len(primer_experiment.library_schema_groups)} library schema groups"
        )
        for j, _sub_library_set in enumerate(primer_experiment.library_schema_groups):
            _call_cycle_start = datetime.datetime.now()
            logger.info(
                f"running calling cycle {j + 1} of "
                f"{len(primer_experiment.library_schema_groups)}"
            )
            logger.debug(
                f"detected library call region "
                f"{_sub_library_set.library_call_schema.full_barcode}"
            )

            caller = BarcodeCaller(
                libraries=_sub_library_set,
                indexes=primer_experiment.indexes,
                **primer_experiment.caller_settings,
            )
            logger.debug(f"index calling turned {'off' if caller.skip_calling_index else 'on'}")
            logger.debug(f"library calling turned {'off' if caller.skip_calling_lib else 'on'}")
            if _sub_library_set.requires_multistep_calling:
                logger.debug("turned on multistep (library first) calling")
                logger.debug(
                    f"calling library with barcode schema'"
                    f"{_sub_library_set.library_call_schema.full_barcode}'"
                )

            calls = caller.call_tags(_uncalled_matches)

            _uncalled_matches = [
                _uncalled_matches[idx]
                for idx, _call in enumerate(calls)
                if not _call.called_successfully()
            ]

            _seq_with_calls = [
                _call.parent_match.sequence.read_id
                for _call in calls
                if _call.called_successfully()
            ]
            for _seq in _seq_with_calls:
                _sub_reads_with_calls.add(_seq)
            _passed_calls = len(_seq_with_calls)
            _call_attempts = len(calls)
            _sub_total_call_count += _call_attempts
            _sub_total_valid_call_count += _passed_calls

            logger.debug(
                f"calling cycle {j} called {_passed_calls:,} "
                f"matches out of {_call_attempts:,} attempts"
            )

            with open(call_file_path, "a") as csv_file:
                csv_writer = DictWriter(csv_file, fieldnames=_headers_fieldnames)
                for _call in calls:
                    csv_writer.writerow(_call.to_row_dict())
            logger.debug(
                f"wrote calls {_call_attempts:,} for calling cycle {j} to {call_file_path}"
            )
            logger.debug(
                f"call cycle {i + 1} completed in {datetime.datetime.now() - _call_cycle_start}"
            )

        logger.debug(
            f"calling for experiment {i} attempted {_sub_total_call_count:,} calls in"
            f"{datetime.datetime.now() - _call_start}"
        )
        logger.info(f"calling for experiment {i} made {_sub_total_valid_call_count:,} valid calls")

        _total_valid_call_count += _sub_total_valid_call_count
        _total_call_count += _sub_total_call_count
        _reads_with_calls = _reads_with_calls.union(_sub_reads_with_calls)

    logger.debug(
        f"made {_total_valid_match_count:,} matches out of {_total_match_count:,} attempts"
    )
    logger.debug(f"made {_total_valid_call_count:,} calls out of {_total_call_count:,} attempts")

    logger.info(f"called {len(_reads_with_calls):,} out of {len(sequences):,} reads")
    logger.info(f"calls written to {call_file_path}")
    logger.info(
        f"DELi decoding for experiment {_experiment_name} "
        f"completed in {datetime.datetime.now() - _start}"
    )
