"""command line functions for deli"""

import datetime
import json
import os
from collections import Counter
from csv import DictWriter
from pathlib import Path

import click

from deli.configure import init_deli_config_dir, init_deli_data_directory, load_deli_config
from deli.decode import (
    BarcodeCaller,
    BarcodeMatcher,
    DecodeReportStats,
    DELExperiment,
    build_decoding_report,
)
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
def cli():
    """Main command group entry"""
    pass


@cli.group()
def config():
    """Group for config related commands"""
    pass


@config.command(name="init-config")
@click.argument(
    "path",
    type=click.Path(exists=True),
    required=False,
    default=os.path.join(os.getcwd(), ".deli"),
)
@click.option(
    "--include-data-dir",
    is_flag=True,
    help="include a default deli data directory in the new config directory",
)
@click.option(
    "--include-hamming", is_flag=True, help="include a default hamming files (if include-data-dir)"
)
@click.option(
    "--use-extra-parity-hamming",
    is_flag=True,
    help="include a extra parity hamming files (if include-hamming)",
)
def init_config(path, include_data_dir, include_hamming, use_extra_parity_hamming):
    """Initialize the configuration directory"""
    init_deli_config_dir(
        path,
        fail_on_exist=True,
        include_deli_data_dir=include_data_dir,
        create_default_hamming_files=include_hamming,
        use_extra_parity=use_extra_parity_hamming,
    )


@config.command(name="init-deli-data-dir")
@click.argument(
    "path",
    type=click.Path(exists=True),
    required=False,
    default=os.path.join(os.getcwd(), ".deli"),
)
@click.option(
    "--include-hamming", is_flag=True, help="include a default hamming files (if include-data-dir)"
)
@click.option(
    "--use-extra-parity-hamming",
    is_flag=True,
    help="include a extra parity hamming files (if include-hamming)",
)
def init_deli_data_dir(path, include_hamming, use_extra_parity_hamming):
    """Initialize the configuration directory"""
    init_deli_data_directory(
        path,
        fail_on_exist=True,
        create_default_hamming_files=include_hamming,
        use_extra_parity=use_extra_parity_hamming,
    )


@cli.command()
@click.argument("fastq_file", type=click.Path(exists=True), required=True)
@click.argument("experiment_file", type=click.Path(exists=True), required=True)
@click.option(
    "--out_dir", "-o", type=click.Path(), required=False, default="", help="Output directory"
)
@click.option(
    "--prefix", "-p", type=click.STRING, required=False, default="", help="Prefix for output files"
)
@click.option("--debug", is_flag=True, help="Enable debug mode")
@click.option("--save_report_data", is_flag=True, help="Save data required for reporting")
@click.option("--skip_report", is_flag=True, help="Do not render a decoding report")
@click.option(
    "--save_failed_calls", is_flag=True, help="save failed calls to a separate failed call file"
)
@click.option(
    "--deli-config",
    "-o",
    type=click.Path(),
    required=False,
    default="",
    help="Path to DELi config file to use",
)
def decode(
    fastq_file,
    experiment_file,
    out_dir,
    prefix,
    debug,
    save_report_data,
    skip_report,
    save_failed_calls,
    deli_config,
):
    """
    run decoding on a given fastq file of DEL sequences

    Parameters
    ----------
    fastq_file: Path
        the path to a fastq file
    experiment_file: Path
        the path to a DEL experiment outline
    out_dir: Path
        the path to a directory to write the output files
    prefix: str, default=""
        a prefix to append to the output files
    debug: bool, default=False
        enable debug logging mode
    save_report_data: bool, default=False
        save data required for reporting
    skip_report: bool, default=False
        do not render a decoding report
    save_failed_calls: bool, default=False
        save failed calls to a separate failed call file
    deli_config: str, default=""
        Path to DELi config file to use
    """
    logger = setup_logger("deli-decode", debug=debug)
    out_dir = _setup_outdir(out_dir)

    if deli_config != "":
        load_deli_config(deli_config)

    _start = datetime.datetime.now()

    # load in the experiment data
    logger.debug(f"reading experiment settings from {experiment_file}")
    _experiment = DELExperiment.load_experiment(experiment_file)
    _experiment_name = os.path.basename(experiment_file).split(".")[0]
    if _experiment.indexes:
        logger.info(f"loaded experiment {_experiment_name}")
        logger.debug("detected indexes; turning on de-multiplexing")
        logger.debug(f"detected indexes {[idx.index_id for idx in _experiment.indexes]}")
    logger.debug(f"detected libraries {[lib.library_id for lib in _experiment.libraries]}")

    primer_experiments = _experiment.break_into_matching_experiments()
    logger.debug(f"detected {len(primer_experiments)} primer experiments")

    # check for hamming encoding (for logging)
    _hamming_encoded_dict = {
        library.library_id: library.barcode_schema.is_hamming_encoded()
        for primer_experiment in primer_experiments
        for library in primer_experiment.libraries
    }
    for _lib_id, _is_hamming_enc in _hamming_encoded_dict.items():
        logger.debug(f"library {_lib_id} {'IS' if _is_hamming_enc else 'IS NOT'} hamming encoded")
    _has_any_hamming_encoding = any(_hamming_encoded_dict.values())  # overall hamming in use

    # write the sub experiment files
    _sub_experiment_filenames = []
    for i, primer_experiment in enumerate(primer_experiments):
        _sub_experiment_filename = out_dir / f"{_experiment_name}_sub_{i}.txt"
        _sub_experiment_filenames.append(_sub_experiment_filename)
        primer_experiment.to_experiment_file(_sub_experiment_filename)
        logger.debug(f"saving experiment {i} file at {_sub_experiment_filename}")

    logger.info(f"reading sequences from {fastq_file}")
    sequences = read_fastq(fastq_file)
    _num_reads = len(sequences)
    seq_lengths = Counter([len(seq) for seq in sequences])
    logger.info(f"read in {_num_reads:,} reads")

    # prep call file
    call_file_name = f"{prefix}_{_experiment_name}_{_timestamp()}_calls.csv"
    failed_call_file_name = f"{prefix}_{_experiment_name}_{_timestamp()}_calls_FAILED.csv"
    call_file_path = out_dir / call_file_name
    failed_call_file_path = out_dir / failed_call_file_name
    logger.debug(f"writing calls to {call_file_path}")
    if save_failed_calls:
        logger.debug(f"saving failed calls to {failed_call_file_path}")
    else:
        logger.debug("ignoring failed called")

    _headers_fieldnames = list()
    for lib in _experiment.libraries:
        _headers_fieldnames.extend(lib.barcode_schema.to_csv_header())
    _headers_fieldnames = list(set(_headers_fieldnames))
    _headers_fieldnames.extend(["from_read", "from_match", "match_seq"])

    with open(call_file_path, "w") as csv_file:
        csv_writer = DictWriter(csv_file, fieldnames=_headers_fieldnames)
        csv_writer.writeheader()

    _total_match_count = 0
    _total_match_too_big = 0
    _total_match_too_small = 0
    _total_valid_match_count = 0
    _reads_with_matches = set()
    _total_call_count = 0
    _total_valid_call_count = 0
    _reads_with_calls = set()

    _lib_calls = Counter()
    _index_calls = Counter()

    # Experiment matching loop
    for i, primer_experiment in enumerate(primer_experiments):
        _match_start = datetime.datetime.now()
        logger.info(f"running matching experiment {i + 1} of {len(primer_experiments)}")

        matcher = BarcodeMatcher(primer_experiment, **primer_experiment.matching_settings.__dict__)

        logger.debug(f"detected primer sequence {primer_experiment.primer}")
        if matcher.primer_match_length != -1:
            logger.debug("matching to truncated primer region")
        logger.debug(f"matching to pattern {matcher.pattern}")
        logger.debug(f"using error tolerance of {matcher.error_tolerance}")

        matches = matcher.match(sequences)

        _num_valid_matches = sum([_match.passed for _match in matches])
        _seqs_matched = set([_match.sequence.read_id for _match in matches if _match.passed])
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
        _total_match_too_big += sum([m.match_type() == "SequenceTooBig" for m in matches])
        _total_match_too_small += sum([m.match_type() == "SequenceTooSmall" for m in matches])

        # calling loop
        _sub_total_call_count = 0
        _sub_total_valid_call_count = 0
        _sub_reads_with_calls = set()
        _uncalled_matches = [m for m in matches if m.passed]
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
                **primer_experiment.caller_settings.__dict__,
            )
            logger.debug(f"index calling turned {'off' if caller.skip_calling_index else 'on'}")
            logger.debug(f"library calling turned {'off' if caller.skip_calling_lib else 'on'}")
            _hamming_turned_on = primer_experiment.caller_settings.__dict__["hamming"]
            logger.debug(f"hamming decoding turned {'on' if _hamming_turned_on else 'off'}")
            if not _has_any_hamming_encoding and _hamming_turned_on:
                logger.warning(
                    "hamming decoding turned on but no libraries have hamming encoding; "
                    "set 'hamming' to 'false' for 'call_settings' in experiment file"
                )
            if _has_any_hamming_encoding and not _hamming_turned_on:
                logger.warning(
                    "hamming decoding turned off but some libraries have hamming encoding; "
                    "set 'hamming' to 'true' for 'call_settings' in experiment file"
                )
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

            _lib_calls += Counter([str(c.library_call) for c in calls])
            if _experiment.has_index():
                _index_calls += Counter([str(c.index_call) for c in calls])

            for _seq in _seq_with_calls:
                _sub_reads_with_calls.add(_seq)

            _num_passed_calls = len(_seq_with_calls)
            _call_attempts = len(calls)
            _sub_total_call_count += _call_attempts
            _sub_total_valid_call_count += _num_passed_calls

            logger.debug(
                f"calling cycle {j} called {_num_passed_calls:,} "
                f"matches out of {_call_attempts:,} attempts"
            )

            # split out calls
            _passed_calls = []
            _failed_calls = []
            for c in calls:
                if c.called_successfully():
                    _passed_calls.append(c)
                else:
                    _failed_calls.append(c)

            with open(call_file_path, "a") as csv_file:
                csv_writer = DictWriter(csv_file, fieldnames=_headers_fieldnames)
                csv_writer.writeheader()
                for _call in _passed_calls:
                    csv_writer.writerow(_call.to_row_dict())

            logger.debug(
                f"wrote calls {_call_attempts:,} for calling cycle {j} to {call_file_path}"
            )

            # write failed calls if asked
            if save_failed_calls:
                with open(failed_call_file_path, "a") as csv_file:
                    csv_writer = DictWriter(csv_file, fieldnames=_headers_fieldnames)
                    csv_writer.writeheader()
                    for _call in _failed_calls:
                        csv_writer.writerow(_call.to_row_dict())

                logger.debug(
                    f"wrote {len(_failed_calls)} failed calls for "
                    f"calling cycle {j} to {failed_call_file_path}"
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

    logger.info(f"called {len(_reads_with_calls):,} out of {_num_reads:,} reads")
    logger.info(f"calls written to {call_file_path}")

    # handle report generation
    # only need to do this if report is turned on
    if (not skip_report) or save_report_data:
        # need to add in the missing libraries/index if 0 calls made
        for lib in _experiment.libraries:
            if lib.library_id not in _lib_calls.keys():
                _lib_calls[lib.library_id] = 0

        for idx in _experiment.indexes:
            if idx.index_id not in _index_calls.keys():
                _index_calls[idx.index_id] = 0

        _report_stats = DecodeReportStats(
            num_reads=_num_reads,
            read_length=seq_lengths,
            num_match_attempts=_total_match_count,
            num_call_attempts=_total_call_count,
            num_valid_matches=_total_valid_match_count,
            num_valid_calls=_total_valid_call_count,
            num_reads_with_match=len(_reads_with_matches),
            num_reads_with_calls=len(_reads_with_calls),
            num_match_too_big=_total_match_too_big,
            num_match_too_small=_total_match_too_small,
            experiment_name=_experiment_name,
            libraries=_lib_calls,
            indexes=_index_calls,
        )

        if save_report_data:
            logger.debug("saving report data files")
            json.dump(
                _report_stats.__dict__,
                open(
                    out_dir / f"{prefix}_{_experiment_name}_{_timestamp()}_report_stats.json", "w"
                ),
            )

        if not skip_report:
            logger.debug("generating html decoding report")
            report_path = (
                out_dir / f"{prefix}_{_experiment_name}_{_timestamp()}_decode_report.html"
            )
            build_decoding_report(report_stats=_report_stats, out_path=report_path)
            logger.info(f"decoding report written to {report_path}")
    else:
        logger.debug("reporting turned off")

    logger.info(
        f"DELi decoding for experiment {_experiment_name} "
        f"completed in {datetime.datetime.now() - _start}"
    )


@cli.group()
def report():
    """Report command group"""
    pass


@report.command()
@click.argument("report_stats", nargs=-1, type=click.Path(exists=True, dir_okay=False))
@click.option(
    "-n",
    "--name",
    type=click.STRING,
    required=False,
    default="merged_decoding_report",
    help="name to give the merged report file",
)
@click.option("--render_report", is_flag=True, help="render merged report as html file")
def merge(report_stats, name, render_report):
    """
    Merge a set of decode reports into a single report stat file.

    Optionally renders an html report if --render_report is used
    Will always created the merged report in the current working directory
    """
    report_stat = DecodeReportStats.load_report_file(report_stats[0])
    for file in report_stats[1:]:
        report_stat = report_stat + DecodeReportStats.load_report_file(file)

    print(render_report)

    if render_report:
        build_decoding_report(report_stat, os.path.join(os.getcwd(), f"{name}.html"))

    json.dump(report_stat.__dict__, open(os.path.join(os.getcwd(), f"{name}.json"), "w"))


@report.command()
@click.argument("report_stat", nargs=1, type=click.Path(exists=True, dir_okay=False))
def render(report_stat):
    """
    Given a report stat json file, render the decoding html report from it

    Will create the html report in the current working directory
    """
    name = os.path.basename(os.path.abspath(report_stat)).split(".")[0] + ".html"
    _report_stat = DecodeReportStats.load_report_file(report_stat)
    build_decoding_report(_report_stat, os.path.join(os.getcwd(), name))
