"""command line functions for deli"""

import configparser
import datetime
import os
from pathlib import Path

import click

from deli.configure import (
    fix_deli_data_directory,
    init_deli_config_dir,
    init_deli_data_directory,
    set_deli_data_dir,
    validate_deli_data_dir,
)
from deli.decode import (
    DecodeStatistics,
    DecodingRunner,
    build_decoding_report,
)


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


@config.command(name="init")
def init_config():
    """
    Initialize the configuration directory

    PATH is the path to the deli config directory to initialize.
    If not provided, defaults to ~/.deli/.deli
    """
    init_deli_config_dir(
        None,
        fail_on_exist=True,
        include_deli_data_dir=False,
        create_default_hamming_files=True,
        use_extra_parity=True,
    )


@cli.group()
def data():
    """Group for config related commands"""
    pass


@data.command(name="init")
@click.argument(
    "path",
    type=click.Path(exists=True),
    required=False,
    default="./deli_data",
)
def init_deli_data_dir(path):
    """
    Initialize the configuration directory

    PATH is the path to the deli data directory to initialize.
    If not provided, defaults to CWD './deli_data'.
    """
    init_deli_data_directory(path, fail_on_exist=True)


@data.command(name="fix")
@click.argument(
    "path",
    type=click.Path(exists=True),
    required=False,
    default="./",
)
@click.option(
    "--overwrite-hamming",
    is_flag=True,
    help="overwrite/fix the hamming files in the data directory",
)
def fix_deli_data_dir(path, overwrite_hamming):
    """
    Fix a deli data directory that has missing subfolders or files

    PATH is the path to the deli data directory to fix.
    """
    path = Path(path)
    if not path.exists():
        raise click.ClickException(f"Path {path} does not exist")
    if not path.is_dir():
        raise click.ClickException(f"Path {path} is not a directory")
    fix_deli_data_directory(path, overwrite_hamming=overwrite_hamming)


@data.command(name="set")
@click.argument("path", type=click.Path(exists=True), required=True)
@click.option(
    "--update-config",
    is_flag=True,
    help="Update the DELi config to use this data directory as default",
)
def set_deli_data_dir_command(path, update_config):
    """
    Set the DELi data directory to use for decoding

    PATH is the path to the deli data directory to set.
    """
    validate_deli_data_dir(path)
    os.environ["DELI_DATA_DIR"] = str(path)

    # update the config with the new data directory if requested
    if update_config:
        if os.path.exists(os.path.join(os.path.expanduser("~"), ".deli", ".deli")):
            deli_config_path = os.path.join(os.path.expanduser("~"), ".deli", ".deli")
            _config = configparser.RawConfigParser()
            _config.read(os.path.normpath(deli_config_path))
            _config["SETTINGS"]["deli_data_dir"] = str(path)
            _config.write(open(deli_config_path, "w"), True)
        else:
            raise FileNotFoundError("cannot find DELi config file for user")


@cli.group()
def decode():
    """Decode command group"""
    pass


@decode.command(name="run")
@click.argument("decode", type=click.Path(exists=True), required=True)
@click.option(
    "--out-dir",
    "-o",
    type=click.Path(),
    required=False,
    default="./",
    help="Output directory",
)
@click.argument("fastq_files", nargs=-1, type=click.Path(exists=True), required=False)
@click.option(
    "--ignore-decode-seqs",
    "-i",
    is_flag=True,
    help="Ignore the fastq sequence files in the decode file",
)
@click.option(
    "--prefix", "-p", type=click.STRING, required=False, default="", help="Prefix for output files"
)
@click.option(
    "--save-failed", is_flag=True, help="Save failed decoding results to a separate file"
)
@click.option("--tqdm", "-t", is_flag=True, help="Show tqdm progress")
@click.option("--debug", is_flag=True, help="Enable debug mode")
@click.option("--disable-logging", is_flag=True, help="Turn off DELi logging")
@click.option("--skip-report", is_flag=True, help="Skip generating the decoding report at the end")
@click.option(
    "--deli-data-dir",
    type=click.Path(),
    required=False,
    default=None,
    help="Path to DELi data directory to read libraries from",
)
def run_decode(
    decode_,
    fastq_files,
    ignore_decode_seqs,
    out_dir,
    prefix,
    save_failed,
    tqdm,
    debug,
    disable_logging,
    skip_report,
    deli_data_dir,
):
    """
    Run decoding on a given fastq file of DEL sequences

    DECODE is the path to a YAML file describing the decoding run settings.
    FASTQ_FILES is a path, list of paths or glob of FASTQ files to decode.

    NOTE: if the DECODE file contains a `selection` field, it will be used to select the
    """
    if deli_data_dir is not None:
        set_deli_data_dir(deli_data_dir)

    runner = DecodingRunner.from_file(
        decode_,
        fastq_files,
        ignore_decode_seqs=ignore_decode_seqs,
        debug=debug,
        disable_logging=disable_logging,
    )
    save_failed_to = out_dir if save_failed else None
    runner.run(save_failed_to=save_failed_to, use_tqdm=tqdm)

    runner.logger.info(f"Saving outputs to {out_dir}")
    runner.write_cube(out_dir=out_dir, prefix=prefix)
    runner.write_decode_statistics(out_dir=out_dir, prefix=prefix)
    if not skip_report:
        runner.write_decode_report(out_dir=out_dir, prefix=prefix)


@decode.group()
def statistics():
    """Decode statistics command group"""
    pass


@statistics.command(name="merge")
@click.argument("statistics", type=click.Path(exists=True), required=True, nargs=-1)
@click.option(
    "-o",
    "--out-path",
    type=click.STRING,
    required=False,
    default="merged_decoding_statistics.json",
    help="location to save merged statistics file",
)
def merge_statistics_file(statistics_, out_path):
    """
    Merge multiple decode statistics files into one

    STATISTICS is a list of paths to decode statistics files to merge.
    """
    loaded_statistics = [DecodeStatistics.from_file(p) for p in statistics_]
    merged_stats = sum(loaded_statistics, DecodeStatistics())
    merged_stats.to_file(out_path)


@decode.group()
def report():
    """Decode report command group"""
    pass


@report.command("merge")
@click.argument("decode", type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument("statistics", nargs=-1, type=click.Path(exists=True, dir_okay=False))
@click.option(
    "-o",
    "--out-path",
    type=click.STRING,
    required=False,
    default="merged_decoding_report.html",
    help="location to dave merged report",
)
def merge(decode_, statistics_, out_path):
    """
    Given decode settings and set of decode stats, merge into a single report

    Should only be used to merge stats generated from the same decode experiment
    Helpful for parallelization stuff

    DECODE is the path to a YAML file describing the decoding run.
    STATISTICS is a list of paths to decode statistics files to merge.
    """
    loaded_decode_settings = DecodingRunner.from_file(decode_)

    loaded_statistics: list[DecodeStatistics] = [
        DecodeStatistics.from_file(p) for p in statistics_
    ]
    merged_stats = sum(loaded_statistics, DecodeStatistics())

    build_decoding_report(
        selection=loaded_decode_settings.selection,
        stats=merged_stats,
        out_path=out_path,
    )


@report.command(name="generate")
@report.command("merge")
@click.argument("decode", type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument("statistic", type=click.Path(exists=True, dir_okay=False), required=True)
@click.option(
    "-o",
    "--out-path",
    type=click.STRING,
    required=False,
    default="decoding_report.html",
    help="location to dave merged report",
)
def generate(decode_, statistic_, out_path):
    """
    Generate a decoding report from a decode run config and statistic file

    DECODE is the path to a YAML file describing the decoding experiment.
    STATISTIC is the path to a decode statistics file to use for the report.
    """
    loaded_decode_settings = DecodingRunner.from_file(decode_)
    loaded_statistic = DecodeStatistics.from_file(statistic_)

    build_decoding_report(
        selection=loaded_decode_settings.selection,
        stats=loaded_statistic,
        out_path=out_path,
    )
