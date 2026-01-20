"""command line functions for deli"""

import ast
import configparser
import csv
import datetime
import logging
import os
import sys
from pathlib import Path

import click
from tqdm import tqdm

from deli import __version__
from deli.configure import DeliDataDirError, DeliDataNotFound
from deli.dels.combinatorial import CombinatorialLibrary
from deli.selection import DELSelection, Selection, SequencedSelection


# Suppress RDKit warnings at module level
try:
    from rdkit import RDLogger

    rd_logger = RDLogger.logger()
    rd_logger.setLevel(RDLogger.ERROR)
except ImportError:
    pass


def _custom_excepthook(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        # Don't log KeyboardInterrupt as an error
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logging.exception("Uncaught exception: ", exc_info=(exc_type, exc_value, exc_traceback))


def _timestamp() -> str:
    """Get the time as a timestamp string"""
    return datetime.datetime.now().strftime("%m_%d_%Y_%H%M%S%f")


def _load_any_selection(selection: os.PathLike) -> Selection:
    """
    Load in a selection file trying all known selection types

    Notes
    -----
    Prefer SequencedSelection > DELSelection > Selection
    """
    logger = logging.getLogger("deli")
    try:
        try:
            return SequencedSelection.from_yaml(selection)
        except KeyError as e_:
            if ("sequence_files" in str(e_)) or ("libraries" in str(e_)):
                logger.debug("selection file lacks sequencing info")
            else:
                raise e_
        try:
            return DELSelection.from_yaml(selection)
        except (DeliDataNotFound, FileNotFoundError):
            logger.debug(f"failed to load libraries from selection file '{selection}'")
        except KeyError as e_:
            if "libraries" not in str(e_):
                raise e_
            logger.debug("selection file lacks DEL info")
        return Selection.from_yaml(selection)
    except Exception as e_:
        logger.exception(f"failed to parse selection from '{selection}': {e_}'")
        click.echo(f"failed to parse selection file {selection}; is this a valid DELi selection YAML file?")
        sys.exit(1)


# set up root logger
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logging.captureWarnings(True)
sys.excepthook = _custom_excepthook


@click.group()
@click.version_option(
    str(__version__),
    "--version",
    "-v",
)
@click.option("--debug", is_flag=True, help="Enable debug mode")
@click.option("--disable-logging", is_flag=True, help="Turn off DELi logging")
@click.option(
    "--deli-data-dir",
    type=click.Path(),
    required=False,
    default=None,
    help="The DELi data directory to use; overrides all other settings",
)
@click.option(
    "--config-file",
    type=click.Path(exists=True, dir_okay=False),
    required=False,
    default=None,
    help="Path to DELi config file to use; if not provided, will use default at ~/.deli",
)
@click.pass_context
def cli(ctx, debug, disable_logging, deli_data_dir, config_file):
    """Main command group entry"""
    from deli.configure import get_deli_config

    ctx.ensure_object(dict)

    # prepare logger
    logger = logging.getLogger("deli")
    if debug:
        logger.setLevel(logging.DEBUG)
    if disable_logging:
        logging.disable(logging.CRITICAL + 1)  # disable all logging
    ctx.obj["logger"] = logger

    # set the deli data directory if provided
    if deli_data_dir is not None:
        from deli.configure import set_deli_data_dir

        try:
            set_deli_data_dir(deli_data_dir)
        except Exception as e:
            logger.error(f"error setting DELi data directory to '{deli_data_dir}': {e}")
            click.echo(f"failed to set DELi data directory to '{deli_data_dir}': {e}")
            sys.exit(1)

    if config_file is not None:
        from deli.configure import load_deli_config

        try:
            load_deli_config(config_file)
        except Exception as e:
            logger.error(f"error loading DELi config file from '{config_file}': {e}")
            click.echo(f"failed to load DELi config file from '{config_file}': {e}")
            sys.exit(1)

    try:
        deli_config = get_deli_config()
    except Exception as e:
        logger.error(f"error loading DELi configuration: {e}")
        click.echo(f"failed to load DELi configuration: {e}")
        sys.exit(1)

    ctx.obj["deli_config"] = deli_config
    logger.debug(f"using DELi config file at: '{deli_config.location}'")
    logger.debug(f"using DELi Data Directory at: '{deli_config.deli_data_dir}'")


@cli.group(name="config")
@click.pass_context
def config_group(ctx):
    """Group for config related commands"""
    pass


@config_group.command(name="init")
@click.argument("path", type=click.Path(), required=False, default=None)
@click.option("--overwrite", "-o", is_flag=True, help="Overwrite any existing config file")
@click.pass_context
def click_init_deli_config(ctx, path, overwrite):
    """
    Create a default DELi configuration file

    PATH is the path to the deli config directory to initialize.
    If not provided, defaults to ~/.deli
    """
    from deli.configure import init_deli_config

    _path = Path(path).absolute() if path is not None else Path.home() / ".deli"

    try:
        init_deli_config(_path, fail_on_exist=not overwrite)
    except FileExistsError:
        click.echo(
            f"'{_path}' already exists; config file not created\n"
            f"Use `deli config init --overwrite` to overwrite existing config file"
        )
        sys.exit(1)


@cli.group(name="data")
@click.pass_context
def data(ctx):
    """Group for config related commands"""
    pass


@data.command(name="init")
@click.argument(
    "path",
    type=click.Path(),
    required=False,
    default="./deli_data",
)
@click.option("--fix-missing", "-f", is_flag=True, help="Add missing sub-directories to a DELi Data Directory")
@click.option("--overwrite", "-o", is_flag=True, help="Overwrite any existing data directories")
def click_init_deli_data_dir(path, fix_missing, overwrite):
    """
    Initialize the configuration directory

    PATH is the path to the deli data directory to initialize.
    Will initialize in the CWD if not provided.

    NOTE: fix-missing will not overwrite existing subdirectories, while overwrite will.
    The 'overwrite' option will add any missing subdirectories, but also replace existing ones.
    """
    from deli.configure import init_deli_data_directory

    _path = Path(path).resolve()
    try:
        init_deli_data_directory(_path, fail_on_exist=not (fix_missing or overwrite), overwrite=overwrite)
    except FileExistsError:
        click.echo(
            f"'{_path}' already exists\n"
            f"you can create a new DELi data directory with this name using "
            f"'deli data init --overwrite {_path}'"
        )
        sys.exit(1)
    except NotADirectoryError:
        click.echo(
            f"'{_path}' is not a directory; cannot be fixed\n"
            f"you can create a new DELi data directory with this name using "
            f"'deli data init --overwrite {_path}'"
        )
        sys.exit(1)


@data.command(name="which")
@click.pass_context
def click_which_deli_data_dir(ctx):
    """
    Print the current DELi data directory

    If the DELi data directory is not set, will print a message
    and exit with a non-zero status code.
    """
    try:
        _path = ctx.obj["deli_config"].deli_data_dir
    except DeliDataDirError:
        click.echo("DELi data directory is not set")
        sys.exit(1)
    click.echo(f"Current DELi data directory: {_path}")


@data.command(name="set")
@click.argument("path", type=click.Path(), required=True)
@click.option(
    "--update-config",
    "-u",
    is_flag=True,
    help="Update the DELi config to use this data directory as default",
)
@click.pass_context
def click_set_deli_data_dir(ctx, path, update_config):
    """
    Set the DELi data directory to use for decoding

    PATH is the path to the deli data directory to set.

    NOTE: if not using --update-config, you will need to set the DELI_DATA_DIR environment variable
    manually; the command required will be printed after running.
    """
    from deli.configure import validate_deli_data_dir

    _path = Path(path).resolve()
    try:
        validate_deli_data_dir(_path)
    except FileNotFoundError:
        click.echo(
            f"Directory at '{_path}' does not exist.\n"
            f"You can create a new DELi data directory using 'deli data init {_path}'"
        )
        sys.exit(1)
    except NotADirectoryError:
        click.echo(
            f"'{_path}' is not a directory.\n"
            f"You can create a new DELi data directory with this name using "
            f"'deli data init --overwrite {_path}'"
        )
        sys.exit(1)
    except DeliDataDirError:
        click.echo(
            f"DELi data directory '{_path}' is missing required sub-directories\n"
            f"use 'deli data init --fix-missing {_path}' to add missing sub-directories"
        )
        sys.exit(1)

    # update the config with the new data directory if requested
    if update_config:
        logger = ctx.obj["logger"]
        _config_path = ctx.obj["deli_config"].location
        if _config_path.exists():
            _config = configparser.RawConfigParser()
            try:
                _config.read(os.path.normpath(_config_path))
                _config["deli.data"]["deli_data_dir"] = str(_path.resolve())
            except Exception as e:
                logger.error(e)
                click.echo(f"Failed to parse DELi config file at '{_config_path}'\nIs this a valid DELi config file?")
                sys.exit(1)
            _config.write(open(_config_path, "w"), True)
            logger.debug(f"updated DELi config with deli data directory at '{_path}'")
            click.echo(f"DELi config file at {_config_path} updated")
            sys.exit(0)
        else:
            click.echo(f"Cannot find DELi config file at {_config_path}")
            sys.exit(1)

    # set the environment variable
    import platform

    if platform.system() == "Windows":
        _command = f"set DELI_DATA_DIR={_path}"
    else:
        _command = f"export DELI_DATA_DIR={_path}"
    click.echo(
        f"To set the deli data directory for all processes in this shell, run '{_command}'\n"
        f"You can also specify a DELi data directory for any command: "
        f"'deli --deli-data-dir {_path} <REST OF COMMAND>'"
    )


@cli.group(name="decode")
@click.pass_context
def decode_group(ctx):
    """Group for decoding related commands"""
    pass


@decode_group.command(name="run")
@click.argument("decode-file", type=click.Path(exists=True), required=True)
@click.option(
    "--out-dir",
    "-o",
    type=click.Path(),
    required=False,
    default="./",
    help="Output directory to save results to",
)
@click.option("--prefix", "-p", type=click.STRING, required=False, default="", help="Prefix for output files")
@click.option("--show-tqdm", "-t", is_flag=True, help="Show tqdm progress")
@click.option("--save-fastq-info", "-q", is_flag=True, help="Save fastq file info to output")
@click.option("--save-failed", "-f", is_flag=True, help="Save failed decoding results to a separate file")
@click.option("--split-by-lib", "-s", is_flag=True, help="Save decoded results split by library")
@click.option("--exclude-score", is_flag=True, help="Exclude scores from decoding results")
@click.option("--skip-report", is_flag=True, help="Skip generating the decoding report at the end")
@click.pass_context
def run_decode(
    ctx,
    decode_file,
    out_dir,
    prefix,
    show_tqdm,
    save_fastq_info,
    save_failed,
    split_by_lib,
    exclude_score,
    skip_report,
):
    """
    Decode DNA sequences from a DEL selection

    DECODE-FILE is the path to a YAML configure file describing the DEL selection,
    sequences, and decoding settings.

    Some outputs are generated on the fly, stopping a job mid run can result
    in partial output files.

    See the docs for a detailed description of output files generated.
    """
    from deli.decode.base import FailedDecodeAttempt
    from deli.decode.decoder import DecodingSettings, SelectionDecoder
    from deli.dna.io import get_reader
    from deli.selection import SequencedSelection

    logger = ctx.obj["logger"]

    # validate output directory
    out_dir_path = Path(out_dir).absolute()
    out_dir_path.mkdir(parents=True, exist_ok=True)
    logger.debug(f"using output directory at: '{out_dir_path}'")

    # load in the sequenced selection file
    selection: SequencedSelection = SequencedSelection.from_yaml(decode_file)
    logger.info(f"loaded selection '{selection.selection_id}' from '{decode_file}'")
    logger.debug(
        f"detected {len(selection.library_collection)} libraries: "
        f"{[lib.library_id for lib in selection.library_collection.libraries]}"
    )
    if len(selection.tool_compounds) > 0:
        logger.debug(
            f"loaded {len(selection.tool_compounds)} tool compounds for decoding: "
            f"{[tool.compound_id for tool_lib in selection.tool_compounds for tool in tool_lib.compounds]}"
        )

    # make sequence reader
    reader = get_reader(selection.sequence_files)
    logger.debug(f"detected {len(reader.sequence_files)} sequencing files: {reader.sequence_files}")

    # deal with prefix
    if prefix is None or prefix == "":
        prefix = selection.selection_id
    logger.debug(f"using output file prefix: '{prefix}'")

    # determine output file fields
    header: list[str] = []
    if save_fastq_info:
        header.extend(["fastq_file", "read_name"])
    header.append("library_id")
    if not exclude_score:
        header.append("library_score")
    header.append("bb_ids")
    if not exclude_score:
        header.append("bb_scores")
    header.append("umi")
    if not exclude_score:
        header.append("overall_score")

    # handle output file locations
    if split_by_lib:
        decode_out_dir = out_dir_path / "decodes_by_library"
        decode_out_dir.mkdir(parents=False, exist_ok=False)
        decoded_out_writers = dict()
        for library in selection.library_collection.libraries:
            writer = csv.DictWriter(
                open(decode_out_dir / f"{prefix}_{library.library_id}.tsv", "w", newline=""),
                fieldnames=header,
                delimiter="\t",
                extrasaction="ignore",
            )
            writer.writeheader()
            decoded_out_writers[library.library_id] = writer
        logger.info(f"saving decoded sequences per library to '{decode_out_dir}'")
    else:
        decoded_out_path = out_dir_path / f"{prefix}_decoded.tsv"
        writer = csv.DictWriter(
            open(decoded_out_path, "w", newline=""), fieldnames=header, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        decoded_out_writers = {"all": writer}
        logger.info(f"writing decoded sequences to: '{decoded_out_path}'")

    # handle failed decoding output file
    failed_out_file = None
    if save_failed:
        failed_out_path = out_dir_path / f"{prefix}_failed_decoding.tsv"
        failed_out_file = open(failed_out_path, "w")
        failed_out_file.write("fastq_file\tread_name\tissue\treason\n")  # add header
        logger.info(f"writing failed decoding sequences to: '{failed_out_path}'")

    # load decode settings
    decode_settings = DecodingSettings.from_file(decode_file)
    logger.debug(f"loaded decoding settings: {decode_settings}")

    # load in decoder
    decoder = SelectionDecoder(selection=selection, decode_settings=decode_settings)

    # loop through sequences and decode
    prev_file = ""
    curr_seq_count = 0
    for i, (sequence_file, sequence_record) in tqdm(
        enumerate(reader.iter_seqs_with_filenames()),
        disable=not show_tqdm,
        desc=f"decoding reads for selection {selection.selection_id}",
    ):
        if prev_file != sequence_file.name:
            if curr_seq_count != 0:
                logger.debug(f"decoded {i - curr_seq_count} sequences from file: '{prev_file}'")
            logger.debug(f"decoding sequences from file: '{sequence_file.name}'")
            prev_file = sequence_file.name
            curr_seq_count = i

        decoded_read = decoder.decode_read(sequence_record)

        if isinstance(decoded_read, FailedDecodeAttempt):
            if failed_out_file is not None:
                line = (
                    f"{sequence_file}\t{sequence_record.name}\t{decoded_read.__class__.__name__}\t{decoded_read.reason}"
                )
                failed_out_file.write(f"{line}\n")

        else:  # write the decoded compound
            if split_by_lib:
                writer = decoded_out_writers[decoded_read.get_library_id()]
            else:
                writer = decoded_out_writers["all"]

            output_dict = decoded_read.to_decode_res_row_dict()
            if save_fastq_info:
                output_dict["fastq_file"] = str(sequence_file)
                output_dict["read_name"] = sequence_record.name
            writer.writerow(output_dict)

        if (i + 1) % 500000 == 0:
            logger.debug(f"decoded {i + 1:,} sequences...")

    logger.info(f"finished decoding for selection '{selection.selection_id}'")

    for file in decoded_out_writers.values():
        file.close()
        logger.debug(f"closed file '{file}'")

    if failed_out_file is not None:
        failed_out_file.close()
        logger.debug(f"closed file '{failed_out_file}'")

    # handle writing statistics and report file
    statistics_file = out_dir_path / "decode_statistics.json"
    decoder.decode_stats.to_file(statistics_file)
    logger.info(f"wrote decoding statistics to: '{statistics_file}'")

    # generate report if not skipped
    if not skip_report:
        from deli.decode.report import build_decoding_report

        report_file = out_dir_path / f"{prefix}_decode_report.html"
        build_decoding_report(stats=decoder.decode_stats, selection=selection, out_path=report_file)
        logger.info(f"wrote decoding report to: '{report_file}'")


@decode_group.command(name="aggregate")
@click.argument("decoded-reads", type=click.Path(exists=True), nargs=-1, required=True)
@click.option(
    "--score-threshold",
    "-s",
    type=click.INT,
    required=False,
    default=100,
    help="Reject decoded compounds above this score",
)
@click.option(
    "--count-threshold",
    "-c",
    type=click.INT,
    required=False,
    default=1,
    help="Minimum count threshold for including compounds in output",
)
@click.option(
    "--out-loc",
    "-o",
    type=click.Path(),
    required=False,
    default="./aggregated_decodes.json",
    help="Output location to save aggregated decodes to; will add .json suffix if missing",
)
@click.option("--compress", "-z", is_flag=True, help="Compress output JSON with gzip")
@click.pass_context
def aggregate_decodes(ctx, decoded_reads, score_threshold, count_threshold, out_loc, compress):
    """
    Aggregate decoded reads from decoded TSV file(s) into a JSON count format

    Aggregation is used to determine how many times each compound was decoded
    and the unique UMIs (and their counts) associated with each compound.
    """
    import polars as pl

    logger = ctx.obj["logger"]

    # validate output directory
    out_path = Path(out_loc).absolute()
    if out_path.exists() and out_path.is_dir():
        print(f"output location '{out_path}' is a directory; please provide a file path")
        sys.exit(1)
    if compress:
        if out_path.suffixes[-2:] != [".json", ".gz"]:
            if out_path.suffixes[-1] != ".json":
                out_path = out_path.with_suffix(".json.gz")
            else:
                out_path = out_path.with_suffix(".gz")
    else:
        if out_path.suffix != ".json":
            out_path = out_path.with_suffix(".json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    logger.debug(f"writing aggregated decoded sequences to: '{out_path}'")

    aggregated_decodes = (
        pl.scan_csv(decoded_reads, has_header=True, separator="\t", ignore_errors=True)
        .select(["library_id", "bb_ids", "umi", "overall_score"])
        .filter(pl.col("overall_score") <= score_threshold)
        .groupby(["library_id", "bb_ids"])
        .agg(pl.col("umi"))
        .filter(pl.col("umi").list.len() >= count_threshold)
        .with_columns(
            umi_counts=pl.col("umi").list.eval(pl.element().value_counts().struct.rename_fields(["k", "c"])),
        )
        .collect()
    )
    logger.info(f"found {len(aggregated_decodes)} compounds in decoded reads after filtering")

    if compress:
        import gzip

        logger.debug("compressing output JSON with gzip")
        with gzip.open(out_path, "wt") as out_file:
            out_file.write(aggregated_decodes.write_ndjson())
    else:
        aggregated_decodes.write_ndjson(str(out_path))
    logger.info(f"wrote aggregated decoded sequences to: '{out_path}'")


@decode_group.command(name="report")
@click.argument("decode_stats_file", nargs=-1, type=click.Path(exists=True), required=True)
@click.option(
    "--out-loc",
    "-o",
    type=click.Path(),
    required=False,
    default="./decode_report.html",
    help="Output location to save report to",
)
@click.option(
    "--selection",
    "-s",
    type=click.Path(exists=True),
    required=False,
    help="Selection file contain info about conditions for decoding",
    default=None,
)
@click.pass_context
def generate_report(ctx, decode_stats_file, out_loc, selection):
    """
    Generate an HTML decoding report from decoding statistics file(s)

    If multiple statistic files are provided, will merge them into a single report.
    This is useful for dealing with parallel decoding runs for the same selection.

    Providing the selection file used for the decoding will enable more detailed reporting,
    such as listing libraries that were present in the selection but had no decoded sequences.

    DECODE_STATS_FILE are the paths to the decoding statistics JSON files.
    """
    from deli.decode.decoder import DecodeStatistics
    from deli.decode.report import build_decoding_report

    logger = ctx.obj["logger"]

    out_loc_path = Path(out_loc).absolute()
    if out_loc_path.suffix == "":
        out_loc_path = out_loc_path / "decode_report.html"
    elif out_loc_path.suffix != ".html":
        out_loc_path = out_loc_path.with_suffix(".html")

    logger.debug(f"writing decoding report to: '{out_loc_path}'")
    out_loc_path.parent.mkdir(parents=True, exist_ok=True)

    selection_obj: Selection | None = None
    if selection is not None:
        selection_obj = _load_any_selection(selection)
        logger.debug(f"loaded selection '{selection_obj.selection_id}' from '{selection}'")

    overall_stats = DecodeStatistics()

    for stats_file in decode_stats_file:
        stats_path = Path(stats_file).absolute()
        try:
            stats = DecodeStatistics.from_file(stats_path)
            logger.debug(f"loaded decoding stats file: '{stats_path}'")
        except Exception as e:
            logger.exception(f"failed to load decoding stats file '{stats_path}': {e}'")
            click.echo(
                f"failed to parse decoding statistics file {stats_file}; is this a valid DELi decoding statistics file?"
            )
            sys.exit(1)

        overall_stats += stats

    build_decoding_report(
        stats=overall_stats,
        out_path=out_loc_path,
        selection=selection_obj,
    )
    logger.info(f"wrote decoding report to: '{out_loc_path}'")


@cli.group(name="degen")
@click.pass_context
def degen_group(ctx):
    """Group for degenerating decoded sequences commands"""
    pass


@degen_group.command(name="run")
@click.argument("aggregated-compounds", type=click.Path(exists=True), required=True)
@click.option(
    "--out-loc",
    "-o",
    type=click.Path(),
    required=False,
    default="./degenerated_compounds.tsv",
    help="Location to save results to",
)
@click.option("--raw-counts", "-r", is_flag=True, help="Include raw counts in output")
@click.option("--cpd-ids", "-d", is_flag=True, help="Include full DEL compound IDs in output")
@click.option(
    "--count-threshold",
    "-q",
    type=click.INT,
    required=False,
    default=0,
    help="Minimum score threshold for including decoded sequences in degeneration",
)
@click.option("--use-tqdm", "-t", is_flag=True, help="Show tqdm progress bar")
@click.pass_context
def run_degen(ctx, aggregated_compounds, out_loc, raw_counts, cpd_ids, count_threshold, use_tqdm):
    """
    Degenerate decoded sequences into compound counts based on UMIs

    AGGREGATED_COMPOUNDS is the path to the aggregated decoded JSON file to degenerate.

    Output will be a TSV file mapping compounds to their degenerated counts in a
    See the Docs for more information on degeneration and the meaning of certain count values.

    Note: Degeneration assumes that all decoded sequences have already been aggregated into unique
    compound + UMI pairs with counts. Use the `deli decode aggregate` command to generate
    the required aggregated decoded file from decoded TSV files.
    """
    from deli.dels.compound import generate_del_compound_id

    logger = ctx.obj["logger"]

    out_loc_path = Path(out_loc).absolute()
    if out_loc_path.exists() and out_loc_path.is_dir():
        print(f"output location '{out_loc_path}' is a directory; please provide a file path")
        sys.exit(1)
    if out_loc_path.suffix != ".tsv":
        out_loc_path = out_loc_path.with_suffix(".tsv")
    out_loc_path.parent.mkdir(parents=True, exist_ok=True)
    logger.debug(f"writing degenerate sequences to: '{out_loc_path}'")

    compounds_json_path = Path(aggregated_compounds)
    logger.info(f"loading aggregated decoded compounds from '{compounds_json_path}'")
    if Path.suffix == ".gz":
        import gzip

        logger.debug("decompressing file with gzip")
        file = gzip.open(compounds_json_path, "rt")
    else:
        file = open(compounds_json_path, "r")

    header = ["lib_id", "bb_ids", "count"]
    if raw_counts:
        header.append("raw_count")
    if cpd_ids:
        header = ["compound_id"] + header

    out_file = csv.DictWriter(
        open(out_loc_path, "w", newline=""), fieldnames=header, delimiter="\t", extrasaction="ignore"
    )

    skipped = 0
    passed = 0
    for i, line in tqdm(enumerate(file), desc="degenerating compounds", disable=not use_tqdm):
        cpd_info = ast.literal_eval(line)
        # calculate degen counts
        if cpd_ids:
            cpd_info["compound_id"] = generate_del_compound_id(cpd_info["library_id"], cpd_info["bb_ids"].split(","))
        if raw_counts:
            cpd_info["raw_count"] = sum([count_struct["c"] for count_struct in cpd_info["umi_counts"]])
        cpd_info["count"] = len(cpd_info["umi_counts"])
        # check for count threshold
        if cpd_info["count"] > count_threshold:
            out_file.writerow(cpd_info)
            passed += 1
        else:
            skipped += 1

        if ((i + 1) % 10000) == 0:
            logger.debug(f"degenerated {i + 1:,} compounds...")

    logger.info(f"degenerated {passed} compounds")
    if skipped > 0:
        logger.info(
            f"skipped {skipped} out of {skipped + passed} ({round((skipped / (skipped + passed)) * 100, 2)}%) "
            f"compounds below count threshold of {count_threshold}"
        )
    logger.info(f"wrote degenerated compounds to '{out_loc_path}'")


@cli.command(name="enumerate")
@click.argument("library_file", type=click.Path(exists=True), required=True)
@click.option("--out_path", "-o", type=click.Path(), required=False, default="", help="Output CSV file path")
@click.option("--tqdm", "-t", is_flag=True, help="Enable TQDM progress bar")
@click.option("--fail-on-error", "-f", is_flag=True, help="Fail on first error during enumeration")
@click.option("--drop-failed", "-d", is_flag=True, help="Drop compounds with failed enumerations")
def enumerate_(library_file, out_path, tqdm, fail_on_error, drop_failed):
    """
    Enumerates compounds from a given library

    If out_path is not provided, will save to the current working directory
    as a CSV file named <library_id>_enumerated.csv

    LIBRARY_FILE is the path to a DELi library file to enumerate.
    """
    library_id = os.path.basename(library_file).split(".")[0]
    output_file = out_path if out_path != "" else os.path.join(os.getcwd(), f"{library_id}_enumerated.csv")

    _start = datetime.datetime.now()

    enumerator = CombinatorialLibrary.load(library_file)
    enumerator.enumerate_to_file(
        output_file,
        separator=",",
        use_tqdm=tqdm,
        fail_on_error=fail_on_error,
        drop_failed=drop_failed,
    )
