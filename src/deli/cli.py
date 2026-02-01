"""command line functions for deli"""

import ast
import configparser
import csv
import datetime
import gzip
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

# custom exception hook to log uncaught exceptions
def _custom_excepthook(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        # Don't log KeyboardInterrupt as an error
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logging.exception("Uncaught exception: ", exc_info=(exc_type, exc_value, exc_traceback))


def _timestamp() -> str:
    """Get the time as a timestamp string"""
    return datetime.datetime.now().strftime("%m_%d_%Y_%H%M%S%f")


def _standardize_and_validate_output_loc(
        path: Path | str,
        is_file: bool,
        can_exist: bool = True,
        overwrite: bool = False,
        extension: str = None,
        is_compressed: bool = False,
) -> Path:
    """
    Standardize and validate the output location based on the given parameters.

    Will Standardize to absolute path.
    Will add ".gz" suffix if `is_compressed` is True and not already present.
    Will *not* add any other suffixes automatically, will just warn.

    Will exit and print message to CLI if validation fails.

    Parameters
    ----------
    path : Path | str
        The output path to validate.
    is_file : bool
        Whether the output path should be a file (True) or a directory (False).
    can_exist : bool, default = True
        Whether the output path is allowed to exist (even if `override` is `False`).
    overwrite : bool, default = False
        Whether to allow overwriting existing files/directories.
    extension : str, optional
        The required output file extension (e.g., '.txt').
        Only used if `is_file` is `True`.
    is_compressed : bool, default = False
        Whether the output file is compressed with gzip.
        Only used if `is_file` is `True`.

    Returns
    -------
    Path
        The validated (and possibly modified) output path.
    """

    logger = logging.getLogger("deli")

    def _communicate_error_and_exit(message: str):
        """helper func to communicate error and exit"""
        click.echo()
        if logger:
            logger.error(message)
        sys.exit(1)

    def _communicate_warning(message: str):
        """helper func to communicate warning"""
        click.echo(f"\nWARNING: {message}\n")
        if logger:
            logger.warning(message)

    path = Path(path).absolute()  # standardize to absolute path

    # add the compressed suffix if needed
    if is_compressed and is_file:
        if path.suffix != ".gz":
            msg = f"Output file '{path}' is going to be gzip compressed but is missing a .gz suffix"
            path = path.with_suffix(path.suffix + ".gz")
            msg += f"\n'.gz' added automatically; new file path: '{path}'"
            _communicate_warning(msg)

    if path.exists():
        # mismatch between expected type and actual type
        # these messages come first so a user known that --overwrite won't help here
        if not is_file and path.is_file():
            msg = (f"Path '{path}' is an existing, not a directory\n"
                   f"--overwrite will not overwrite a file with a directory")
            _communicate_error_and_exit(msg)
        elif is_file and path.is_dir():
            msg = (f"Path '{path}' is an existing directory, not a file\n"
                   f"--overwrite will not overwrite a directory with a file")
            _communicate_error_and_exit(msg)
        # files in existing directory are not deleted with --overwrite
        elif not is_file and path.is_dir() and overwrite and any(path.iterdir()):
            msg = (f"Path '{path}' is an existing directory;\n"
                   f"--overwrite will not remove existing files, though they many be overwritten")
            _communicate_warning(msg)

        # file exists, can be overwritten, but is not allowed too
        if (not can_exist) and (not overwrite):
            msg = f"Path '{path}' already exists; use deli --overwrite <your-command> to overwrite"
            _communicate_error_and_exit(msg)

        # check for reserved file names on Windows
        if os.path.isreserved(path):
            msg = f"Path '{path}' uses a reserved name on Windows; please choose a different name"
            _communicate_error_and_exit(msg)

        # warning if path is a symlink
        if path.is_symlink():
            _communicate_warning(f"Path '{path}' is a symbolic link; original reference file will *not* be overwritten")

    # handle the extension check
    if is_file and extension:
        if len(path.suffixes) == is_compressed:
            msg = f"Output file '{path}' is missing a file extension {extension}"
            _communicate_warning(msg)
        else:
            suffix = path.suffixes[-1 - is_compressed]  # second extension from the end if compressed
            if suffix != extension:
                msg = f"Output file '{path}' has extension '{suffix}' but expected '{extension}'"
                _communicate_warning(msg)

    return path


def _open_text_file(path: Path):
    """
    Open a plain text file, handling gzip if needed

    Parameters
    ----------
    path : Path
        The path to the file to open.

    Returns
    -------
    file object
        The opened file object.
    """
    import gzip

    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    else:
        try:
            with open(path, "r") as f:
                # try reading a small chunk to check if it's a text file
                f.read(1024)
            return open(path, "r")
        except UnicodeDecodeError:
            msg = f"Cannot read file '{path}' as a text file; is this a valid text file?"
            logging.getLogger("deli").error(msg)
            click.echo(msg)


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
        _config_path = ctx.obj["deli_config"].location
        if _config_path.exists():
            _config = configparser.RawConfigParser()
            try:
                _config.read(os.path.normpath(_config_path))
                _config["deli.data"]["deli_data_dir"] = str(_path.resolve())
            except Exception as e:
                click.echo(
                    f"Failed to parse DELi config file at '{_config_path}'\nIs this a valid DELi config file?"
                    f"\nError details: {e}"
                )
                sys.exit(1)
            _config.write(open(_config_path, "w"), True)
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
@click.option("--out-dir", "-o", type=click.Path(), required=False, default="./", help="Output directory to save results to")
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
    Convert DNA sequences into DEL compound identities

    DECODE-FILE is the path to a YAML configure file describing the DEL selection,
    sequences, and decoding settings.

    DELi treats this process as "embarrassingly parallel"; if you were to split
    the input sequences across N separate processes, the output files can be trivially
    combined and will be identical to a single process with all sequences.


    NOTE: Some outputs are generated on the fly, stopping a job mid run can result
    in partial output files.

    See the docs for a detailed description of output files generated.
    See the `deli decode merge` command for more info on combining output files.

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
        logger.debug(f"closed file '{file}'")

    if failed_out_file is not None:
        failed_out_file.close()
        logger.debug(f"closed file '{failed_out_file}'")

    # handle writing statistics and report file
    statistics_file = out_dir_path / f"{prefix}_decode_statistics.json"
    decoder.decode_stats.to_file(statistics_file)
    logger.info(f"wrote decoding statistics to: '{statistics_file}'")

    # generate report if not skipped
    if not skip_report:
        from deli.decode.report import build_decoding_report

        report_file = out_dir_path / f"{prefix}_decode_report.html"
        build_decoding_report(stats=decoder.decode_stats, selection=selection, out_path=report_file)
        logger.info(f"wrote decoding report to: '{report_file}'")


@decode_group.command(name="collect")
@click.argument("decoded-reads", type=click.Path(exists=True), nargs=-1, required=True)
@click.option("--score-threshold", "-s", type=click.INT, required=False, default=100, help="Reject decoded compounds above this score")
@click.option("--count-threshold", "-c", type=click.INT, required=False, default=1, help="Minimum (dedup) count threshold for including compounds in output")
@click.option("--out-loc", "-o", type=click.Path(), required=False, default="./aggregated_decodes.json", help="Output location to save aggregated decodes to; will add .json suffix if missing")
@click.option("--compress", "-z", is_flag=True, help="Compress output with gzip")
@click.pass_context
def collect_decodes(ctx, decoded_reads, score_threshold, count_threshold, out_loc, compress):
    """
    Collect decoded reads from decoded TSV file(s) into a JSON count format

    This process requires that all decoded that are part of a given selection are provided.
    If they are not, count values can be overestimated in future jobs

    NOTE: this is a memory intensive operation; ensure sufficient memory is available

    JSON will be new line delimited (NDJSON), with each line representing a unique compound
    as a dictionary with the following fields:
    - library_id: str
        The library ID the compound belongs to
    - bb_ids: str
        The building block IDs of the compound, joined by commas
    - umi_counts: list[dict]
        A list of dictionaries, each with keys 'k' and 'c' representing the unique
        UMIs ('k') and their counts ('c') associated with this compound
    """
    import polars as pl

    logger = ctx.obj["logger"]

    # validate output directory
    out_path = Path(out_loc).absolute()
    if out_path.exists() and out_path.is_dir():
        msg = f"output location '{out_path}' is a directory; please provide a file path"
        logger.error(msg)
        click.echo(msg)
        sys.exit(1)
    if compress:
        if out_path.suffixes[-2:] != [".ndjson", ".gz"]:
            if out_path.suffixes[-1] != ".ndjson":
                out_path = out_path.with_suffix(".ndjson.gz")
            else:
                out_path = out_path.with_suffix(".gz")
    else:
        if out_path.suffix != ".ndjson":
            out_path = out_path.with_suffix(".ndjson")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    logger.debug(f"writing aggregated decoded sequences to: '{out_path}'")

    # validate input files
    for f in decoded_reads:
        f_path = Path(f).absolute()
        if not f_path.exists() or not f_path.is_file():
            msg = f"decoded reads file '{f_path}' does not exist or is not a file"
            logger.error(msg)
            click.echo(msg)
            sys.exit(1)
    logger.debug(f"aggregating decoded reads from {len(decoded_reads)} input files")

    if compress:
        logger.debug("using gzip compression for output JSON")

    # stream in all the decodes and aggregate
    (
        pl.scan_csv(decoded_reads, has_header=True, separator="\t", ignore_errors=True)
        .select(["library_id", "bb_ids", "umi", "overall_score"])
        .filter(pl.col("overall_score") <= score_threshold)
        .groupby(["library_id", "bb_ids"])
        .agg(pl.col("umi"))
        .filter(pl.col("umi").list.len() >= count_threshold)
        .with_columns(
            umi_counts=pl.col("umi").list.eval(pl.element().value_counts().struct.rename_fields(["k", "c"])),
        )
        .sink_ndjson(out_path, compression="gzip" if compress else "uncompressed")
    )
    logger.info(f"aggregated decoded sequences from {len(decoded_reads)} input files to {out_path}")


@decode_group.command(name="count")
@click.argument("aggregated-compounds", type=click.Path(exists=True), required=True)
@click.option("--out-loc", "-o", type=click.Path(), required=False, default=None, help="Path to output file")
@click.option("--cluster-umis", "-u", is_flag=True, help="Cluster UMIs to determine final count")
@click.option("--keep-raw-count", "-r", is_flag=True, help="Keep raw count in output")
@click.option("--keep-dedup-count", "-d", is_flag=True, help="Keep deduplicated count in output; ignored unless --cluster-umis is used")
@click.option("--output-type" "-t", type=click.Choice(["tsv", "gzip", "parquet", "avro"]), default="tsv", help="Output file type")
@click.pass_context
def count_compounds(ctx, aggregated_compounds, out_loc, cluster_umis, keep_raw_count, keep_dedup_count, output_type):
    """
    Count compounds from aggregated decoded compounds file

    AGGREGATED-COMPOUNDS is the path to the aggregated decoded compounds NDJSON file
    generated using `deli decode collect`

    NOTE: "gzip" output type will generate a gzip compressed TSV file. Parquet and AVRO files are written at the
    end using a temperary file
    """
    from deli.decode.count import raw_count, dedup_count, corrected_count

    # check conditional imports
    if output_type == "parquet":
        try:
            import pyarrow as pa # noqa: F401
            import pyarrow.parquet as pq # noqa: F401
        except ImportError:
            msg = "pyarrow is required for parquet output; please install with 'pip install pyarrow'"
            click.echo(msg)
            ctx.obj["logger"].error(msg)
            sys.exit(1)
    elif output_type == "avro":
        try:
            import fastavro # noqa: F401
        except ImportError:
            msg = "fastavro is required for avro output; please install with 'pip install fastavro'"
            click.echo(msg)
            ctx.obj["logger"].error(msg)
            sys.exit(1)

    logger = ctx.obj["logger"]

    # determine output type
    compress = False
    if output_type == "gzip":
        output_type = "tsv"
        compress = True
    elif output_type == "parquet":
        output_type = "parquet"
    elif output_type == "avro":
        output_type = "avro"
    elif output_type == "tsv":
        output_type = "tsv"
    else:
        msg = f"unknown output type '{output_type}'"
        logger.error(msg)
        click.echo(msg)
        sys.exit(1)

    # prepare output path
    if out_loc is None:
        out_loc = Path(f"./counted_compounds.{output_type}")
    out_loc = _standardize_and_validate_output_loc(
        out_loc,
        is_file=True,
        can_exist=True,
        overwrite=True,
        extension=f".{output_type}",
        is_compressed=compress,
    )

    # initialize writer
    if output_type in ["tsv", "gzip"]:
        if output_type == "tsv":
            out_file = open(out_loc, "w")
        else:
            out_file = gzip.open(out_loc, "wt")
        def _write_tsv(row_dict):
            out_file.write("\t".join([str(val) for val in row_dict.values()]).strip() + "\n")
        writer = _write_tsv
    elif output_type == "avro":
        schema = {
            'type': 'record',
            'name': 'CompoundCounts',
            'fields': [
                {'name': 'library_id', 'type': 'string'},
                {'name': 'bb_ids', 'type': 'string'},
                {'name': 'count', 'type': 'int'},
                {'name': 'raw_count', 'type': ['int', 'null']} if keep_raw_count else None,
                {'name': 'dedup_count', 'type': ['int', 'null']} if keep_dedup_count and not cluster_umis else None,
            ]
        }
        schema['fields'] = [field for field in schema['fields'] if field is not None]  # Remove None fields
        out_file = open(out_loc, "wb")
        writer = fastavro.writer(out_file, schema)
    elif output_type == "parquet":
        table_schema = pa.schema([
            pa.field("library_id", pa.string()),
            pa.field("bb_ids", pa.string()),
            pa.field("count", pa.int32()),
            pa.field("raw_count", pa.int32()) if keep_raw_count else None,
            pa.field("dedup_count", pa.int32()) if keep_dedup_count and not cluster_umis else None,
        ])
        table_schema = table_schema.remove_metadata()  # Remove None fields
        pq_writer = pq.ParquetWriter(out_loc, table_schema)
        def _write_parquet(batch):
            pa_batch = pa.RecordBatch.from_pylist(batch, schema=table_schema)
            pq_writer.write_batch(pa_batch)
        writer = _write_parquet
    else:
        msg = f"unknown output type '{output_type}'"
        logger.error(msg)
        click.echo(msg)
        sys.exit(1)

    _batch_size = 5000
    _batch = []
    _ticker = 0
    with _open_text_file(Path(aggregated_compounds)) as in_file:
        for line in tqdm(in_file, desc="Counting compounds"):
            _ticker += 1
            cpd_info = ast.literal_eval(line)

            # calculate counts
            raw_count = sum([count_struct["c"] for count_struct in cpd_info["umi_counts"]])
            dedup_count = len(cpd_info["umi_counts"])
            count = dedup_count  # default to dedup count if not clustering
            if cluster_umis:
                count = corrected_count(cpd_info["umi_counts"])

            # prepare output info
            row = {
                "library_id": cpd_info["library_id"],
                "bb_ids": cpd_info["bb_ids"],
                "count": count,
            }
            if keep_raw_count:
                row["raw_count"] = raw_count
            if keep_dedup_count and not cluster_umis:
                row["dedup_count"] = dedup_count
            _batch.append(row)

            if len(_batch) >= _batch_size:
                logger.info(f"writing {len(_batch)} records to output")
                writer(_batch)
                _batch = []

    if _batch:
        writer(_batch)

    # close out files
    if output_type == "avro":
        writer.flush()
        logger.debug("flushed avro writer")
        out_file.close()
    elif output_type == "parquet":
        pq_writer.fl
        pq_writer.close()
    else:
        out_file.close()
    logger.debug("closed output file")

    logger.info(f"wrote f{_ticker} compounds to: '{out_loc}'")


@decode_group.command(name="report")
@click.argument("decode_stats_file", nargs=-1, type=click.Path(exists=True), required=True)
@click.argument("decode_file", type=click.Path(exists=True), required=False, default=None)
@click.option("--out-loc", "-o", type=click.Path(), required=False, default="./decode_report.html", help="Output location to save report to")
@click.pass_context
def generate_report(ctx, decode_stats_file, out_loc, selection):
    """
    Generate an HTML decoding report from decoding statistics file(s)

    DECODE_STATS_FILE is the path(s) to the decoding statistics JSON files to generate a report with .

    If multiple statistic files are provided, will merge them into a single report.
    This is useful for dealing with parallel decoding runs for the same selection.

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


@cli.command(name="cubify")
@click.argument("aggregated-compounds", type=click.Path(exists=True), required=True)
@click.option("--output", "-o", type=click.Path(), required=False, default="./cube", help="Location to save results to")
@click.option("--cube-format", "-f", type=click.STRING, default="lbc", help="Cube format string; see DELi documentation for details")
@click.option("--overwrite", "-w", is_flag=True, help="Overwrite existing cube file if it exists")
@click.option("--corrected-count-threshold", "-E", type=click.INT, default=0, help="Threshold to include compounds in Cube based on error-corrected count")
@click.option("--dedup-count-threshold", "-D", type=click.INT, default=0, help="Threshold to include compounds in Cube based on deduped count")
@click.option( "--raw-count-threshold","-R", type=click.INT, default=0, help="Threshold to include compounds in Cube based on raw count")
@click.option("--use-tqdm", "-t", is_flag=True, help="Show tqdm progress bar")
@click.option("--selection", "-s", type=click.Path(exists=True), required=False, default=None, help="Path to DELi selection file to use for compound ID generation")
@click.pass_context
def run_cubify(ctx, aggregated_compounds, output, cube_format, overwrite, corrected_count_threshold, dedup_count_threshold, raw_count_threshold, use_tqdm, selection):
    """
    Covert aggregated decoded compounds (in NDJSON format) into a formated "Cube" file.

    Cube formats can include:
    - i: compound id
    - l: library id
    - s: compound SMILES
    - b: building block ids (one column per cycle)
    - B: building block SMILES (one column per cycle)
    - e: error-corrected count
    - d: deduplicated count
    - r: raw count



    For example; to just include the compound ID and its error-corrected count, use the format string "ie".

    NOTE: DELi considers Compounds IDs to be irreversible; that is, the compound ID cannot be used to
    retrieve the original building block IDs or library information. Be aware that excluding library ID
    and building block IDs from the Cube file will prevent DELi from being able to map back to the original
    compound / DEL information in any future uses of the Cube file.

    See the DELi documentation for more information on how to specify Cube files.
    """
    pass


@cli.command(name="enumerate")
@click.argument("library_file", type=click.Path(exists=True), required=True)
@click.option("--out_path", "-o", type=click.Path(), required=False, default="", help="Output CSV file path")
@click.option("--use-tqdm", "-t", is_flag=True, help="Enable TQDM progress bar")
@click.option("--fail-on-error", "-f", is_flag=True, help="Fail on first error during enumeration")
@click.option("--drop-failed", "-d", is_flag=True, help="Drop compounds with failed enumerations")
def enumerate_(library_file, out_path, use_tqdm, fail_on_error, drop_failed):
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
        use_tqdm=use_tqdm,
        fail_on_error=fail_on_error,
        drop_failed=drop_failed,
    )
