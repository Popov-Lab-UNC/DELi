"""command line functions for deli"""

import configparser
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
            if "sequence_files" in e_ or "libraries" in e_:
                logger.debug("selection file lacks sequencing info")
            else:
                raise e_
        try:
            return DELSelection.from_yaml(selection)
        except (DeliDataNotFound, FileNotFoundError):
            logger.debug(f"failed to load libraries from selection file '{selection}'")
        except KeyError as e_:
            if "libraries" not in e_:
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
    skip_report,
):
    """
    Decode DNA sequences from a DEL selection

    DECODE is the path to a YAML file describing the decoding run settings.

    See the docs for a detailed description of output files generated.

    """
    from deli.decode.base import FailedDecodeAttempt
    from deli.decode.decoder import DecodedDELCompound, DecodedToolCompound, DecodingSettings, SelectionDecoder
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

    # handle output file locations
    if split_by_lib:
        decode_out_dir = out_dir_path / "decodes_by_library"
        decode_out_dir.mkdir(parents=False, exist_ok=False)
        decoded_out_paths = {
            library.library_id: open(decode_out_dir / f"{prefix}_{library.library_id}.tsv", "w")
            for library in selection.library_collection.libraries
        }
        logger.info("saving decoded sequences per library to 'decode_out_dir'")
    else:
        decoded_out_path = out_dir_path / f"{prefix}_decoded.tsv"
        decoded_out_paths = {"all": open(decoded_out_path, "w")}
        logger.info(f"writing decoded sequences to: '{decoded_out_path}'")

    # write headers
    for file in decoded_out_paths.values():
        header = ""
        if save_fastq_info:
            header += "fastq_file\tread_name\t"
        header += "library_id\tbb_ids\tumi\tscore"
        file.write(f"{header}\n")

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
                out_file = decoded_out_paths[decoded_read.get_library_id()]
            else:
                out_file = decoded_out_paths["all"]

            line = ""
            if save_fastq_info:
                line += f"{sequence_file}\t{sequence_record.name}\t"
            if isinstance(decoded_read, DecodedDELCompound):
                line += f"{decoded_read.get_library_id()}\t{','.join(decoded_read.building_block_ids)}\t"
            elif isinstance(decoded_read, DecodedToolCompound):
                line += f"{decoded_read.tool_compound.compound_id}\tnull\t"
            line += f"{decoded_read.umi_str()}\t{decoded_read.get_score()}"

            out_file.write(f"{line}\n")

        if (i + 1) % 100000 == 0:
            logger.debug(f"decoded {i + 1} sequences...")

    logger.info(f"finished decoding for selection '{selection.selection_id}'")

    for file in decoded_out_paths.values():
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
@click.argument("decode_files", type=click.Path(exists=True), nargs=-1, required=True)
@click.option("--out-loc", "-o", type=click.Path(), required=False, default="./", help="Location to save results to")
@click.option("--use-tqdm", "-t", is_flag=True, help="Show tqdm progress")
@click.option("--counts", "-c", is_flag=True, help="Return counts only, no UMI information (CSV output)")
@click.option("--raw-counts", "-r", is_flag=True, help="Include raw counts in output (only valid with --counts flag)")
@click.option(
    "--skip-missing-umi", "-m", is_flag=True, help="Skip sequences without/missing UMIs rather than counting them"
)
@click.option(
    "--score-threshold",
    "-q",
    type=click.INT,
    required=False,
    default=100,
    help="Minimum score threshold for including decoded sequences in degeneration",
)
@click.pass_context
def run_degen(ctx, decode_files, out_loc, use_tqdm, counts, raw_counts, skip_missing_umi, score_threshold):
    """
    Degenerate decoded sequences into compound counts based on UMIs

    DECODE_FILES are the paths to the decoded TSV files to degenerate.

    If --counts flag is used, the output will be a TSV file with compound counts only.
    Otherwise, the output will be a JSON file with full UMI count information.

    Note: Degeneration outputs must be in the UMI JSON format to be merged, using the
    count CSV format will prevent future merging of degeneration outputs.

    See the docs for a detailed description of output files generated.

    """
    from deli.decode.degen import DegenSettings, SelectionDegenerator

    logger = ctx.obj["logger"]

    suffix = ".csv" if counts else ".json"

    # validate output directory
    out_dir_path = Path(out_loc).absolute()
    if out_dir_path.suffix == "":
        out_dir_path = out_dir_path / f"degen_results{suffix}"
    elif out_dir_path.suffix != suffix:
        out_dir_path = out_dir_path.with_suffix(suffix)
    out_dir_path.mkdir(parents=True, exist_ok=True)
    logger.debug(f"writing degenerate sequences to '{out_dir_path}'")

    settings = DegenSettings(ignore_missing_umi=skip_missing_umi)
    degenerator = SelectionDegenerator(settings)

    pg = tqdm(desc="degenerating decoded sequences", total=len(decode_files), disable=not use_tqdm)
    _ticker = 0
    for decode_file in decode_files:
        logger.info(f"processing '{decode_file}'")
        with open(decode_file, "r") as in_file:
            header = in_file.readline().strip().split("\t")
            bb_ids_index = header.index("bb_ids")
            umi_index = header.index("umi")
            lib_id_index = header.index("lib_id")
            score_index = header.index("score")

            for i, decode_line in enumerate(in_file):
                splits = decode_line.strip().split("\t")

                score = int(splits[score_index])
                if score > score_threshold:
                    logger.debug(f"skipping decoded sequence with score {score} below threshold {score_threshold}")
                    continue

                umi = splits[umi_index]
                if skip_missing_umi and ((umi == "null") or (umi == "")):
                    logger.debug("skipping decoded sequence with missing UMI")
                    continue

                lib_id = splits[lib_id_index]
                bb_ids = splits[bb_ids_index].split(",") if splits[bb_ids_index] != "null" else []
                del_id = tuple([lib_id] + bb_ids)

                degenerator.degen_decoded_compound((del_id, umi))
                pg.update(1)
                _ticker += 1

                if (i + 1) % 100000 == 0:
                    logger.debug(f"processed {i + 1} decoded sequences from file '{decode_file}'")
    pg.close()
    logger.info(f"degeneration completed for {_ticker} decoded sequences from {len(decode_files)} files")

    if counts:
        logger.debug("writing output in counts TSV format")

        with open(out_dir_path) as out_file:
            header = "library_id\tbb_ids\tcount"
            if raw_counts:
                header += "\traw_count"
            out_file.write(f"{header}\n")

            for lib_id, lib_counter in degenerator.counter.items():
                for bb_id_tuple, counter in lib_counter.counter.items():
                    line = f"{lib_id}\t{','.join(bb_id_tuple)}\t{counter.degen_count}"
                    if raw_counts:
                        line += f"\t{counter.raw_count}"
                    out_file.write(f"{line}\n")
    else:
        logger.debug("writing output in full UMI JSON format")
        degenerator.write_json(out_dir_path)

    logger.info(f"wrote degenerate sequences to '{out_dir_path}'")


@degen_group.command(name="merge")
@click.argument("degen_files", type=click.Path(exists=True), nargs=-1, required=True)
@click.option(
    "--out_loc",
    "-o",
    type=click.Path(),
    required=False,
    default="./",
    help="Location to save merged degeneration file to",
)
@click.option("--counts", "-c", is_flag=True, help="Output merged file in counts TSV format instead of UMI JSON format")
@click.option("--raw-counts", "-r", is_flag=True, help="Include raw counts in output (only valid with --counts flag)")
@click.pass_context
def merge_degen(ctx, degen_files, out_loc, counts, raw_counts):
    """
    Merge multiple degeneration output files into a single file

    DEGEN_FILES are the paths to the degeneration output files to merge.

    All input files must be in the UMI JSON format to be merged.

    The merged output will be saved to 'merged_degen_<timestamp>.json' in the current working directory.

    """
    from deli.decode.degen import SelectionDegenerator, load_degenerator_from_json

    logger = ctx.obj["logger"]

    suffix = ".csv" if counts else ".json"

    # validate output directory
    out_loc_path = Path(out_loc).absolute()
    if out_loc_path.suffix == "":
        out_loc_path = out_loc_path / f"degen_results{suffix}"
    elif out_loc_path.suffix != suffix:
        out_loc_path = out_loc_path.with_suffix(suffix)
    out_loc_path.mkdir(parents=True, exist_ok=True)
    logger.debug(f"writing degenerate sequences to '{out_loc_path}'")

    merged_degen = SelectionDegenerator()

    logger.info(f"merging {len(degen_files)} degeneration files")
    for degen_file in degen_files:
        logger.debug(f"loading degeneration file '{degen_file}'")
        try:
            degenerator = load_degenerator_from_json(degen_file)
        except FileNotFoundError:
            logger.error(f"degeneration file '{degen_file}' not found; skipping")
            click.echo("cannot find degeneration file '{degen_file}'; skipped")
            continue
        except MemoryError:
            logger.exception(f"ran out of memory loading degeneration file '{degen_file}'; dumping and aborting merge")
            click.echo(
                f"ran out of memory loading degeneration file '{degen_file}' "
                f"dumping current merge and aborting the rest...\ntry running with more memory resources"
            )
            merged_degen.write_json(out_loc_path)
            logger.info(f"wrote partial merged degeneration data to '{out_loc_path}'")
            sys.exit(1)
        except Exception as e:
            logger.exception(f"failed to load degeneration file '{degen_file}': {e}'")
            click.echo(
                f"failed to parse degeneration file {degen_file}; is this a valid DELi degeneration output file?"
            )
            sys.exit(1)
        merged_degen = merged_degen.merge_degenerator(degenerator)
        logger.debug(f"merged degeneration data from '{degen_file}'")

    logger.info(f"merged degeneration data from {len(degen_files)} files")

    if counts:
        logger.debug("writing merged output in counts TSV format")

        with open(out_loc_path, "w") as out_file:
            header = "library_id\tbb_ids\tcount"
            if raw_counts:
                header += "\traw_count"
            out_file.write(f"{header}\n")

            for lib_id, lib_counter in merged_degen.counter.items():
                for bb_id_tuple, counter in lib_counter.counter.items():
                    line = f"{lib_id}\t{','.join(bb_id_tuple)}\t{counter.degen_count}"
                    if raw_counts:
                        line += f"\t{counter.raw_count}"
                    out_file.write(f"{line}\n")
    else:
        logger.debug("writing merged output in full UMI JSON format")
        merged_degen.write_json(out_loc_path)
    logger.info(f"wrote merged degeneration data to '{out_loc_path}'")


@degen_group.command(name="count")
@click.argument("degen_file", type=click.Path(exists=True, dir_okay=False), required=True)
@click.option("--raw-counts", "-r", type=click.STRING, required=True, help="Include raw counts in output")
@click.pass_context
def count_degen(ctx, degen_file, raw_counts):
    """
    Count the total number of unique compounds in a degeneration output file

    DEGEN_FILE is the path to the degeneration output file to count unique compounds from.

    """
    from deli.decode.degen import load_degenerator_from_json

    logger = ctx.obj["logger"]

    logger.debug(f"loading degeneration file '{degen_file}'")
    try:
        degenerator = load_degenerator_from_json(degen_file)
    except MemoryError:
        logger.exception(f"ran out of memory loading degeneration file '{degen_file}'")
        click.echo(
            f"ran out of memory loading degeneration file '{degen_file}'; try running with more memory resources"
        )
        sys.exit(1)
    except Exception as e:
        logger.exception(f"failed to load degeneration file '{degen_file}': {e}'")
        click.echo(f"failed to parse degeneration file {degen_file}; is this a valid DELi degeneration output file?")
        sys.exit(1)

    cur_degen_path = Path(degen_file)
    out_loc_path = (cur_degen_path.parent / cur_degen_path.stem).with_suffix(".tsv")
    logger.debug(f"writing output to '{out_loc_path}'")

    with open(out_loc_path, "w") as out_file:
        header = "library_id\tbb_ids\tcount"
        if raw_counts:
            header += "\traw_count"
        out_file.write(f"{header}\n")

        _ticker = 0
        for lib_id, lib_counter in degenerator.counter.items():
            for bb_id_tuple, counter in lib_counter.counter.items():
                line = f"{lib_id}\t{','.join(bb_id_tuple)}\t{counter.degen_count}"
                if raw_counts:
                    line += f"\t{counter.raw_count}"
                out_file.write(f"{line}\n")
                _ticker += 1

    logger.info(f"counted {_ticker} unique compounds in degeneration file '{degen_file}'")
    logger.info(f"wrote compound counts to '{out_loc_path}'")


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
