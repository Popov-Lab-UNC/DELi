"""command line functions for deli"""

import ast
import configparser
import datetime
import os
import sys
from pathlib import Path
import logging

# Suppress RDKit warnings at module level
try:
    from rdkit import RDLogger
    rd_logger = RDLogger.logger()
    rd_logger.setLevel(RDLogger.ERROR)
except ImportError:
    pass

import click
import pandas as pd
import yaml

from deli.analysis.analysis_report_gen import generate_report
from deli.analysis.cube_class import DELi_Cube
from deli.configure import (
    DeliDataDirError,
    get_deli_config,
    init_deli_config,
    init_deli_data_directory,
    set_deli_data_dir,
    validate_deli_data_dir,
)
from deli.decode import (
    DecodingRunner,
)
from deli.dels.library import Library


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
@click.argument("path", required=False, default=None)
@click.option("--overwrite", "-o", is_flag=True, help="Overwrite any existing config directory")
def click_init_deli_config(path, overwrite):
    """
    Create a default DELi configuration file

    PATH is the path to the deli config directory to initialize.
    If not provided, defaults to ~/.deli
    """
    _path = Path(path) if path is not None else Path.home() / ".deli"
    try:
        init_deli_config(_path, fail_on_exist=not overwrite)
    except FileExistsError:
        print(
            f"'{_path}' already exists; config file not created\n"
            f"Use `deli config init --overwrite` to overwrite existing config file"
        )
        sys.exit(1)


@cli.group()
def data():
    """Group for config related commands"""
    pass


@data.command(name="init")
@click.argument(
    "path",
    type=click.Path(),
    required=False,
    default="./deli_data",
)
@click.option(
    "--fix-missing", "-f", is_flag=True, help="Fix a deli data directory missing sub-directories"
)
@click.option("--overwrite", "-o", is_flag=True, help="Overwrite any existing data directories")
def click_init_deli_data_dir(path, fix_missing, overwrite):
    """
    Initialize the configuration directory

    PATH is the path to the deli data directory to initialize.

    NOTE: fix-missing will not overwrite existing sub-directories, while overwrite will.
    'overwrite' will also add any missing sub-directories, but also replace existing ones.
    """
    _path = Path(path).resolve()
    try:
        init_deli_data_directory(
            _path, fail_on_exist=not (fix_missing or overwrite), overwrite=overwrite
        )
    except FileExistsError:
        print(
            f"'{_path}' already exists\n"
            f"you can create a new DELi data directory with this name using "
            f"'deli data init --overwrite {_path}'"
        )
        sys.exit(1)
    except NotADirectoryError:
        print(
            f"'{_path}' is not a directory; cannot be fixed\n"
            f"you can create a new DELi data directory with this name using "
            f"'deli data init --overwrite {_path}'"
        )
        sys.exit(1)


@data.command(name="which")
def click_which_deli_data_dir():
    """
    Print the current DELi data directory

    If the DELi data directory is not set, will print a message
    and exit with a non-zero status code.
    """
    try:
        _path = get_deli_config().deli_data_dir
        print(f"Current DELi data directory: {_path}")
    except DeliDataDirError:
        print("DELi data directory is not set")
        sys.exit(1)


@data.command(name="set")
@click.argument("path", type=click.Path(), required=True)
@click.option(
    "--update-config",
    "-u",
    is_flag=True,
    help="Update the DELi config to use this data directory as default",
)
def click_set_deli_data_dir(path, update_config):
    """
    Set the DELi data directory to use for decoding

    PATH is the path to the deli data directory to set.

    NOTE: if not using --update-config, you will need to set the DELI_DATA_DIR environment variable
    manually; the command required will be printed after running.
    """
    _path = Path(path).resolve()
    try:
        validate_deli_data_dir(_path)
    except FileNotFoundError:
        print(
            f"directory at '{_path}' does not exist\n"
            f"you can create a new DELi data directory using 'deli data init {_path}'"
        )
        sys.exit(1)
    except NotADirectoryError:
        print(
            f"'{_path}' is not a directory\n"
            f"you can create a new DELi data directory with this name using "
            f"'deli data init --overwrite {_path}'"
        )
        sys.exit(1)
    except DeliDataDirError:
        print(
            f"DELi data directory '{_path}' is missing required sub-directories\n"
            f"use 'deli data init --fix-missing {_path}' to add missing sub-directories"
        )
        sys.exit(1)

    # update the config with the new data directory if requested
    if update_config:
        _config_path = Path.home() / ".deli"
        if _config_path.exists():
            _config = configparser.RawConfigParser()
            _config.read(os.path.normpath(_config_path))
            _config["deli.data"]["deli_data_dir"] = str(_path.resolve())
            _config.write(open(_config_path, "w"), True)
        else:
            print(f"cannot find DELi config file at {_config_path}")
            sys.exit(1)
    else:
        import platform

        if platform.system() == "Windows":
            _command = f"set DELI_DATA_DIR={_path}"
        else:
            _command = f"export DELI_DATA_DIR={_path}"
        print(f"to set the deli data directory, run '{_command}'")


@cli.command(name="decode")
@click.argument("decode-file", type=click.Path(exists=True), required=True)
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
    help="Ignore the fastq observed_seq files in the decode file",
)
@click.option(
    "--prefix", "-p", type=click.STRING, required=False, default="", help="Prefix for output files"
)
@click.option("--tqdm", "-t", is_flag=True, help="Show tqdm progress")
@click.option(
    "--save-failed", is_flag=True, help="Save failed decoding results to a separate file"
)
@click.option("--save-counter", is_flag=True, help="Save raw decoding counters as JSON file")
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
    decode_file,
    fastq_files,
    ignore_decode_seqs,
    out_dir,
    prefix,
    tqdm,
    save_failed,
    save_counter,
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
    out_dir = os.path.abspath(out_dir)

    if deli_data_dir is not None:
        set_deli_data_dir(deli_data_dir)

    runner = DecodingRunner.from_file(
        decode_file,
        fastq_files,
        ignore_decode_seqs=ignore_decode_seqs,
        debug=debug,
        disable_logging=disable_logging,
    )
    save_failed_to = out_dir if save_failed else None

    os.makedirs(out_dir, exist_ok=True)
    if prefix is None or prefix == "":
        prefix = runner.selection.selection_id

    results = runner.run(save_failed_to=save_failed_to, use_tqdm=tqdm)

    runner.logger.info(f"Saving outputs to {out_dir}")

    statistics_out_path = os.path.join(out_dir, f"{prefix}_decode_statistics.json")
    runner.logger.debug(f"Saving decode statistics to {statistics_out_path}")
    results.write_decode_statistics(statistics_out_path)

    cube_out_path = os.path.join(out_dir, f"{prefix}_cube.csv")
    runner.logger.debug(f"Saving cube to {cube_out_path}")
    results.write_cube(cube_out_path)

    if save_counter:
        counter_out_path = os.path.join(out_dir, f"{prefix}_counter.json.gz")
        runner.logger.debug(f"Saving counter to {counter_out_path}")
        runner.degen.to_json(counter_out_path, compress=True)

    if not skip_report:
        report_out_path = os.path.join(out_dir, f"{prefix}_report.html")
        runner.logger.debug(f"Saving cube to {report_out_path}")
        results.write_decode_report(report_out_path)


@cli.command(name="enumerate")
@click.argument("library_file", type=click.Path(exists=True), required=True)
@click.option(
    "--out_path", "-o", type=click.Path(), required=False, default="", help="Output CSV file path"
)
@click.option("--tqdm", "-t", is_flag=True, help="Enable TQDM progress bar")
def enumerate_(library_file, out_path, tqdm):
    """
    Enumerates compounds from a given library

    If out_path is not provided, will save to the current working directory
    as a CSV file named <library_id>_enumerated.csv

    LIBRARY_FILE is the path to a DELi library file to enumerate.
    """
    library_id = os.path.basename(library_file).split(".")[0]
    output_file = (
        out_path if out_path != "" else os.path.join(os.getcwd(), f"{library_id}_enumerated.csv")
    )

    _start = datetime.datetime.now()

    enumerator = Library.load(library_file)
    enumerator.enumerate_to_file(output_file, separator=",", use_tqdm=tqdm)


@cli.command(name="analyze")
@click.option(
    "--config",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Path to the YAML config file",
)
def analyze(config):
    """
    Perform DEL analysis based on the provided YAML configuration file.

    CONFIG is the path to the YAML configuration file.
    """
    print("Analysis started with config:", config)
    
    def create_output_dir(output_dir):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        return output_dir

    def create_dated_output_dir(base_output_dir, name_suffix=""):
        output_dir = create_output_dir(base_output_dir)
        date_str = datetime.datetime.now().strftime("%Y%m%d")
        dated_output_dir = os.path.join(output_dir, f"{date_str}_{name_suffix}")
        if not os.path.exists(dated_output_dir):
            os.makedirs(dated_output_dir)
        return dated_output_dir

    def load_config(config_path):
        with open(config_path, "r") as file:
            return yaml.safe_load(file)

    def read_dict_from_file(file_path):
        with open(file_path, "r") as f:
            content = f.read()
            return ast.literal_eval(content)

    config = load_config(config)
    output_dir_base = config["general"].get("output_dir", "output")
    output_dir = create_dated_output_dir(output_dir_base, "analysis")
    
    # Set up log file to capture all output
    log_file = os.path.join(output_dir, "deli.log")
    log_file_handle = open(log_file, 'w')
    
    # Redirect stdout and stderr to log file
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    sys.stdout = log_file_handle
    sys.stderr = log_file_handle
    
    print("Analysis started with config:", config)
    print(f"Logging all output to: {log_file}")

    data_file = config["general"].get("data", "")
    if not data_file:
        raise ValueError("Data file must be specified in the YAML config.")
    df = pd.read_csv(data_file)

    id_col = config["general"].get("ID_col", df.columns[0])
    indexes = ast.literal_eval(str(config.get("indexes", {})))
    control_cols = ast.literal_eval(str(config.get("control_cols", {})))
    raw_indexes = ast.literal_eval(str(config.get("raw_indexes", {})))

    indexes_path = os.path.join(output_dir, "indexes.txt")
    control_cols_path = os.path.join(output_dir, "control_cols.txt")
    raw_indexes_path = os.path.join(output_dir, "raw_indexes.txt")

    with open(indexes_path, "w") as file:
        file.write(str(indexes))
    with open(control_cols_path, "w") as file:
        file.write(str(control_cols))
    with open(raw_indexes_path, "w") as file:
        file.write(str(raw_indexes))

    indexes = read_dict_from_file(indexes_path)
    control_cols = read_dict_from_file(control_cols_path)
    raw_indexes = read_dict_from_file(raw_indexes_path)
    if not indexes:
        indexes = raw_indexes

    cube = DELi_Cube(
        df,
        id_col,
        indexes,
        control_cols,
        int(config["general"].get("lib_size", 0)),
        raw_indexes,
    )

    if "flags" in config:
        flags = config["flags"]
        if flags.get("SD_min", False):
            nsc_max_dict, sd_min_dict, sampling_depth_dict = cube.SD_min()
        if flags.get("NSC_values", False):
            cube.NSC_values()
        if flags.get("MLE", False):
            cube.maximum_likelihood_enrichment_ratio()
        if flags.get("Z_score", False):
            cube.z_score()
        if flags.get("z_score_log_data", False):
            cube.z_score_log_data()
        if flags.get("disynthon_data", False):
            disynthon_data, disynth_exp_dict = cube.disynthonize()
            cube.data = disynthon_data
        if flags.get("polyO", False):
            cube.PolyO()
        if flags.get("top_disynthons", False):
            comparison_type = flags["top_disynthons"].get("comparison", "control")
            exp_name = flags["top_disynthons"].get("exp_name", "None")
            exp2_name = flags["top_disynthons"].get("exp2_name", "None")
            control_name = flags["top_disynthons"].get("control_name", "None")
            top_count = int(flags["top_disynthons"].get("top_count", 10))
            comparison_metric = flags["top_disynthons"].get("comparison_metric", "avg")
            top_disynthons_dir = create_output_dir(os.path.join(output_dir, "top_disynthons"))
            cube.get_top_disynthons(
                disynthon_data=disynthon_data,
                exp_name1=exp_name,
                comparison_type=comparison_type,
                exp_name2=exp2_name,
                control_name=control_name,
                comparison_metric=comparison_metric,
                top_count=top_count,
                output_dir=top_disynthons_dir,
            )
        if flags.get("trisynthon_overlap", False):
            trisynthon_dir = create_output_dir(os.path.join(output_dir, "trisynthon"))
            cube.trisynthon_overlap(output_dir=trisynthon_dir)
        if flags.get("disynthon_overlap", False):
            disynthon_dir = create_output_dir(os.path.join(output_dir, "disynthon"))
            cube.disynthon_overlap(
                output_dir=disynthon_dir,
                disynthon_data=disynthon_data,
                disynth_exp_dict=disynth_exp_dict,
                threshold=int(flags.get("disynthon_threshold", 20)),
            )
        if flags.get("normalized_data", False):
            cube.normalize()
        if flags.get("simple_spotfire_version", False):
            spotfire = cube.simple_spotfire_version()
            today_date = datetime.datetime.now().strftime("%Y%m%d")
            spotfire.to_csv(os.path.join(output_dir, f"spotfire_{today_date}.csv"), index=False)
        if flags.get("ml_fingerprints_to_RF_reg", False):
            ml_fingerprints_to_RF_dir = create_output_dir(
                os.path.join(output_dir, "ml_fingerprints_to_RF")
            )
            cube.ml_fingerprints_to_RF(output_dir=ml_fingerprints_to_RF_dir)
        if flags.get("ml_fingerprints_to_RF_clf", False):
            ml_fingerprints_to_RF_clf_dir = create_output_dir(
                os.path.join(output_dir, "ml_fingerprints_to_clf")
            )
            cube.ml_fingerprints_to_classifier(
                output_dir=ml_fingerprints_to_RF_clf_dir,
                threshold=int(flags.get("clf_thresh", 10)),
            )
        if flags.get("gnn_classifier", False):
            gnn_dir = create_output_dir(os.path.join(output_dir, "gnn"))
            cube.gnn_classifier(
                output_dir=gnn_dir,
                threshold=int(flags.get("gnn_threshold", 10)),
                arch=flags.get("gnn_arch", "GAT"),
            )
        if "top_hits" in flags:
            top_hits_dir = create_output_dir(os.path.join(output_dir, "top_hits"))
            cube.top_n_compounds(
                int(flags["top_hits"]),
                flags.get("top_hits_metric", "sum"),
                output_dir=top_hits_dir,
            )
        if flags.get("monosynthon_chemical_space", False):
            cube.monosynthon_chemical_space(output_dir=output_dir)
        if flags.get("report", False):
            nsc_max_dict = nsc_max_dict if flags.get("SD_min", False) else None
            sd_min_dict = sd_min_dict if flags.get("SD_min", False) else None
            sampling_depth_dict = sampling_depth_dict if flags.get("SD_min", False) else None
            generate_report(
                output_dir_base,
                indexes,
                control_cols,
                nsc_max_dict,
                sd_min_dict,
                sampling_depth_dict,
            )
            today_date = datetime.datetime.now().strftime("%Y%m%d")
            cube.data.to_csv(
                os.path.join(output_dir_base, f"cube_data_{today_date}.csv"),
                index=False,
            )
    
    # Cleanup: restore stdout/stderr and close log file
    sys.stdout = original_stdout
    sys.stderr = original_stderr
    log_file_handle.close()
    
    print(f"Analysis completed! Check the log file at: {log_file}")
    print(f"Results saved to: {output_dir}")
