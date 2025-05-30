"""command line functions for deli"""

import configparser
import datetime
import enum
import os
import pickle
from collections import defaultdict
from pathlib import Path

import click

from deli.configure import (
    init_deli_config,
    init_deli_data_directory,
    set_deli_data_dir,
    validate_deli_data_dir,
)
from deli.decode import (
    DecodeStatistics,
    DecodingRunner,
    DecodingRunnerResults,
    DELCollectionCounter,
    build_decoding_report,
)
from deli.dels import Selection
from deli.dels.enumerator import DELEnumerator


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
@click.argument("path", type=click.Path(exists=False), required=False, default=None)
@click.option("--overwrite", "-o", is_flag=True, help="Overwrite any existing config directory")
def init_config(path, overwrite):
    """
    Initialize the configuration directory

    PATH is the path to the deli config directory to initialize.
    If not provided, defaults to ~/.deli/.deli
    """
    init_deli_config(path, fail_on_exist=not overwrite)


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
@click.option("--fix-missing", "-f", help="Fix a deli data directory missing sub-directories")
@click.option("--overwrite", "-o", is_flag=True, help="Overwrite any existing data directories")
def init_deli_data_dir(path, fix_missing, overwrite):
    """
    Initialize the configuration directory

    PATH is the path to the deli data directory to initialize.
    """
    init_deli_data_directory(
        path, fail_on_exist=not (fix_missing or overwrite), overwrite=overwrite
    )


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
        _config_path = Path.home() / ".deli"
        if _config_path.exists():
            _config = configparser.RawConfigParser()
            _config.read(os.path.normpath(_config_path))
            _config["deli.data"]["deli_data_dir"] = str(path.resolve())
            _config.write(open(_config_path, "w"), True)
        else:
            raise FileNotFoundError("cannot find DELi config file at {}")


@cli.group()
def decode():
    """Decode command group"""
    pass


@decode.command(name="run")
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
    help="Ignore the fastq sequence files in the decode file",
)
@click.option(
    "--prefix", "-p", type=click.STRING, required=False, default="", help="Prefix for output files"
)
@click.option("--tqdm", "-t", is_flag=True, help="Show tqdm progress")
@click.option(
    "--save-failed", is_flag=True, help="Save failed decoding results to a separate file"
)
@click.option(
    "--save-degen", is_flag=True, help="Save decoding degen counters instead of cube files"
)
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
    save_degen,
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

    if save_degen:
        degen_out_path = os.path.join(out_dir, f"{prefix}_counters.pkl")
        runner.logger.debug(f"Saving cube to {degen_out_path}")
        pickle.dump(results.degen, open(degen_out_path, "wb"))
    else:
        cube_out_path = os.path.join(out_dir, f"{prefix}_cube.csv")
        runner.logger.debug(f"Saving cube to {cube_out_path}")
        results.write_cube(cube_out_path)

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

    enumerator = DELEnumerator.load(library_file)
    enumerator.enumerate_to_csv_file(output_file, use_tqdm=tqdm)


@decode.group()
def statistics():
    """Decode statistics command group"""
    pass


@statistics.command(name="merge")
@click.argument("statistic-files", type=click.Path(exists=True), required=True, nargs=-1)
@click.option(
    "--counter-file",
    type=click.Path(exists=True),
    required=False,
    help="counter file to extract corrected degen counts from",
)
@click.option(
    "-o",
    "--out-path",
    type=click.STRING,
    required=False,
    default="merged_decoding_statistics.json",
    help="location to save merged statistics file",
)
def merge_statistics_file(statistic_files, counter_file, out_path):
    """
    Merge multiple decode statistics files into one

    WARNING: Merging statistics files from different decode runs that use
    UMI corrected counts will result in incorrect counts, as the UMIs are
    not stored in the statistics file. In this case, you should use
    the --save-degen option during decoding, merge the counter files
    and then pass the merged counter file to this function.

    STATISTICS is a list of paths to decode statistics files to merge.
    """
    merged_statistics = DecodeStatistics.from_file(statistic_files[0])
    for statistic_file in statistic_files[1:]:
        loaded_statistics = DecodeStatistics.from_file(statistic_file)
        merged_statistics += loaded_statistics

    # if a counter file is provided, use it to set the degen count
    if counter_file:
        loaded_counter: DELCollectionCounter = pickle.load(open(counter_file, "rb"))
        _updated_num_seqs_degen_per_lib = {
            lib_id: sum([del_counter.get_degen_count() for del_counter in lib_counters.values()])
            for lib_id, lib_counters in loaded_counter.del_counter.items()
        }

        merged_statistics.num_seqs_degen_per_lib = defaultdict(
            int, _updated_num_seqs_degen_per_lib
        )

    merged_statistics.to_file(out_path)


@decode.group()
def counter():
    """Decode degen command group"""
    pass


@counter.command(name="merge")
@click.argument("counter-files", type=click.Path(exists=True), required=True, nargs=-1)
@click.option(
    "-o",
    "--out-path",
    type=click.STRING,
    required=False,
    default="merged_counter.pkl",
    help="location to save merged counter file",
)
def merge_counter_file(counter_files, out_path):
    """
    Merge multiple decode counters into one

    WARNING: If using umi, you must need to merge these counter to get
    accurate umi corrected counts.

    COUNTER-FILES is a list of paths to decode counter pickles to merge.
    """
    merged_counter = pickle.load(open(counter_files[0], "rb"))
    for counter_file in counter_files[1:]:
        loaded_counter = pickle.load(open(counter_file, "rb"))
        merged_counter += loaded_counter
    pickle.dump(merged_counter, open(out_path, "wb"))


@decode.group()
def report():
    """Decode report command group"""
    pass


@report.command("merge")
@click.argument("decode-file", type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument("statistic-files", nargs=-1, type=click.Path(exists=True, dir_okay=False))
@click.option(
    "-o",
    "--out-path",
    type=click.STRING,
    required=False,
    default="merged_decoding_report.html",
    help="location to dave merged report",
)
def merge(decode_file, statistic_files, out_path):
    """
    Given decode settings and set of decode stats, merge into a single report

    Should only be used to merge stats generated from the same decode experiment
    Helpful for parallelization stuff

    DECODE is the path to a YAML file describing the decoding run.
    STATISTICS is a list of paths to decode statistics files to merge.
    """
    # decode_settings = DecodingSettings.from_file(decode_file)
    selection_info = Selection.from_yaml(decode_file)

    loaded_statistics: list[DecodeStatistics] = [
        DecodeStatistics.from_file(p) for p in statistic_files
    ]
    merged_stats = sum(loaded_statistics, DecodeStatistics())

    build_decoding_report(
        selection=selection_info,
        stats=merged_stats,
        out_path=out_path,
    )


@report.command(name="generate")
@click.argument("decode-file", type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument("statistic-file", type=click.Path(exists=True, dir_okay=False), required=True)
@click.option(
    "-o",
    "--out-path",
    type=click.STRING,
    required=False,
    default="./",
    help="path to save generated report",
)
def generate(decode_file, statistic_file, out_path):
    """
    Generate a decoding report from a decode run config and statistic file

    DECODE is the path to a YAML file describing the decoding experiment.
    STATISTIC is the path to a decode statistics file to use for the report.
    """
    # decode_settings = DecodingSettings.from_file(decode_file)
    selection_info = Selection.from_yaml(decode_file)
    loaded_statistic = DecodeStatistics.from_file(statistic_file)

    build_decoding_report(
        selection=selection_info,
        stats=loaded_statistic,
        out_path=out_path,
    )


@cli.group(name="cube")
def cube():
    """Cube command group"""
    pass


class FileType(enum.StrEnum):
    """Enum for file types"""

    CSV = "csv"
    TSV = "tsv"


@cube.command(name="from-counter")
@click.argument("counter-file", type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument("decode-file", type=click.Path(exists=True, dir_okay=False), required=True)
@click.option(
    "-o",
    "--out-path",
    type=click.Path(dir_okay=False),
    required=False,
    default=None,
    help="path to save cube file",
)
@click.option(
    "--file-format",
    "-f",
    type=click.Choice(["csv", "tsv"]),
    required=False,
    default="csv",
    help="format of the cube file; will append file type to path if not present",
)
@click.option(
    "--enumerate-smiles",
    "-e",
    is_flag=True,
    help="Enumerate the SMILES of the compounds; WARNING: can have dramatic effect on runtime",
)
@click.option(
    "--include-library-id", is_flag=True, help="Include the Library ID column in the cube file"
)
@click.option(
    "--include-bb-id", is_flag=True, help="Include the building block ID columns in the cube file"
)
@click.option(
    "--include-bb-smi",
    is_flag=True,
    help="Include the building block SMILES columns in the cube file",
)
@click.option(
    "--exclude-raw-counts", is_flag=True, help="Exclude the raw counts columns in the cube file"
)
def cube_from_counter(
    counter_file,
    decode_file,
    out_path,
    enumerate_smiles,
    file_format,
    include_library_id,
    include_bb_id,
    include_bb_smi,
    exclude_raw_counts,
):
    """
    Generate a cube file from a counter pickle

    COUNTER_FILE is the path to a counter pickle file to use for the cube.
    DECODE_FILE is the path to the decode YAML file used to generate the counters.

    Note: Enumerating the SMILES of the compounds can have a dramatic increase on runtime.
    """
    loaded_counter: DELCollectionCounter = pickle.load(open(counter_file, "rb"))
    selection = Selection.from_yaml(decode_file)

    runner_results = DecodingRunnerResults(selection=selection, degen=loaded_counter)

    if out_path is None:
        _out_path = f"{selection.selection_id}_cube.{file_format}"
    else:
        _out_path = out_path
        if not _out_path.endswith(f".{file_format}"):
            _out_path += f".{file_format}"

    runner_results.write_cube(
        out_path=_out_path,
        file_format=file_format,
        include_library_id_col=include_library_id,
        include_bb_id_cols=include_bb_id,
        include_bb_smi_cols=include_bb_smi,
        include_raw_count_col=not exclude_raw_counts,
        enumerate_smiles=enumerate_smiles,
    )


# this is broken still because the analysis code has import issues
# removing more now will fix in new branch
# @cli.command(name="analyze")
# @click.option(
#     "--config",
#     type=click.Path(exists=True, dir_okay=False),
#     required=True,
#     help="Path to the YAML config file",
# )
# def analyze(config):
#     """
#     Perform DEL analysis based on the provided YAML configuration file.
#
#     CONFIG is the path to the YAML configuration file.
#     """
#     print("Analysis started with config:", config)
#
#     def create_output_dir(output_dir):
#         if not os.path.exists(output_dir):
#             os.makedirs(output_dir)
#         return output_dir
#
#     def create_dated_output_dir(base_output_dir, name_suffix=""):
#         output_dir = create_output_dir(base_output_dir)
#         date_str = datetime.datetime.now().strftime("%Y%m%d")
#         dated_output_dir = os.path.join(output_dir, f"{date_str}_{name_suffix}")
#         if not os.path.exists(dated_output_dir):
#             os.makedirs(dated_output_dir)
#         return dated_output_dir
#
#     def load_config(config_path):
#         with open(config_path, "r") as file:
#             return yaml.safe_load(file)
#
#     def read_dict_from_file(file_path):
#         with open(file_path, "r") as f:
#             content = f.read()
#             return ast.literal_eval(content)
#
#     config = load_config(config)
#     output_dir_base = config["general"].get("output_dir", "output")
#     output_dir = create_dated_output_dir(output_dir_base, "analysis")
#
#     data_file = config["general"].get("data", "")
#     if not data_file:
#         raise ValueError("Data file must be specified in the YAML config.")
#     df = pd.read_csv(data_file)
#
#     id_col = config["general"].get("ID_col", df.columns[0])
#     indexes = ast.literal_eval(str(config.get("indexes", {})))
#     control_cols = ast.literal_eval(str(config.get("control_cols", {})))
#     raw_indexes = ast.literal_eval(str(config.get("raw_indexes", {})))
#
#     indexes_path = os.path.join(output_dir, "indexes.txt")
#     control_cols_path = os.path.join(output_dir, "control_cols.txt")
#     raw_indexes_path = os.path.join(output_dir, "raw_indexes.txt")
#
#     with open(indexes_path, "w") as file:
#         file.write(str(indexes))
#     with open(control_cols_path, "w") as file:
#         file.write(str(control_cols))
#     with open(raw_indexes_path, "w") as file:
#         file.write(str(raw_indexes))
#
#     indexes = read_dict_from_file(indexes_path)
#     control_cols = read_dict_from_file(control_cols_path)
#     raw_indexes = read_dict_from_file(raw_indexes_path)
#     if not indexes:
#         indexes = raw_indexes
#
#     cube = DELi_Cube(
#         df,
#         id_col,
#         indexes,
#         control_cols,
#         int(config["general"].get("lib_size", 0)),
#         raw_indexes,
#     )
#
#     if "flags" in config:
#         flags = config["flags"]
#         if flags.get("SD_min", False):
#             nsc_max_dict, sd_min_dict, sampling_depth_dict = cube.SD_min()
#         if flags.get("NSC_values", False):
#             cube.NSC_values()
#         if flags.get("MLE", False):
#             cube.maximum_likelihood_enrichment_ratio()
#         if flags.get("Z_score", False):
#             cube.z_score()
#         if flags.get("z_score_log_data", False):
#             cube.z_score_log_data()
#         if flags.get("disynthon_data", False):
#             disynthon_data, disynth_exp_dict = cube.disynthonize()
#             cube.data = disynthon_data
#         if flags.get("polyO", False):
#             cube.PolyO()
#         if flags.get("top_disynthons", False):
#             comparison_type = flags["top_disynthons"].get("comparison", "control")
#             exp_name = flags["top_disynthons"].get("exp_name", "None")
#             exp2_name = flags["top_disynthons"].get("exp2_name", "None")
#             control_name = flags["top_disynthons"].get("control_name", "None")
#             top_count = int(flags["top_disynthons"].get("top_count", 10))
#             comparison_metric = flags["top_disynthons"].get("comparison_metric", "avg")
#             top_disynthons_dir = create_output_dir(os.path.join(output_dir, "top_disynthons"))
#             cube.get_top_disynthons(
#                 disynthon_data=disynthon_data,
#                 exp_name1=exp_name,
#                 comparison_type=comparison_type,
#                 exp_name2=exp2_name,
#                 control_name=control_name,
#                 comparison_metric=comparison_metric,
#                 top_count=top_count,
#                 output_dir=top_disynthons_dir,
#             )
#         if flags.get("trisynthon_overlap", False):
#             trisynthon_dir = create_output_dir(os.path.join(output_dir, "trisynthon"))
#             cube.trisynthon_overlap(output_dir=trisynthon_dir)
#         if flags.get("disynthon_overlap", False):
#             disynthon_dir = create_output_dir(os.path.join(output_dir, "disynthon"))
#             cube.disynthon_overlap(
#                 output_dir=disynthon_dir,
#                 disynthon_data=disynthon_data,
#                 disynth_exp_dict=disynth_exp_dict,
#                 threshold=int(flags.get("disynthon_threshold", 20)),
#             )
#         if flags.get("normalized_data", False):
#             cube.normalize()
#         if flags.get("simple_spotfire_version", False):
#             spotfire = cube.simple_spotfire_version()
#             today_date = datetime.datetime.now().strftime("%Y%m%d")
#             spotfire.to_csv(os.path.join(output_dir, f"spotfire_{today_date}.csv"), index=False)
#         if flags.get("ml_fingerprints_to_RF_reg", False):
#             ml_fingerprints_to_RF_dir = create_output_dir(
#                 os.path.join(output_dir, "ml_fingerprints_to_RF")
#             )
#             cube.ml_fingerprints_to_RF(output_dir=ml_fingerprints_to_RF_dir)
#         if flags.get("ml_fingerprints_to_RF_clf", False):
#             ml_fingerprints_to_RF_clf_dir = create_output_dir(
#                 os.path.join(output_dir, "ml_fingerprints_to_clf")
#             )
#             cube.ml_fingerprints_to_classifier(
#                 output_dir=ml_fingerprints_to_RF_clf_dir,
#                 threshold=int(flags.get("clf_thresh", 10)),
#             )
#         if flags.get("gnn_classifier", False):
#             gnn_dir = create_output_dir(os.path.join(output_dir, "gnn"))
#             cube.gnn_classifier(
#                 output_dir=gnn_dir,
#                 threshold=int(flags.get("gnn_threshold", 10)),
#                 arch=flags.get("gnn_arch", "GAT"),
#             )
#         if "top_hits" in flags:
#             top_hits_dir = create_output_dir(os.path.join(output_dir, "top_hits"))
#             cube.top_n_compounds(
#                 int(flags["top_hits"]),
#                 flags.get("top_hits_metric", "sum"),
#                 output_dir=top_hits_dir,
#             )
#         if flags.get("report", False):
#             nsc_max_dict = nsc_max_dict if flags.get("SD_min", False) else None
#             sd_min_dict = sd_min_dict if flags.get("SD_min", False) else None
#             sampling_depth_dict = sampling_depth_dict if flags.get("SD_min", False) else None
#             generate_report(
#                 output_dir_base,
#                 indexes,
#                 control_cols,
#                 nsc_max_dict,
#                 sd_min_dict,
#                 sampling_depth_dict,
#             )
#             today_date = datetime.datetime.now().strftime("%Y%m%d")
#             cube.data.to_csv(
#                 os.path.join(output_dir_base, f"cube_data_{today_date}.csv"),
#                 index=False,
#             )
