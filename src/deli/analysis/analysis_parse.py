"""Run automated DEL analysis from a YAML configuration file."""

import argparse
import ast
import os
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from deli.analysis import analysis_report_gen as report
from deli.analysis.cube_class import DELi_Cube


def create_output_dir(output_dir: str | os.PathLike) -> str:
    """Create an output directory if it does not already exist."""
    os.makedirs(output_dir, exist_ok=True)
    return str(output_dir)


def create_dated_output_dir(base_output_dir: str | os.PathLike, name_suffix: str = "") -> str:
    """Create a dated output directory under the base output directory."""
    output_dir = create_output_dir(base_output_dir)
    date_str = datetime.now().strftime("%Y%m%d")
    dated_output_dir = os.path.join(output_dir, f"{date_str}_{name_suffix}")
    return create_output_dir(dated_output_dir)


def load_config(config_path: str | os.PathLike) -> dict:
    """Load a YAML analysis configuration file."""
    with open(config_path, "r") as file:
        return yaml.safe_load(file)


def read_dict_from_file(file_path: str | os.PathLike) -> dict:
    """Read a Python literal dict written to disk."""
    with open(file_path, "r") as file:
        content = file.read()
        return ast.literal_eval(content)


def _normalize_comparison_block(block: Any, defaults: dict) -> list[dict]:
    """
    Normalize config that can be either:
    - a single comparison dict, or
    - a dict with `comparisons: [...]`.
    """
    if not block:
        return []
    if isinstance(block, dict) and isinstance(block.get("comparisons"), list):
        out = []
        for item in block["comparisons"]:
            if not isinstance(item, dict):
                continue
            merged = defaults.copy()
            merged.update(item)
            out.append(merged)
        return out
    if isinstance(block, dict):
        merged = defaults.copy()
        merged.update(block)
        return [merged]
    return []


def run_analysis(config_path: str | os.PathLike) -> Path:
    """
    Perform DEL analysis based on the provided YAML configuration file.

    Parameters
    ----------
    config_path : str | os.PathLike
        Path to the YAML configuration file.

    Returns
    -------
    Path
        The dated analysis output directory.
    """
    config_path = Path(config_path).resolve()
    config = load_config(config_path)
    config_dir = config_path.parent

    output_dir_base = config["general"].get("output_dir", "output")
    if not Path(output_dir_base).is_absolute():
        output_dir_base = str((config_dir / output_dir_base).resolve())
    output_dir = create_dated_output_dir(output_dir_base, "analysis")

    data_file = config["general"].get("data", "")
    if not data_file:
        raise ValueError("Data file must be specified in the YAML config.")
    data_path = Path(data_file)
    if not data_path.is_absolute():
        data_path = (config_dir / data_path).resolve()
    df = pd.read_csv(data_path)

    if "ID_col" in config["general"]:
        id_col = config["general"]["ID_col"]
    else:
        id_col = df.columns[0]

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

    disynthon_data = None
    disynth_exp_dict = None
    nsc_max_dict = None
    sd_min_dict = None
    sampling_depth_dict = None
    top_disynthon_specs: list[dict] = []
    top_delta_specs: list[dict] = []

    print(cube.data.head())
    if "flags" not in config:
        return Path(output_dir)

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
        print(disynth_exp_dict)
    if flags.get("polyO", False):
        cube.PolyO()
    if flags.get("top_disynthons", False):
        top_disynthon_specs = _normalize_comparison_block(
            flags.get("top_disynthons", {}),
            {
                "comparison": "control",
                "exp_name": "None",
                "exp2_name": "None",
                "control_name": "None",
                "top_count": 10,
                "comparison_metric": "avg",
            },
        )
        top_disynthons_dir = create_output_dir(os.path.join(output_dir, "top_disynthons"))
        for spec in top_disynthon_specs:
            cube.get_top_disynthons(
                disynthon_data=disynthon_data,
                exp_name1=spec["exp_name"],
                comparison_type=spec["comparison"],
                exp_name2=spec.get("exp2_name"),
                control_name=spec.get("control_name"),
                comparison_metric=spec.get("comparison_metric", "avg"),
                top_count=int(spec.get("top_count", 10)),
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
            threshold=flags.get("disynthon_threshold", 20),
        )
    if flags.get("normalized_data", False):
        cube.normalize()
    if flags.get("simple_spotfire_version", False):
        spotfire = cube.simple_spotfire_version()
        today_date = datetime.now().strftime("%Y%m%d")
        spotfire.to_csv(os.path.join(output_dir, f"spotfire_{today_date}.csv"), index=False)
    if flags.get("ml_fingerprints_to_RF_reg", False):
        ml_fingerprints_to_rf_dir = create_output_dir(os.path.join(output_dir, "ml_fingerprints_to_RF"))
        cube.ml_fingerprints_to_RF(output_dir=ml_fingerprints_to_rf_dir)
    if flags.get("ml_fingerprints_to_RF_clf", False):
        ml_fingerprints_to_rf_clf_dir = create_output_dir(
            os.path.join(output_dir, "ml_fingerprints_to_clf")
        )
        cube.ml_fingerprints_to_classifier(
            output_dir=ml_fingerprints_to_rf_clf_dir,
            threshold=int(flags.get("clf_thresh", 10)),
        )
    if flags.get("gnn_classifier", False):
        gnn_dir = create_output_dir(os.path.join(output_dir, "gnn"))
        try:
            cube.gnn_classifier(
                output_dir=gnn_dir,
                threshold=int(flags.get("gnn_threshold", 10)),
                arch=flags.get("gnn_arch", "GAT"),
            )
        except ImportError as exc:
            raise ImportError(
                "GNN analysis requires optional ML dependencies. Install with: pip install 'deli-chem[ml]'"
            ) from exc
    if "top_hits" in flags:
        top_hits_dir = create_output_dir(os.path.join(output_dir, "top_hits"))
        cube.top_n_compounds(
            int(flags["top_hits"]),
            flags.get("top_hits_metric", "sum"),
            output_dir=top_hits_dir,
        )
    if flags.get("top_delta_compounds", False):
        top_delta_specs = _normalize_comparison_block(
            flags.get("top_delta_compounds", {}),
            {
                "exp_name": "None",
                "control_name": "None",
                "metric": "avg",
                "top_count": 20,
            },
        )
        top_delta_dir = create_output_dir(os.path.join(output_dir, "top_deltas"))
        cube.top_delta_compounds(top_delta_specs, output_dir=top_delta_dir)
    if flags.get("monosynthon_chemical_space", False):
        mono_cfg = flags.get("monosynthon_chemical_space")
        if isinstance(mono_cfg, dict):
            mono_experiments = mono_cfg.get("experiments")
            cube.monosynthon_chemical_space(
                output_dir=output_dir,
                experiments=mono_experiments,
            )
        else:
            cube.monosynthon_chemical_space(output_dir=output_dir)
    if flags.get("report", False):
        report.generate_report(
            output_dir_base,
            indexes,
            control_cols,
            nsc_max_dict,
            sd_min_dict,
            sampling_depth_dict,
            top_disynthon_specs=top_disynthon_specs,
            top_delta_specs=top_delta_specs,
            disynthon_thresholds=flags.get("disynthon_threshold"),
        )
        print("Report generation completed!")
        today_date = datetime.now().strftime("%Y%m%d")
        cube.data.to_csv(
            os.path.join(output_dir_base, f"cube_data_{today_date}.csv"),
            index=False,
        )
        print(f"Cube data saved to {os.path.join(output_dir, 'cube_data.csv')}")

    return Path(output_dir)


def parse_args():
    """Parse command-line arguments for standalone analysis runs."""
    parser = argparse.ArgumentParser(
        description="Parse the necessary arguments for DEL analysis of one library"
    )
    parser.add_argument("--config", type=str, help="Path to the YAML config file")
    return parser.parse_args()


def main():
    """Entry point for running analysis from the command line."""
    args = parse_args()
    if not args.config:
        raise ValueError("Config file path must be provided with --config")

    output_dir = run_analysis(args.config)
    print(f"Analysis completed. Results saved to: {output_dir}")


if __name__ == "__main__":
    main()
