"""Tests for DEL analysis helpers, DELi_Cube methods, and the `deli analyze` workflow."""

from __future__ import annotations

import shutil
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pandas as pd
import pytest
import yaml
from click.testing import CliRunner

from deli.analysis import analysis_report_gen as report
from deli.analysis.analysis_parse import (
    _normalize_comparison_block,
    load_config,
    read_dict_from_file,
    run_analysis,
)
from deli.analysis.cube_class import DELi_Cube
from deli.analysis.poly_o import PolyO
from deli.cli import analyze

ANALYSIS_FIXTURE_DIR = Path(__file__).parent / "data" / "analysis"
INDEXES = {"sel": ["corrected_idx1", "corrected_idx2", "corrected_idx3"]}
RAW_INDEXES = {"sel": ["raw_idx1", "raw_idx2", "raw_idx3"]}


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def minimal_cube() -> DELi_Cube:
    df = pd.read_csv(ANALYSIS_FIXTURE_DIR / "minimal_cube.csv")
    return DELi_Cube(df, "DEL_ID", INDEXES, control_cols={}, lib_size=1000, raw_indexes=RAW_INDEXES)


@pytest.fixture
def analysis_workspace(tmp_path: Path) -> Path:
    """Copy bundled analysis fixtures into an isolated working directory."""
    shutil.copy(ANALYSIS_FIXTURE_DIR / "minimal_cube.csv", tmp_path / "minimal_cube.csv")
    shutil.copy(ANALYSIS_FIXTURE_DIR / "minimal_config.yaml", tmp_path / "minimal_config.yaml")
    return tmp_path


@pytest.mark.parametrize(
    ("block", "expected_count"),
    [
        (None, 0),
        (
            {
                "comparison": "none",
                "exp_name": "sel",
                "top_count": 5,
            },
            1,
        ),
        (
            {
                "comparisons": [
                    {"comparison": "none", "exp_name": "sel", "top_count": 3},
                    {"comparison": "control", "exp_name": "sel", "control_name": "ctrl"},
                ]
            },
            2,
        ),
    ],
)
def test_normalize_comparison_block(block, expected_count):
    defaults = {"comparison": "control", "exp_name": "None", "top_count": 10}
    specs = _normalize_comparison_block(block, defaults)
    assert len(specs) == expected_count
    if specs:
        assert specs[0]["comparison"] in {"none", "control"}


def test_load_config_round_trip(analysis_workspace: Path):
    config = load_config(analysis_workspace / "minimal_config.yaml")
    assert config["general"]["ID_col"] == "DEL_ID"
    assert config["flags"]["disynthon_data"] is True


def test_read_dict_from_file_round_trip(tmp_path: Path):
    payload = {"sel": ["corrected_idx1", "corrected_idx2"]}
    path = tmp_path / "indexes.txt"
    path.write_text(str(payload))
    assert read_dict_from_file(path) == payload


def test_deli_cube_splits_del_id_into_synthon_columns(minimal_cube_df: pd.DataFrame):
    cube = DELi_Cube(minimal_cube_df, "DEL_ID", INDEXES, lib_size=1000)
    assert list(cube.data.columns[:4]) == ["DEL_ID", "ID_A", "ID_B", "ID_C"]
    assert cube.data.loc[0, "ID_A"] == "A001"
    assert cube.data.loc[0, "ID_B"] == "B001"
    assert cube.data.loc[0, "ID_C"] == "C001"


def test_deli_cube_sd_min_returns_qc_metrics(minimal_cube: DELi_Cube):
    nsc_max, sd_min, sampling_depth = minimal_cube.SD_min()
    assert nsc_max["sel_NSC_max"] > 0
    assert sd_min["sel_SD_min"] > 0
    assert sampling_depth["sel_sampling_depth"] > 0


def test_deli_cube_nsc_values_adds_enrichment_columns(minimal_cube: DELi_Cube):
    minimal_cube.NSC_values()
    assert "sel_NSC" in minimal_cube.data.columns
    assert "sel_sum" in minimal_cube.data.columns
    assert minimal_cube.data["sel_NSC"].notna().all()


def test_deli_cube_disynthonize_adds_pair_columns_and_counts(minimal_cube: DELi_Cube):
    disynthon_data, exp_dict = minimal_cube.disynthonize()
    assert {"AB", "AC", "BC"}.issubset(disynthon_data.columns)
    assert "AB_sel" in exp_dict
    assert any(col.startswith("count_AB_sel_") for col in disynthon_data.columns)


def test_deli_cube_top_n_compounds_writes_structure_plot(tmp_path: Path, minimal_cube: DELi_Cube):
    minimal_cube.NSC_values()
    minimal_cube.data["SMILES"] = "CCO"
    minimal_cube.top_n_compounds(n=2, metric="sum", output_dir=str(tmp_path))
    assert list(tmp_path.glob("sel_top_2_compounds.svg"))


def test_deli_cube_z_score_adds_normalized_columns(ml_cube: DELi_Cube):
    ml_cube.z_score()
    assert "sel_z_score" in ml_cube.data.columns
    assert "sel_norm_z_score" in ml_cube.data.columns
    assert ml_cube.data["sel_z_score"].notna().all()


def test_deli_cube_z_score_requires_lib_size(minimal_cube_df: pd.DataFrame):
    cube = DELi_Cube(minimal_cube_df, "DEL_ID", INDEXES, lib_size=0)
    with pytest.raises(ValueError, match="Positive library size"):
        cube.z_score()


def test_deli_cube_normalize_subtracts_control(ml_cube: DELi_Cube):
    before = ml_cube.data["corrected_idx1"].copy()
    ml_cube.normalize()
    assert not ml_cube.data["corrected_idx1"].equals(before)
    assert "sel_avg" in ml_cube.data.columns
    assert ml_cube._exclude_control_cols_from_spotfire is True


def test_deli_cube_mle_adds_ratio_column(ml_cube: DELi_Cube):
    ml_cube.maximum_likelihood_enrichment_ratio()
    assert "sel_MLE" in ml_cube.data.columns
    assert ml_cube.data["sel_MLE"].notna().all()


def test_deli_cube_nsc_enrichment_intervals(ml_cube: DELi_Cube):
    nsc_df = ml_cube.NSC_values()
    ml_cube.NSC_enrichment_intervals(nsc_df)
    assert "sel_NSC+" in ml_cube.data.columns
    assert "sel_NSC-" in ml_cube.data.columns


def test_deli_cube_polyo_adds_score_columns(ml_cube: DELi_Cube):
    ml_cube.disynthonize()
    ml_cube.PolyO()
    assert any(col.endswith("_PolyO_score") for col in ml_cube.data.columns)


def test_polyo_pipeline_on_disynthonized_data(ml_cube: DELi_Cube):
    disynthon_data, _ = ml_cube.disynthonize()
    polyo = PolyO(disynthon_data, INDEXES, lib_size=1000, feature_mode="disynthon")
    polyo.calculate_polyOraw()
    ccpd = polyo.find_Ccpd()
    cread = polyo.calculate_Cread()
    polyo.poly_o_base()
    result = polyo.calculate_polyO_score()
    assert ccpd
    assert cread
    assert any("PolyO_score" in col for col in result.columns)


def test_deli_cube_simple_spotfire_version_keeps_summary_columns(ml_cube: DELi_Cube):
    ml_cube.NSC_values()
    spotfire = ml_cube.simple_spotfire_version()
    assert "SMILES" in spotfire.columns
    assert "sel_sum" in spotfire.columns


def test_deli_cube_ml_fingerprints_to_rf_writes_plot(ml_cube: DELi_Cube, tmp_path: Path):
    ml_cube.ml_fingerprints_to_RF(output_dir=str(tmp_path))
    assert list(tmp_path.glob("sel_RF_regression.svg"))


def test_deli_cube_ml_fingerprints_classifier_writes_plot(ml_cube: DELi_Cube, tmp_path: Path):
    ml_cube.ml_fingerprints_to_classifier(threshold=2, output_dir=str(tmp_path))
    assert list(tmp_path.glob("sel_classifier.svg"))


def test_deli_cube_top_delta_compounds_writes_plot(ml_cube: DELi_Cube, tmp_path: Path):
    specs = [
        {
            "exp_name": "sel",
            "control_name": "ntc",
            "metric": "avg",
            "top_count": 5,
        }
    ]
    ml_cube.top_delta_compounds(specs, output_dir=str(tmp_path))
    assert list(tmp_path.glob("sel_vs_ntc_top_5_avg_deltas.svg"))


def test_deli_cube_monosynthon_chemical_space_writes_outputs(ml_cube: DELi_Cube, tmp_path: Path):
    clusters_dir = tmp_path / "clusters"
    ml_cube.monosynthon_chemical_space(output_dir=str(tmp_path), experiments=["sel"])
    assert clusters_dir.is_dir()
    assert list(clusters_dir.glob("monosynthon_A_sel_chemical_space.html"))
    assert list(clusters_dir.glob("monosynthon_A_sel_chemical_space.svg"))


def test_deli_cube_z_score_log_data_with_controls(ml_cube: DELi_Cube):
    ml_cube.z_score_log_data()
    assert "sel_z_score_log" in ml_cube.data.columns
    assert ml_cube.data["sel_z_score_log"].notna().all()


def test_get_experiment_info():
    info = report.get_experiment_info(INDEXES, {})
    assert info == [
        {
            "name": "sel",
            "index": ["corrected_idx1", "corrected_idx2", "corrected_idx3"],
            "control_columns": [],
        }
    ]


def test_run_analysis_missing_data_raises(tmp_path: Path):
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "general": {"data": "", "output_dir": "output"},
                "indexes": INDEXES,
                "control_cols": {},
                "raw_indexes": RAW_INDEXES,
                "flags": {},
            }
        )
    )
    with pytest.raises(ValueError, match="Data file must be specified"):
        run_analysis(config_path)


def test_run_analysis_minimal_config(analysis_workspace: Path):
    output_dir = run_analysis(analysis_workspace / "minimal_config.yaml")

    assert output_dir.is_dir()
    assert (output_dir / "indexes.txt").exists()
    assert (output_dir / "top_disynthons").is_dir()
    assert list((output_dir / "top_disynthons").glob("*.svg"))
    assert list((output_dir / "trisynthon").glob("*.svg"))
    assert list((output_dir / "disynthon").glob("*.svg"))

    today = datetime.now().strftime("%Y%m%d")
    report_path = analysis_workspace / "output" / f"{today}_analysis_report.html"
    assert report_path.exists()
    assert "experiment_info" not in report_path.read_text(encoding="utf-8")


def test_deli_analyze_cli_minimal(runner: CliRunner, analysis_workspace: Path):
    result = runner.invoke(
        analyze,
        ["--config", str(analysis_workspace / "minimal_config.yaml")],
    )
    assert result.exit_code == 0, result.output
    assert "Analysis completed" in result.output
