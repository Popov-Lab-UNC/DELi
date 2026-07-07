"""configure testing environment for deli"""

import os
from pathlib import Path

import pandas as pd
import pytest

from deli.analysis.cube_class import DELi_Cube
from deli.configure import init_deli_config, set_deli_data_dir

ANALYSIS_INDEXES = {"sel": ["corrected_idx1", "corrected_idx2", "corrected_idx3"]}
ANALYSIS_RAW_INDEXES = {"sel": ["raw_idx1", "raw_idx2", "raw_idx3"]}
ANALYSIS_CONTROL_COLS = {"sel": ["control_sel"], "ntc": ["control_sel"]}
ML_SMILES = ["CCO", "CCC", "CC(C)O", "CCCC", "c1ccccc1", "CCN", "CCOC", "CC(C)C"]


@pytest.fixture
def minimal_cube_df() -> pd.DataFrame:
    return pd.read_csv(Path(__file__).parent / "data" / "analysis" / "minimal_cube.csv")


@pytest.fixture
def ml_cube_df() -> pd.DataFrame:
    """Synthetic cube with SMILES and enough rows for ML/GNN sampling."""
    rows = []
    for i in range(350):
        row_idx = i % 40
        rows.append(
            {
                "DEL_ID": f"A{row_idx:03d}-B{row_idx:03d}-C{row_idx:03d}",
                "corrected_idx1": max(1, 25 - i // 12) + (i % 4),
                "corrected_idx2": max(1, 15 - i // 18) + (i % 3),
                "corrected_idx3": max(1, 8 - i // 24) + (i % 2),
                "raw_idx1": max(1, 25 - i // 12) + (i % 4),
                "raw_idx2": max(1, 15 - i // 18) + (i % 3),
                "raw_idx3": max(1, 8 - i // 24) + (i % 2),
                "control_sel": max(1, 4 - i // 80),
                "SMILES": ML_SMILES[i % len(ML_SMILES)],
                "SMILES_A": "CCO",
                "SMILES_B": "CCC",
                "SMILES_C": "CC(C)O",
            }
        )
    return pd.DataFrame(rows)


@pytest.fixture
def ml_cube(ml_cube_df: pd.DataFrame) -> DELi_Cube:
    return DELi_Cube(
        ml_cube_df,
        "DEL_ID",
        ANALYSIS_INDEXES,
        control_cols=ANALYSIS_CONTROL_COLS,
        lib_size=1000,
        raw_indexes=ANALYSIS_RAW_INDEXES,
    )


@pytest.fixture(autouse=True, scope="function")
def set_up_test_env():
    """Set up the test environment by setting the necessary environment variables."""
    os.environ["DELI_CONFIG"] = os.path.join(os.path.dirname(__file__), "data", ".deli")
    os.environ["DELI_DATA_DIR"] = os.path.join(os.path.dirname(__file__), "data", "test_deli_data_dir")

    init_deli_config(Path(os.path.join(os.path.dirname(__file__), "data", ".deli")), fail_on_exist=False)
    set_deli_data_dir(os.path.join(os.path.dirname(__file__), "data", "test_deli_data_dir"))
