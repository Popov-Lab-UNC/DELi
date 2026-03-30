"""Tests for the deli.cube package."""

import polars as pl
import pytest

from deli.cube import Cube, CubeColumnError


def _example_polars_df() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "lib_col": ["L1", "L1"],
            "compound_col": ["L1-A-B", "L1-C-D"],
            "bb_a": ["A", "C"],
            "bb_b": ["B", "D"],
            "count_1": [10, 3],
            "count_2": [4, 1],
            "enrich_1": [2.1, 0.2],
            "syn_a": ["A", "C"],
            "syn_b": ["B", "D"],
        }
    )


def test_cube_initializes_with_valid_mappings():
    """Cube accepts valid role mappings and keeps the data."""
    data = _example_polars_df()

    cube = Cube.from_dataframe(
        data,
        library_id="lib_col",
        bbid=["bb_a", "bb_b"],
        count=["count_1", "count_2"],
        enrichment=["enrich_1"],
        synthon={1: "syn_a", 2: "syn_b"},
    )

    assert cube.data.shape == data.shape
    assert cube.primary_id_column == "lib_col"
    assert cube.synthon_types == (1, 2)
    assert cube.synthon_column(1) == "syn_a"


def test_cube_from_pandas_converts_to_polars():
    """Cube can be loaded from pandas and converted to polars."""
    pd = pytest.importorskip("pandas")
    data = pd.DataFrame({"compound_col": ["CMP-1", "CMP-2"], "count_1": [1, 2]})

    cube = Cube.from_pandas(
        data,
        compound_id="compound_col",
        count=["count_1"],
    )

    assert isinstance(cube.data, pl.DataFrame)
    assert cube.primary_id_column == "compound_col"


def test_cube_rejects_both_primary_ids():
    """Exactly one primary ID role is allowed."""
    data = _example_polars_df()

    with pytest.raises(CubeColumnError, match="Exactly one of 'library_id' or 'compound_id' must be configured"):
        Cube.from_dataframe(
            data,
            library_id="lib_col",
            compound_id="compound_col",
        )


def test_cube_rejects_missing_primary_ids():
    """Exactly one primary ID role is required."""
    data = _example_polars_df()

    with pytest.raises(CubeColumnError, match="Exactly one of 'library_id' or 'compound_id' must be configured"):
        Cube.from_dataframe(data)


def test_cube_rejects_missing_columns():
    """Mapped columns must exist in the dataframe."""
    data = _example_polars_df()

    with pytest.raises(CubeColumnError, match="missing columns"):
        Cube.from_dataframe(
            data,
            library_id="lib_col",
            count=["count_not_present"],
        )


def test_cube_rejects_duplicate_column_assignment():
    """A single column cannot be mapped to multiple roles."""
    data = _example_polars_df()

    with pytest.raises(CubeColumnError, match="assigned to multiple roles"):
        Cube.from_dataframe(
            data,
            library_id="lib_col",
            bbid=["bb_a"],
            synthon={1: "bb_a"},
        )


def test_cube_rejects_non_contiguous_synthon_types():
    """Synthon types must be contiguous integers 1..N."""
    data = _example_polars_df()

    with pytest.raises(CubeColumnError, match="contiguous integers"):
        Cube.from_dataframe(
            data,
            library_id="lib_col",
            synthon={1: "syn_a", 3: "syn_b"},
        )


def test_cube_accepts_synthon_iterable_input():
    """Synthon list input is interpreted as types 1..N."""
    data = _example_polars_df()

    cube = Cube.from_dataframe(
        data,
        library_id="lib_col",
        synthon=["syn_a", "syn_b"],
    )

    assert cube.synthon_types == (1, 2)
    assert cube.synthon_column(2) == "syn_b"
