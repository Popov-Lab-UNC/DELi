"""test cases for deli data directory configuration"""

import os
from pathlib import Path

import pytest

from deli.configure import (
    DELI_DATA_EXTENSIONS,
    DELI_DATA_SUB_DIRS,
    DeliDataDirError,
    get_deli_config,
    resolve_deli_data_name,
    set_deli_data_dir,
    validate_deli_data_dir,
)


def build_valid_deli_data_dir(tmp_path: Path) -> Path:
    """Build a valid deli data directory structure for testing"""
    valid_data_dir = tmp_path / "valid_deli_data"
    valid_data_dir.mkdir()
    for sub_dir, ext in zip(DELI_DATA_SUB_DIRS, DELI_DATA_EXTENSIONS, strict=True):
        (valid_data_dir / sub_dir).mkdir()
        with open(valid_data_dir / sub_dir / f"placeholder1.{ext}", "w") as f:
            f.write("This is a placeholder file.")
        with open(valid_data_dir / sub_dir / "README.md", "w") as f:
            f.write("This is a placeholder file with the wrong extension.")
        (valid_data_dir / sub_dir / "sub_sub_dir1").mkdir()
        with open(valid_data_dir / sub_dir / "sub_sub_dir1" / f"placeholder2.{ext}", "w") as f:
            f.write("This is a placeholder file.")
        (valid_data_dir / sub_dir / "sub_sub_dir1" / "sub_sub_sub_dir").mkdir()
        with open(valid_data_dir / sub_dir / "sub_sub_dir1" / "sub_sub_sub_dir" / f"placeholder3.{ext}", "w") as f:
            f.write("This is a placeholder file.")
        (valid_data_dir / sub_dir / "sub_sub_dir2").mkdir()
        with open(valid_data_dir / sub_dir / "sub_sub_dir1" / f"placeholder4.{ext}", "w") as f:
            f.write("This is a placeholder file.")

    return valid_data_dir


def build_invalid_deli_data_dir(temp_home_path: Path) -> Path:
    """Build an invalid deli data directory structure for testing"""
    invalid_data_dir = temp_home_path / "invalid_deli_data"
    invalid_data_dir.mkdir()
    for sub_dir, ext in zip(DELI_DATA_SUB_DIRS[1:], DELI_DATA_EXTENSIONS[1:], strict=True):
        (invalid_data_dir / sub_dir).mkdir()
        with open(invalid_data_dir / sub_dir / f"placeholder1.{ext}", "w") as f:
            f.write("This is a placeholder file.")
        with open(invalid_data_dir / sub_dir / "README.md", "w") as f:
            f.write("This is a placeholder file with the wrong extension.")
        (invalid_data_dir / sub_dir / "sub_sub_dir1").mkdir()
        # this file is a duplicate to cause failure
        with open(invalid_data_dir / sub_dir / "sub_sub_dir1" / f"placeholder1.{ext}", "w") as f:
            f.write("This is a placeholder file.")
        (invalid_data_dir / sub_dir / "sub_sub_dir1" / "sub_sub_sub_dir").mkdir()
        with open(invalid_data_dir / sub_dir / "sub_sub_dir1" / "sub_sub_sub_dir" / f"placeholder3.{ext}", "w") as f:
            f.write("This is a placeholder file.")
        (invalid_data_dir / sub_dir / "sub_sub_dir2").mkdir()
        with open(invalid_data_dir / sub_dir / "sub_sub_dir1" / f"placeholder4.{ext}", "w") as f:
            f.write("This is a placeholder file.")
    return invalid_data_dir


def test_set_deli_data_dir(tmp_path: Path):
    """Test setting the deli data directory."""
    test_data_dir = tmp_path / "deli_data"
    test_data_dir.mkdir()

    valid_data_dir = build_valid_deli_data_dir(test_data_dir)
    set_deli_data_dir(valid_data_dir)

    assert get_deli_config().deli_data_dir == valid_data_dir
    assert os.environ["DELI_DATA_DIR"] == str(valid_data_dir)


def test_validation_warning_deli_data_dir(tmp_path: Path):
    """Test setting the deli data directory."""
    test_data_dir = tmp_path / "deli_data"
    test_data_dir.mkdir()

    invalid_dir = build_invalid_deli_data_dir(test_data_dir)

    with pytest.warns(Warning) as record:
        validate_deli_data_dir(invalid_dir)

    assert any("is missing sub-directories" in str(warning.message) for warning in record)
    assert any("duplicate names" in str(warning.message) for warning in record)


def test_resolve_deli_data_name(tmp_path: Path):
    """Test that setting a deli data directory with name conflicts raises an error."""
    test_data_dir = tmp_path / "deli_data"
    test_data_dir.mkdir()

    invalid_dir = build_invalid_deli_data_dir(test_data_dir)
    with pytest.warns(Warning):
        get_deli_config().deli_data_dir = invalid_dir

    @resolve_deli_data_name(DELI_DATA_SUB_DIRS[1], DELI_DATA_EXTENSIONS[1], target_param="path")
    def dummy_load(path):
        pass

    with pytest.raises(DeliDataDirError, match="multiple files named"):
        dummy_load("placeholder1")

    dummy_load("placeholder_not_exist")  # doesn't validate existence of file or name

    dummy_load("placeholder3")  # no duplicates should work
