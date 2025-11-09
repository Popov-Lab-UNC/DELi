"""Test cases for the Deli CLI commands."""

import os
from pathlib import Path

import pytest
from click.testing import CliRunner

from deli.cli import (
    click_init_deli_config,
    click_init_deli_data_dir,
    click_set_deli_data_dir,
    click_which_deli_data_dir,
)
from deli.configure import _DeliConfig, set_deli_data_dir, validate_deli_data_dir


@pytest.fixture()
def runner():
    """Create a CLI runner for testing"""
    return CliRunner()


@pytest.mark.functional
def test_deli_config_init(monkeypatch, tmpdir, runner):
    """Testing the command `deli config init`"""
    # path the home directory to the temporary directory
    temp_home_path = Path(tmpdir)
    monkeypatch.setattr(Path, "home", lambda: temp_home_path)

    result = runner.invoke(click_init_deli_config)

    assert result.exit_code == 0

    file_path = temp_home_path / ".deli"
    assert file_path.exists()  # check file exists
    _DeliConfig.load_config(file_path)  # check config file is valid

    # check for failure when trying init again
    result = runner.invoke(click_init_deli_config)
    assert result.exit_code == 1

    # check that existing config can be overwritten
    result = runner.invoke(click_init_deli_config, ["--overwrite"])
    assert result.exit_code == 0
    assert file_path.exists()
    _DeliConfig.load_config(file_path)

    # check custom path for config file
    custom_config_path = temp_home_path / ".custom_deli_config"
    result = runner.invoke(click_init_deli_config, [str(custom_config_path)])
    assert result.exit_code == 0
    assert custom_config_path.exists()
    _DeliConfig.load_config(custom_config_path)


@pytest.mark.functional
def test_deli_data_init(tmpdir, runner):
    """Testing the command `deli data init`"""
    os.chdir(tmpdir)
    tmpdir_path = Path(tmpdir)

    result = runner.invoke(click_init_deli_data_dir)

    assert result.exit_code == 0

    tmpdir_path = tmpdir_path / "deli_data"
    assert tmpdir_path.exists()  # check file exists
    assert tmpdir_path.is_dir()
    validate_deli_data_dir(tmpdir_path)

    # check for failure when trying init again
    result = runner.invoke(click_init_deli_data_dir)
    assert result.exit_code == 1

    # check ability to overwrite when trying init again
    result = runner.invoke(click_init_deli_data_dir, ["--overwrite"])
    assert result.exit_code == 0
    assert tmpdir_path.exists()
    validate_deli_data_dir(tmpdir_path)

    # check that broken data dir can be fixes
    os.makedirs(tmpdir_path / "broken")
    result = runner.invoke(
        click_init_deli_data_dir, ["--fix-missing", str(tmpdir_path / "broken")]
    )
    assert result.exit_code == 0
    validate_deli_data_dir(tmpdir_path / "broken")

    # check that trying to fix a file, not directory fails
    open(tmpdir_path / "broken_file.txt", "w").close()
    result = runner.invoke(
        click_init_deli_data_dir, ["--fix-missing", str(tmpdir_path / "broken_file.txt")]
    )
    assert result.exit_code == 1


@pytest.mark.functional
def test_deli_data_which(runner):
    """Test the CLI command `deli data which`"""
    result = runner.invoke(click_which_deli_data_dir)
    assert result.exit_code == 0
    assert "test_deli_data_dir" in result.output

    set_deli_data_dir(None)

    result = runner.invoke(click_which_deli_data_dir)
    assert result.exit_code == 1


@pytest.mark.functional
def test_deli_data_set(monkeypatch, tmpdir, runner):
    """Test the CLI command `deli data set`"""
    temp_home_path = Path(tmpdir)
    os.chdir(temp_home_path)

    # make a deli data directory
    runner.invoke(click_init_deli_data_dir)

    # check for example env set output
    result = runner.invoke(click_set_deli_data_dir, ["deli_data"])
    assert result.exit_code == 0
    assert "deli_data" in result.output

    # check for failure on bad data directory
    os.makedirs(temp_home_path / "broken")
    result = runner.invoke(click_set_deli_data_dir, ["broken"])
    assert result.exit_code == 1

    # check for failure on file instead of dir
    open(temp_home_path / "broken_file.txt", "w").close()
    result = runner.invoke(click_set_deli_data_dir, ["broken_file.txt"])
    assert result.exit_code == 1

    # make a config file
    monkeypatch.setattr(Path, "home", lambda: temp_home_path)
    runner.invoke(click_init_deli_config)

    # test if data dir is updated in config
    result = runner.invoke(click_set_deli_data_dir, ["--update-config", "deli_data"])
    assert result.exit_code == 0
    _config = _DeliConfig.load_config(temp_home_path / ".deli")

    assert str(_config.deli_data_dir.resolve()) == str((temp_home_path / "deli_data").resolve())
