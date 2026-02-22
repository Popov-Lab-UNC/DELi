"""Test cases for the Deli CLI commands."""

import gzip
import logging
import os
import shutil
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from deli.cli import (
    click_init_deli_config,
    click_init_deli_data_dir,
    click_set_deli_data_dir,
    click_validate_deli_data_dir,
    click_which_deli_data_dir,
    collect_decodes,
    count_compounds,
    generate_report,
    merge_stats,
    run_decode,
    summarize_decoding,
)
from deli.configure import _DeliConfig, get_deli_config, set_deli_data_dir, validate_deli_data_dir


@pytest.fixture()
def runner():
    """Create a CLI runner for testing"""
    return CliRunner()


@pytest.fixture()
def selection_file_path():
    """Path to the example selection file"""
    return Path(__file__).parent / "data" / "example_decode.yaml"


@pytest.fixture()
def fastq_file_path():
    """Path to the example fastq file"""
    return Path(__file__).parent / "data" / "example.fastq"


@pytest.fixture()
def decoded_file_path():
    """Path to the example fastq file"""
    return Path(__file__).parent / "data" / "example_decoded.tsv"


@pytest.fixture()
def collected_file_path():
    """Path to the example fastq file"""
    return Path(__file__).parent / "data" / "collected_decodes.ndjson"


@pytest.fixture()
def mock_decode_stats_file(tmpdir):
    """Create a mock decode statistics JSON file"""
    import json

    temp_home_path = Path(tmpdir)
    stats_dir = temp_home_path / "stats"
    stats_dir.mkdir(exist_ok=True)

    mock_stats = {
        "num_seqs_read": 1000,
        "num_seqs_decoded_per_lib": {
            "library_A": 800,
            "library_B": 150,
        },
        "num_failed_too_short": 10,
        "num_failed_too_long": 5,
        "num_failed_alignment": 20,
        "num_failed_library_call": 10,
        "num_failed_building_block_call": 5,
        "num_failed_ambiguous_building_block_call": 0,
        "num_failed_umi": 0,
    }

    stats_file = stats_dir / "decode_statistics.json"
    with open(stats_file, "w") as f:
        json.dump(mock_stats, f, indent=4)

    return stats_file


@pytest.fixture()
def mock_counted_decode_tsv(tmpdir):
    """Create a mock counted decode compounds TSV file"""
    temp_home_path = Path(tmpdir)

    # Create a TSV file with mock counted compounds
    counted_file = temp_home_path / "mock_counted_compounds.tsv"
    with open(counted_file, "w") as f:
        # Write header
        f.write("library_id\tbb_ids\tcount\traw_count\n")
        # Write sample data
        f.write("library_A\tBB1,BB2,BB3\t100\t150\n")
        f.write("library_A\tBB1,BB3,BB2\t50\t75\n")
        f.write("library_A\tBB2,BB1,BB3\t75\t110\n")
        f.write("library_B\tBB1,BB2,BB3\t30\t45\n")
        f.write("library_B\tBB2,BB3,BB1\t20\t30\n")

    return counted_file


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
    result = runner.invoke(click_init_deli_data_dir, ["--fix-missing", str(tmpdir_path / "broken")])
    assert result.exit_code == 0
    validate_deli_data_dir(tmpdir_path / "broken")

    # check that trying to fix a file, not directory fails
    open(tmpdir_path / "broken_file.txt", "w").close()
    result = runner.invoke(click_init_deli_data_dir, ["--fix-missing", str(tmpdir_path / "broken_file.txt")])
    assert result.exit_code == 1


def test_deli_data_which(runner):
    """Test the CLI command `deli data which`"""
    result = runner.invoke(click_which_deli_data_dir, obj={"deli_config": get_deli_config()})
    assert result.exit_code == 0
    assert "test_deli_data_dir" in result.output

    set_deli_data_dir(None)

    result = runner.invoke(click_which_deli_data_dir, obj={"deli_config": get_deli_config()})
    assert result.exit_code == 1


def test_deli_data_set(monkeypatch, tmpdir, runner):
    """Test the CLI command `deli data set`"""
    temp_home_path = Path(tmpdir)
    os.chdir(temp_home_path)

    # make a deli data directory
    runner.invoke(
        click_init_deli_data_dir,
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )

    # check for example env set output
    result = runner.invoke(
        click_set_deli_data_dir,
        ["deli_data"],
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )
    assert result.exit_code == 0
    assert "deli_data" in result.output

    # check for failure on bad data directory
    os.makedirs(temp_home_path / "broken")
    result = runner.invoke(
        click_set_deli_data_dir,
        ["broken"],
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )
    assert result.exit_code == 1

    # check for failure on file instead of dir
    open(temp_home_path / "broken_file.txt", "w").close()
    result = runner.invoke(
        click_set_deli_data_dir,
        ["broken_file.txt"],
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )
    assert result.exit_code == 1

    # make a config file
    shutil.copy2(Path(__file__).parent / "data" / ".deli", temp_home_path / ".deli")

    # test if data dir is updated in config
    result = runner.invoke(
        click_set_deli_data_dir,
        ["--update-config", "./deli_data"],
        obj={"deli_config": _DeliConfig.load_config(temp_home_path / ".deli"), "logger": logging.getLogger()},
    )
    assert result.exit_code == 0
    _config = _DeliConfig.load_config(temp_home_path / ".deli")

    assert str(_config.deli_data_dir.resolve()) == str((temp_home_path / "deli_data").resolve())


def test_decode_run(tmpdir, runner, selection_file_path, fastq_file_path):
    """Test the command `deli decode run`"""
    temp_home_path = Path(tmpdir)
    os.chdir(temp_home_path)

    shutil.copy2(fastq_file_path, temp_home_path / "example_reads.fastq")
    shutil.copy2(selection_file_path, temp_home_path / "example_decode.yaml")

    # mock in sequence file to selection
    data = yaml.load(open(temp_home_path / "example_decode.yaml"), Loader=yaml.FullLoader)
    data["sequence_files"] = [str(temp_home_path / "example_reads.fastq")]
    with open(temp_home_path / "example_decode.yaml", "w") as f:
        yaml.dump(data, f)

    result = runner.invoke(
        run_decode,
        [
            str(temp_home_path / "example_decode.yaml"),
            "-o",
            "./DecodeResults",
            "-p",
            "TEST",
            "--save-failed",
        ],
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )

    # check files exist
    assert result.exit_code == 0
    assert (temp_home_path / "DecodeResults" / "TEST_decoded.tsv").exists()
    assert (temp_home_path / "DecodeResults" / "TEST_decode_report.html").exists()
    assert (temp_home_path / "DecodeResults" / "TEST_decode_statistics.json").exists()
    assert (temp_home_path / "DecodeResults" / "TEST_failed_decoding.tsv").exists()


def test_decode_collect(tmpdir, runner, decoded_file_path):
    """Test the command `deli decode collect`"""
    temp_home_path = Path(tmpdir)
    os.chdir(temp_home_path)

    shutil.copy2(decoded_file_path, temp_home_path / "example_decoded.tsv")
    output_file = temp_home_path / "collected_results.ndjson"

    result = runner.invoke(
        collect_decodes,
        [
            str(temp_home_path / "example_decoded.tsv"),
            "--out-loc",
            str(output_file),
        ],
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )
    assert result.exit_code == 0
    assert output_file.exists()

    # rerun with compress
    output_file = temp_home_path / "collected_results.ndjson.gz"
    result = runner.invoke(
        collect_decodes,
        [
            str(temp_home_path / "example_decoded.tsv"),
            "--out-loc",
            str(output_file),
            "--compress",
        ],
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )
    assert result.exit_code == 1
    assert not output_file.exists()


def test_decode_count(tmpdir, runner, collected_file_path):
    """Test the command `deli decode count`"""
    import polars as pl

    temp_home_path = Path(tmpdir)
    os.chdir(temp_home_path)

    shutil.copy2(collected_file_path, temp_home_path / "collected_decodes.ndjson")

    # try tsv input
    output_file = temp_home_path / "count_results.tsv"

    result = runner.invoke(
        count_compounds,
        [
            str(temp_home_path / "collected_decodes.ndjson"),
            "--out-loc",
            str(output_file),
            "--output-format",
            "tsv",
        ],
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )

    assert result.exit_code == 0
    assert output_file.exists()

    # try with umi clustering
    output_file = temp_home_path / "count_results.tsv"

    result = runner.invoke(
        count_compounds,
        [
            str(temp_home_path / "collected_decodes.ndjson"),
            "--out-loc",
            str(output_file),
            "--output-format",
            "tsv",
            "-u",
            "-r",
            "-d",
        ],
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )

    assert result.exit_code == 0
    assert output_file.exists()

    # try tsv.gz input
    output_file = temp_home_path / "count_results.tsv.gz"

    result = runner.invoke(
        count_compounds,
        [
            str(temp_home_path / "collected_decodes.ndjson"),
            "--out-loc",
            str(output_file),
            "--output-format",
            "gzip",
        ],
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )

    assert result.exit_code == 0
    assert output_file.exists()
    with gzip.open(output_file, "rt", encoding="utf-8") as f_in:
        header = f_in.readline()
        assert "count" in header

    # try parquet format
    output_file = temp_home_path / "count_results.parquet"

    result = runner.invoke(
        count_compounds,
        [
            str(temp_home_path / "collected_decodes.ndjson"),
            "--out-loc",
            str(output_file),
            "--output-format",
            "parquet",
        ],
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )

    assert result.exit_code == 0
    assert output_file.exists()
    pl.read_parquet(output_file)


def test_validate_deli_data_dir(tmpdir, runner):
    """Test the validate_deli_data_dir function"""
    from test_deli_data_dir import build_invalid_deli_data_dir, build_valid_deli_data_dir

    temp_home_path = Path(tmpdir)
    os.chdir(temp_home_path)

    valid_data_dir = build_valid_deli_data_dir(temp_home_path)

    result = runner.invoke(
        click_validate_deli_data_dir,
        [str(valid_data_dir)],
    )

    assert result.exit_code == 0

    # Create an invalid deli data directory structure
    invalid_data_dir = build_invalid_deli_data_dir(temp_home_path)

    result = runner.invoke(
        click_validate_deli_data_dir,
        [str(invalid_data_dir)],
    )

    assert result.exit_code == 1
    assert "total issues found" in result.output
    assert "WARNING: Found multiple files with the same name" in result.output


def test_decode_generate_report(tmpdir, runner, mock_decode_stats_file, selection_file_path):
    """Test the command `deli decode report`"""
    temp_home_path = Path(tmpdir)
    os.chdir(temp_home_path)

    # Test generating report with a single statistics file
    output_report = temp_home_path / "test_report.html"
    result = runner.invoke(
        generate_report,
        [
            str(mock_decode_stats_file),
            str(selection_file_path),
            "--out-loc",
            str(output_report),
        ],
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )

    assert result.exit_code == 0
    assert output_report.exists()

    # Test with directory as output (should create report in directory)
    output_dir = temp_home_path / "reports"
    output_dir.mkdir(exist_ok=True)
    result = runner.invoke(
        generate_report,
        [
            str(mock_decode_stats_file),
            str(selection_file_path),
            "--out-loc",
            str(output_dir),
        ],
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )

    assert result.exit_code == 0
    assert (output_dir / "decode_report.html").exists()


def test_decode_merge_stats(tmpdir, runner, mock_decode_stats_file):
    """Test the command `deli decode merge-stats`"""
    import json

    temp_home_path = Path(tmpdir)
    os.chdir(temp_home_path)

    stats_file_1 = mock_decode_stats_file

    # Create a second copy of the stats file
    stats_file_2 = stats_file_1.parent / "TEST2_decode_statistics.json"
    shutil.copy2(stats_file_1, stats_file_2)

    # Test merging statistics files
    output_merged = temp_home_path / "merged_stats.json"
    result = runner.invoke(
        merge_stats,
        [
            str(stats_file_1),
            str(stats_file_2),
            "--out-loc",
            str(output_merged),
        ],
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )

    assert result.exit_code == 0
    assert output_merged.exists()

    # Verify the merged file is valid JSON and can be loaded
    with open(output_merged, "r") as f:
        merged_data = json.load(f)
        assert isinstance(merged_data, dict)

    # Test with default output location
    result = runner.invoke(
        merge_stats,
        [
            str(stats_file_1),
            str(stats_file_2),
        ],
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )

    assert result.exit_code == 0
    assert (temp_home_path / "merged_decode_stats.json").exists()
    import json

    with open(output_merged, "r") as f:
        merged_data = json.load(f)
        assert isinstance(merged_data, dict)

    # Test with default output location
    result = runner.invoke(
        merge_stats,
        [
            str(stats_file_1),
            str(stats_file_2),
        ],
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )

    assert result.exit_code == 0
    assert (temp_home_path / "merged_decode_stats.json").is_file()


def test_decode_summarize_decoding(tmpdir, runner, mock_decode_stats_file, mock_counted_decode_tsv):
    """Test the command `deli decode summarize`"""
    import json

    temp_home_path = Path(tmpdir)
    os.chdir(temp_home_path)

    stats_file = mock_decode_stats_file
    counted_file_tsv = mock_counted_decode_tsv

    # Test summarize with TSV input
    summary_output = temp_home_path / "decode_summary.json"
    result = runner.invoke(
        summarize_decoding,
        [
            str(counted_file_tsv),
            str(stats_file),
            "--out-loc",
            str(summary_output),
        ],
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )

    assert result.exit_code == 0
    assert summary_output.exists()

    # Verify the summary file contains expected keys
    with open(summary_output, "r") as f:
        summary_data = json.load(f)
        assert "total_seqs_read" in summary_data
        assert summary_data["total_seqs_read"] == 1000
        assert "seqs_decoded_per_lib" in summary_data
        assert summary_data["seqs_decoded_per_lib"]["library_A"] == 800
        assert summary_data["seqs_decoded_per_lib"]["library_B"] == 150
        assert "compounds_decoded_per_lib" in summary_data
        assert summary_data["compounds_decoded_per_lib"]["library_A"] == 3
        assert summary_data["compounds_decoded_per_lib"]["library_B"] == 2
        assert "molecules_decoded_per_lib" in summary_data
        assert summary_data["molecules_decoded_per_lib"]["library_A"] == 225
        assert summary_data["molecules_decoded_per_lib"]["library_B"] == 50

    # Test with default output location
    result = runner.invoke(
        summarize_decoding,
        [
            str(counted_file_tsv),
            str(stats_file),
        ],
        obj={
            "deli_config": _DeliConfig.load_config(Path(__file__).parent / "data" / ".deli"),
            "logger": logging.getLogger(),
        },
    )

    assert result.exit_code == 0
    assert (temp_home_path / "decode_summary.json").exists()
