"""tests for decoding classes/functions"""

import os

import pytest

from deli.configure import set_deli_data_dir
from deli.decode import DecodingRunner
from deli.decode.runner import DecodingRunParsingError


set_deli_data_dir(os.path.join(os.path.dirname(__file__), "test_data", "test_deli_data_dir"))

DECODE_FILE = os.path.join(os.path.dirname(__file__), "test_data", "example_decode.yaml")
NO_SEQ_DECODE_FILE = os.path.join(
    os.path.dirname(__file__), "test_data", "example_decode_no_seqs.yaml"
)
FASTQ_FILE = os.path.join(os.path.dirname(__file__), "test_data", "example.fastq")


@pytest.mark.unit
def test_loading_decode_runner():
    """Test loading a DecodingRunner from a file"""
    runner = DecodingRunner.from_file(DECODE_FILE)
    assert isinstance(runner, DecodingRunner)

    runner = DecodingRunner.from_file(DECODE_FILE, [FASTQ_FILE], ignore_decode_seqs=True)
    assert isinstance(runner, DecodingRunner)

    runner = DecodingRunner.from_file(NO_SEQ_DECODE_FILE, [FASTQ_FILE], ignore_decode_seqs=False)
    assert isinstance(runner, DecodingRunner)

    with pytest.raises(DecodingRunParsingError, match="`fastq_files` cannot be provided when"):
        DecodingRunner.from_file(DECODE_FILE, [FASTQ_FILE])

    with pytest.raises(
        DecodingRunParsingError, match="``ignore_decode_seqs`` is True, but no `fastq_files`"
    ):
        DecodingRunner.from_file(DECODE_FILE, ignore_decode_seqs=True)

    with pytest.raises(DecodingRunParsingError, match="found neither"):
        DecodingRunner.from_file(NO_SEQ_DECODE_FILE)


@pytest.fixture
def runner() -> DecodingRunner:
    """A DecodingRunner fixture for testing"""
    return DecodingRunner.from_file(DECODE_FILE, disable_logging=True)


@pytest.mark.unit
def test_decoding_runner_to_file(runner, tmpdir):
    """Test saving a DecodingRunner to a file"""
    runner.to_file(tmpdir.join("decode_test.yaml"))
    runner_2 = DecodingRunner.from_file(tmpdir.join("decode_test.yaml"))

    assert isinstance(runner_2, DecodingRunner)
    assert runner_2.selection.selection_id == runner.selection.selection_id


@pytest.mark.functional
def test_decode_runner_run(runner):
    """Test running a DecodingRunner"""
    runner.run(use_tqdm=False)


@pytest.mark.functional
def test_decode_runner_run_with_save_failed(runner, tmpdir):
    """Test running a DecodingRunner"""
    runner.run(use_tqdm=False, save_failed_to=tmpdir)
    selection_id = runner.selection.selection_id
    assert os.path.exists(tmpdir.join(f"{selection_id}_decode_failed.csv"))


@pytest.fixture
def ran_decode_runner() -> DecodingRunner:
    """A DecodingRunner that has run fixture for testing"""
    runner = DecodingRunner.from_file(DECODE_FILE, disable_logging=True)
    runner.run(use_tqdm=False)
    return runner


@pytest.mark.unit
def test_decode_runner_write_report(ran_decode_runner, tmpdir):
    """Test writing decode report to a file"""
    ran_decode_runner.write_decode_report(tmpdir)
    ran_decode_runner.write_decode_report(tmpdir, prefix="TEST")

    selection_id = ran_decode_runner.selection.selection_id
    assert os.path.exists(tmpdir.join(f"{selection_id}_decode_report.html"))
    assert os.path.exists(tmpdir.join("TEST_decode_report.html"))


@pytest.mark.unit
def test_decode_runner_write_stats(ran_decode_runner, tmpdir):
    """Test writing decode statistics to a file"""
    ran_decode_runner.write_decode_statistics(tmpdir)
    ran_decode_runner.write_decode_statistics(tmpdir, prefix="TEST")

    selection_id = ran_decode_runner.selection.selection_id
    assert os.path.exists(tmpdir.join(f"{selection_id}_decode_statistics.json"))
    assert os.path.exists(tmpdir.join("TEST_decode_statistics.json"))


@pytest.mark.unit
def test_decode_runner_write_cube(ran_decode_runner, tmpdir):
    """Test writing decode cube to a file"""
    ran_decode_runner.write_cube(tmpdir)
    ran_decode_runner.write_cube(tmpdir, prefix="TEST")

    selection_id = ran_decode_runner.selection.selection_id
    assert os.path.exists(tmpdir.join(f"{selection_id}_cube.csv"))
    assert os.path.exists(tmpdir.join("TEST_cube.csv"))
