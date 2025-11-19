"""tests for decoding classes/functions"""

import os
import shutil

import pytest

from deli.decode.barcode_calling import HashMapCollisionError
from deli.decode.runner import DecodingRunner, DecodingRunnerResults, DecodingRunParsingError


DECODE_FILE = os.path.abspath(os.path.join(os.path.dirname(__file__), "test_data", "example_decode.yaml"))
DECODE_FILE_ERR_CORRECT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "test_data", "example_decode_error_correction.yaml")
)
DECODE_FILE_ERR_DISABLE = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "test_data", "example_decode_error_correction_disable.yaml")
)
DECODE_FILE_ERR_CORRECT_ASYM = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "test_data", "example_decode_error_correction_asym.yaml")
)
NO_SEQ_DECODE_FILE = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "test_data", "example_decode_no_seqs.yaml")
)
FASTQ_FILE = os.path.abspath(os.path.join(os.path.dirname(__file__), "test_data", "example.fastq"))


@pytest.mark.unit
def test_loading_decode_runner(tmpdir):
    """Test loading a DecodingRunner from a file"""
    os.chdir(tmpdir)
    tmpdir.mkdir("test_data")
    shutil.copyfile(FASTQ_FILE, tmpdir.join("test_data", "example.fastq"))
    runner = DecodingRunner.from_file(DECODE_FILE, disable_logging=False)
    assert isinstance(runner, DecodingRunner)
    assert os.path.exists(tmpdir.join("deli.log"))

    runner = DecodingRunner.from_file(DECODE_FILE, [FASTQ_FILE], ignore_decode_seqs=True, disable_logging=True)
    assert isinstance(runner, DecodingRunner)

    runner = DecodingRunner.from_file(NO_SEQ_DECODE_FILE, [FASTQ_FILE], ignore_decode_seqs=False, disable_logging=True)
    assert isinstance(runner, DecodingRunner)

    with pytest.raises(DecodingRunParsingError, match="`fastq_files` cannot be provided when"):
        DecodingRunner.from_file(DECODE_FILE, [FASTQ_FILE])

    with pytest.raises(DecodingRunParsingError, match="``ignore_decode_seqs`` is True, but no `fastq_files`"):
        DecodingRunner.from_file(DECODE_FILE, ignore_decode_seqs=True)

    with pytest.raises(DecodingRunParsingError, match="found neither"):
        DecodingRunner.from_file(NO_SEQ_DECODE_FILE)


@pytest.mark.functional
def test_loading_decoder_runner_error_corrections(tmpdir):
    """Test loading a DecodingRunner with error corrections"""
    # check that non-hamming 1 distance fails
    with pytest.raises(HashMapCollisionError, match="collision detected"):
        DecodingRunner.from_file(DECODE_FILE_ERR_CORRECT, disable_logging=True)

    # check that disabling does not raise error
    runner = DecodingRunner.from_file(
        DECODE_FILE_ERR_DISABLE,
        disable_logging=True,
    )
    assert isinstance(runner, DecodingRunner)

    # check that asymmetric error correction works
    runner = DecodingRunner.from_file(
        DECODE_FILE_ERR_CORRECT_ASYM,
        disable_logging=True,
    )
    assert isinstance(runner, DecodingRunner)


@pytest.fixture
def runner() -> DecodingRunner:
    """A DecodingRunner fixture for testing"""
    return DecodingRunner.from_file(DECODE_FILE, disable_logging=True)


@pytest.mark.unit
def test_decoding_runner_to_file(runner, tmpdir):
    """Test saving a DecodingRunner to a file"""
    runner.to_file(tmpdir.join("decode_test.yaml"))
    runner.to_file(os.path.join("C:\\Users\\James\\Downloads", "decode_test.yaml"))
    runner_2 = DecodingRunner.from_file(tmpdir.join("decode_test.yaml"))

    assert isinstance(runner_2, DecodingRunner)
    assert runner_2.selection.selection_id == runner.selection.selection_id


@pytest.mark.functional
@pytest.mark.filterwarnings("ignore:.*UMI length.*:UserWarning")
def test_decode_runner_run(runner):
    """Test running a DecodingRunner"""
    runner.run(use_tqdm=False)


@pytest.mark.functional
@pytest.mark.filterwarnings("ignore:.*UMI length.*:UserWarning")
def test_decode_runner_run_with_save_failed(runner, tmpdir):
    """Test running a DecodingRunner"""
    runner.run(use_tqdm=False, save_failed_to=tmpdir)
    selection_id = runner.selection.selection_id
    assert os.path.exists(tmpdir.join(f"{selection_id}_decode_failed.tsv"))


@pytest.fixture
@pytest.mark.filterwarnings("ignore:.*UMI length.*:UserWarning")
def decode_results() -> DecodingRunnerResults:
    """A DecodingRunner that has run fixture for testing"""
    runner = DecodingRunner.from_file(DECODE_FILE, disable_logging=True)
    return runner.run(use_tqdm=False)


@pytest.mark.unit
@pytest.mark.filterwarnings("ignore:.*UMI length.*:UserWarning")
def test_decode_runner_write_report(decode_results, tmpdir):
    """Test writing decode report to a file"""
    decode_results.write_decode_report(tmpdir.join("TEST_decode_report.html"))
    assert os.path.exists(tmpdir.join("TEST_decode_report.html"))


@pytest.mark.unit
@pytest.mark.filterwarnings("ignore:.*UMI length.*:UserWarning")
def test_decode_runner_write_stats(decode_results, tmpdir):
    """Test writing decode statistics to a file"""
    decode_results.write_decode_statistics(tmpdir.join("TEST_decode_statistics.json"))
    assert os.path.exists(tmpdir.join("TEST_decode_statistics.json"))


@pytest.mark.unit
@pytest.mark.filterwarnings("ignore:.*UMI length.*:UserWarning")
def test_decode_runner_write_cube(decode_results, tmpdir):
    """Test writing decode cube to a file"""
    decode_results.write_cube(tmpdir.join("TEST_cube.csv"))
    assert os.path.exists(tmpdir.join("TEST_cube.csv"))
    _header = open(tmpdir.join("TEST_cube.csv"), "r").readline().strip().split(",")
    assert "DEL_ID" in _header
    assert "UMI_CORRECTED_COUNT" in _header

    with pytest.warns(UserWarning, match="Some libraries are missing enumerators"):
        decode_results.write_cube(tmpdir.join("TEST_cube.csv"), enumerate_smiles=True)
        _header = open(tmpdir.join("TEST_cube.csv"), "r").readline().strip().split(",")
        assert "SMILES" in _header

    with pytest.warns(UserWarning, match="Some libraries are missing building block smiles"):
        decode_results.write_cube(tmpdir.join("TEST_cube.csv"), include_bb_smi_cols=True)
        _header = open(tmpdir.join("TEST_cube.csv"), "r").readline().strip().split(",")
        assert "BB1_SMILES" in _header
