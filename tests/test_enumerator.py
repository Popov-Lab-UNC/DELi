"""tests for enumerator.py classes/functions"""

import os

import pytest

from deli.dels import DELEnumerator, DELibrary
from deli.dels.enumerator import FAILED_ENUMERATION_STR


@pytest.mark.unit
def test_enumerator_load_from_file():
    """Test loading enumerator from file"""
    DELEnumerator.load("DEL006")


@pytest.fixture
def del_enumerator() -> DELEnumerator:
    """An example DEL to load"""
    return DELEnumerator.load("DEL006")


@pytest.fixture
def del_library() -> DELibrary:
    """An example DEL to load"""
    return DELibrary.load("DEL006")


@pytest.mark.unit
def test_enumerate_from_library(del_library: DELibrary, tmpdir):
    """Test is enumerate from library object to file is working"""
    file_path = os.path.join(tmpdir, "enumerate_test.csv")
    del_library.enumerate_library_to_file(file_path, use_tqdm=False)

    lines = open(file_path, "r").readlines()
    smiles_col = lines[0].strip().split(",").index("SMILES")
    compounds = [line.strip().split(",")[smiles_col] for line in lines[1:]]

    assert all([compound != FAILED_ENUMERATION_STR for compound in compounds])

    num_compounds = len(compounds)
    assert num_compounds == del_library.library_size


@pytest.mark.unit
def test_enumerator_enumerate(del_enumerator):
    """Test using the enumerator to enumerate a library"""
    del_enumerator.enumerate(use_tqdm=False)


@pytest.mark.unit
def test_enumerator_enumerate_bb(del_enumerator):
    """Test using the enumerator to enumerate a single compound"""
    bb_map = {"DEL006_BBA": "A003", "DEL006_BBB": "B001", "DEL006_BBC": "C001"}
    compound = del_enumerator.get_enumerated_compound_from_bb_ids(bb_map)

    assert (
        compound.smi == "[H]N(C(=O)c1ccc2c(c1)N=C(c1ccc3cccnc3c1)N2Cc1ccc"
        "(-n2cccn2)cc1)[C@@H](Cc1ccc(I)cc1)C([15NH2])=O"
    )
