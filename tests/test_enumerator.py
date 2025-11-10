"""tests for enumerator.py classes/functions"""

import os

import pytest

from deli.dels.library import DELibrary, Library
from deli.enumeration.enumerator import EnumeratedDELCompound


# filter reaction warning
pytestmark = pytest.mark.filterwarnings("ignore::deli.enumeration.reaction.ReactionWarning")


@pytest.mark.unit
def test_enumerator_load_from_file():
    """Test loading enumerator from file"""
    lib = Library.load("DEL006")
    assert lib.can_enumerate()


@pytest.fixture
def del_enumerator() -> Library:
    """An example DEL to load"""
    return Library.load("DEL006")


@pytest.fixture
def del_library() -> DELibrary:
    """An example DEL to load"""
    return DELibrary.load("DEL006")


@pytest.fixture
def missing_enumerator() -> DELibrary:
    """An example DEL to load"""
    return DELibrary.load("DEL005")


@pytest.mark.unit
def test_enumerate_from_library(del_library: DELibrary, tmpdir):
    """Test is enumerate from library object to file is working"""
    file_path = os.path.join(tmpdir, "enumerate_test.csv")
    del_library.enumerate_to_file(file_path, use_tqdm=False)

    lines = open(file_path, "r").readlines()
    smiles_col = lines[0].strip().split("\t").index("SMILES")
    compounds = [line.strip().split("\t")[smiles_col] for line in lines[1:]]

    num_compounds = len(compounds)
    assert num_compounds == del_library.library_size


@pytest.mark.unit
def test_enumerator_enumerate(del_enumerator, missing_enumerator):
    """Test using the enumerator to enumerate a library"""
    compounds = del_enumerator.enumerate(use_tqdm=False)
    _count = 0
    for compound in compounds:
        assert isinstance(compound, EnumeratedDELCompound)
        _count += 1
    assert _count == del_enumerator.enumerator.get_enumeration_size()

    # test no enumerator fails
    with pytest.raises(RuntimeError):
        missing_enumerator.enumerate(use_tqdm=False).__next__()


@pytest.mark.unit
def test_enumerator_enumerate_bb(del_enumerator):
    """Test using the enumerator to enumerate a single compound"""
    bb_ids = ["A003", "B001", "C001"]
    compound = del_enumerator.enumerate_by_bb_ids(bb_ids)

    assert (
        compound.smi == "[H]N(C(=O)c1ccc2c(c1)N=C(c1ccc3cccnc3c1)N2Cc1ccc"
        "(-n2cccn2)cc1)[C@@H](Cc1ccc(I)cc1)C([15NH2])=O"
    )
