"""test cases for library code"""

import os
import warnings

import pytest

from deli.configure import set_deli_data_dir
from deli.dels import BuildingBlockBarcodeSection
from deli.dels.building_block import BuildingBlockSet
from deli.dels.library import DELibrary, DELibraryPool, LibraryBuildError


set_deli_data_dir(os.path.join(os.path.dirname(__file__), "test_data", "test_deli_data_dir"))


@pytest.mark.functional
@pytest.mark.filterwarnings("ignore:missing 'smiles' column in header")
def test_load_library():
    """Test loading a DELibrary"""
    del_library = DELibrary.load("DEL004")
    assert del_library.library_id == "DEL004"
    assert del_library.scaffold is None
    assert del_library.dna_barcode_on == "DEL004_BBA"
    assert del_library.library_size == 884736
    assert del_library.enumerator is None
    assert len(del_library.bb_sets) == 3
    assert {"DEL004_BBA", "DEL004_BBB", "DEL004_BBC"} == set(
        [bb.bb_set_id for bb in del_library.bb_sets]
    )
    assert del_library.library_tag == "CCTTGGCACCCGAGAATTCCAGCCAGACCAG"
    assert del_library.barcode_schema is not None


@pytest.mark.functional
@pytest.mark.filterwarnings("ignore:missing 'smiles' column in header")
def test_load_invalid_library():
    """Test loading a DELibrary"""
    with pytest.raises(
        LibraryBuildError, match="Number of library cycles does not match barcode schema cycles"
    ):
        DELibrary.load("DEL004_invalid")


@pytest.mark.functional
def test_load_library_no_smiles():
    """Test loading a DELibrary without SMILES"""
    with pytest.warns(UserWarning, match="missing 'smiles' column in header"):
        DELibrary.load("DEL004")


@pytest.mark.functional
def test_load_library_smiles():
    """Test loading a DELibrary with SMILES"""
    with warnings.catch_warnings():
        DELibrary.load("DEL006")


@pytest.fixture
@pytest.mark.filterwarnings("ignore:missing 'smiles' column in header")
def del_library1() -> DELibrary:
    """An example DEL to load"""
    return DELibrary.load("DEL004")


@pytest.fixture
@pytest.mark.filterwarnings("ignore:missing 'smiles' column in header")
def del_library2() -> DELibrary:
    """An example DEL to load"""
    return DELibrary.load("DEL005")


@pytest.mark.unit
def test_bb_set_iter(del_library1: DELibrary):
    """Test iterating over building block sets in a DELibrary"""
    bb_sets = list(del_library1.iter_bb_sets())
    assert len(bb_sets) == 3
    assert all(isinstance(bb_set, BuildingBlockSet) for bb_set in bb_sets)
    assert len([bb.bb_set_id for bb in bb_sets]) == 3


@pytest.mark.unit
def test_iter_bb_barcode_sections_and_sets(del_library1: DELibrary):
    """Test iter_bb_barcode_sections_and_sets function"""
    bb_sets = list(del_library1.iter_bb_barcode_sections_and_sets())
    assert len(bb_sets) == 3
    assert all(
        [
            isinstance(bb_set, BuildingBlockSet)
            and isinstance(section, BuildingBlockBarcodeSection)
            for (section, bb_set) in bb_sets
        ]
    )


@pytest.mark.unit
def test_create_library_pool(del_library1, del_library2):
    """Test creating a DELibraryPool from multiple DELibraries"""
    with pytest.raises(LibraryBuildError, match="multiple libraries share identical `library_id`"):
        DELibraryPool([del_library1, del_library1])

    pool = DELibraryPool([del_library1, del_library2])
    assert pool.pool_size == 1769472
    assert len(pool.libraries) == 2


@pytest.fixture
def library_pool(del_library1, del_library2):
    """Fixture to create a DELibraryPool from two DELibraries"""
    return DELibraryPool([del_library1, del_library2])


@pytest.mark.unit
def test_get_library_from_pool(library_pool):
    """Test creating a DELibraryPool from multiple DELibraries"""
    with pytest.raises(KeyError, match="cannot find library with id"):
        library_pool.get_library("FAKE_LIBRARY_ID")
    library = library_pool.get_library("DEL004")
    assert library.library_id == "DEL004"
