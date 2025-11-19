"""test cases for library code"""

import warnings

import pytest

from deli.dels.barcode import BuildingBlockBarcodeSection
from deli.dels.building_block import BuildingBlockSet
from deli.dels.library import DELibrary, DELibraryCollection, LibraryBuildError


@pytest.mark.functional
def test_load_library():
    """Test loading a DELibrary"""
    del_library = DELibrary.load("DEL004")
    assert del_library.library_id == "DEL004"
    assert del_library.scaffold is None
    assert del_library.dna_barcode_on == "DEL004_BBA"
    assert del_library.library_size == 884736
    assert not del_library.can_enumerate()
    assert len(del_library.bb_sets) == 3
    assert {"DEL004_BBA", "DEL004_BBB", "DEL004_BBC"} == set([bb.bb_set_id for bb in del_library.bb_sets])
    assert del_library.library_tag == "CCTTGGCACCCGAGAATTCCAGCCAGACCAG"
    assert del_library.barcode_schema is not None


@pytest.mark.functional
def test_load_invalid_library():
    """Test loading a DELibrary"""
    with pytest.raises(
        LibraryBuildError, match="Number of library cycles does not match observed_barcode schema cycles"
    ):
        DELibrary.load("DEL004_invalid")


@pytest.mark.functional
def test_load_library_smiles():
    """Test loading a DELibrary with SMILES"""
    with warnings.catch_warnings():
        DELibrary.load("DEL006")


@pytest.fixture
def del_library1() -> DELibrary:
    """An example DEL to load"""
    return DELibrary.load("DEL004")


@pytest.fixture
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
def test_get_section_length(del_library1: DELibrary):
    """Test for getting observed_barcode section lengths"""
    assert del_library1.barcode_schema.get_section_length("umi") == 11
    assert del_library1.barcode_schema.get_section_length("library") == 31
    assert del_library1.barcode_schema.get_section_length("preumi") == 12

    with pytest.raises(KeyError):
        del_library1.barcode_schema.get_section_length("FAKE_SECTION")


@pytest.mark.unit
def test_iter_bb_barcode_sections_and_sets(del_library1: DELibrary):
    """Test iter_bb_barcode_sections_and_sets function"""
    bb_sets = list(del_library1.iter_bb_barcode_sections_and_sets())
    assert len(bb_sets) == 3
    assert all(
        [
            isinstance(bb_set, BuildingBlockSet) and isinstance(section, BuildingBlockBarcodeSection)
            for (section, bb_set) in bb_sets
        ]
    )


@pytest.mark.unit
def test_create_library_pool(del_library1, del_library2):
    """Test creating a DELibraryCollection from multiple DELibraries"""
    with pytest.raises(LibraryBuildError, match="multiple libraries share identical `library_id`"):
        DELibraryCollection([del_library1, del_library1])

    pool = DELibraryCollection([del_library1, del_library2])
    assert pool.collection_size == 1769472
    assert len(pool.libraries) == 2


@pytest.fixture
def library_pool(del_library1, del_library2):
    """Fixture to create a DELibraryCollection from two DELibraries"""
    return DELibraryCollection([del_library1, del_library2])


@pytest.mark.unit
def test_get_library_from_pool(library_pool):
    """Test creating a DELibraryCollection from multiple DELibraries"""
    with pytest.raises(KeyError, match="cannot find library with id"):
        library_pool.get_library("FAKE_LIBRARY_ID")
    library = library_pool.get_library("DEL004")
    assert library.library_id == "DEL004"
