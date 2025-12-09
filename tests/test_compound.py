"""Tests for compound module."""

import pytest

from deli.dels.building_block import BuildingBlock, BuildingBlockSet
from deli.dels.combinatorial import CombinatorialLibrary
from deli.dels.library import LibraryCollection
from deli.dels.compound import DELCompound, DELCompoundRaw


@pytest.fixture
def library():
    """Fixture for a library object."""
    library = CombinatorialLibrary(
        library_id="lib1",
        bb_sets=[
            BuildingBlockSet(
                bb_set_id="set1",
                building_blocks=[
                    BuildingBlock(bb_id="bb1"),
                    BuildingBlock(bb_id="bb2"),
                ],
            ),
            BuildingBlockSet(
                bb_set_id="set2",
                building_blocks=[
                    BuildingBlock(bb_id="bb3"),
                    BuildingBlock(bb_id="bb4"),
                ],
            ),
        ],
    )
    return library


def test_del_compound_initialization(library):
    """Test initialization of DELCompound."""
    building_blocks = [BuildingBlock(bb_id="bb1"), BuildingBlock(bb_id="bb3")]

    compound = DELCompound(library=library, building_blocks=building_blocks)

    assert compound.library.library_id == "lib1"
    assert len(compound.building_blocks) == 2
    assert compound.building_blocks[0].bb_id == "bb1"
    assert compound.building_blocks[1].bb_id == "bb3"
    assert compound.library_id == "lib1"
    assert compound.building_block_ids == ["bb1", "bb3"]


def test_del_compound_raw_load_compound(library):
    """Test loading DELCompound from DELCompoundRaw."""
    raw_compound = DELCompoundRaw(library_id="lib1", building_blocks_ids=["bb1", "bb3"])

    compound = raw_compound.load_compound(collection=LibraryCollection([library]))

    assert isinstance(compound, DELCompound)
    assert compound.library == library
    assert compound.library.library_id == "lib1"
    assert len(compound.building_blocks) == 2
    assert compound.building_blocks[0].bb_id == "bb1"
    assert compound.building_blocks[1].bb_id == "bb3"


def test_to_raw_compound(library):
    """Test converting DELCompound to DELCompoundRaw."""
    building_blocks = [BuildingBlock(bb_id="bb1"), BuildingBlock(bb_id="bb3")]

    compound = DELCompound(library=library, building_blocks=building_blocks)

    raw_compound = compound.to_raw()

    assert isinstance(raw_compound, DELCompoundRaw)
    assert raw_compound.library_id == "lib1"
    assert raw_compound.building_block_ids == ["bb1", "bb3"]


def test_to_dict(library):
    """Test converting DELCompound to dictionary."""
    building_blocks = [BuildingBlock(bb_id="bb1"), BuildingBlock(bb_id="bb3")]

    compound = DELCompound(library=library, building_blocks=building_blocks)

    compound_dict = compound.to_dict()

    assert isinstance(compound_dict, dict)
    assert compound_dict["library_id"] == "lib1"
    assert compound_dict["building_block_ids"] == ["bb1", "bb3"]
