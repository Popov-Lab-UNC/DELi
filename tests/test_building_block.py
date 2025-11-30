"""test for building block related classes and functions"""

import pytest

from deli.dels.building_block import (
    BuildingBlock,
    BuildingBlockSet,
    BuildingBlockSetError,
    TaggedBuildingBlock,
    TaggedBuildingBlockSet,
    generate_building_block_subset_id,
    parse_building_block_subset_id,
)


def test_generate_and_parse_subset_id_roundtrip():
    """Round-trip generate and parse a building-block subset id."""
    s = generate_building_block_subset_id("set1", "subA")
    assert s == "set1:::subA"
    a, b = parse_building_block_subset_id(s)
    assert a == "set1" and b == "subA"

    with pytest.raises(ValueError):
        parse_building_block_subset_id("not-a-valid-id")


def test_building_block_equality_and_subset_behavior():
    """Test equality of BuildingBlock and behavior when not in a subset."""
    bb1 = BuildingBlock("A", smiles="C")
    bb2 = BuildingBlock("A", smiles="C")
    bb3 = BuildingBlock("B", smiles="C")

    assert bb1 == bb2
    assert bb1 != bb3
    assert not bb1.in_subset()
    with pytest.raises(RuntimeError):
        bb1.get_subset_id()


def test_tagged_building_block_tags_list():
    """Ensure TaggedBuildingBlock exposes tags as a list."""
    t = TaggedBuildingBlock("id1", tag="AAAA")
    assert isinstance(t.tags, list)
    assert t.tags == ["AAAA"]


def test_building_blockset_load_from_csv_no_subsets(tmp_path):
    """Load a BuildingBlockSet CSV without subset ids and exercise lookups."""
    csv = "id,smiles\nA,C\nB,CC\n"
    p = tmp_path / "bb_no_subsets.csv"
    p.write_text(csv)

    s = BuildingBlockSet.load_from_csv(str(p))
    assert len(s) == 2
    assert s.get_bb_by_id("A").bb_id == "A"
    # request a missing id returns None or raises based on flag
    assert s.get_bb_by_id("X") is None
    with pytest.raises(KeyError):
        s.get_bb_by_id("X", fail_on_missing=True)

    # since no subsets present, asking for the set id returns all bbs
    all_bbs = s.get_bb_subset(s.bb_set_id)
    assert len(all_bbs) == 2


def test_building_blockset_with_subsets(tmp_path):
    """Load a BuildingBlockSet CSV with subset ids and verify subset helpers."""
    csv = "id,smiles,subset_id\nA,C,s1\nB,CC,s1\nC,CCC,s2\n"
    p = tmp_path / "bb_with_subsets.csv"
    p.write_text(csv)

    s = BuildingBlockSet.load_from_csv(str(p))
    # subsets
    s1 = s.get_bb_subset("s1")
    assert len(s1) == 2
    # get_bb_subsets returns full ids (set:::subset)
    subs = s.get_bb_subsets()
    assert generate_building_block_subset_id(s.bb_set_id, "s1") in subs

    # get_subset_with_bb_id
    assert s.get_subset_with_bb_id("A") == "s1"
    assert s.get_subset_with_bb_id("A", as_bb_subset_id=True).startswith(s.bb_set_id + ":::")


def test_mixed_subset_presence_raises(tmp_path):
    """CSV with mixed presence/absence of subset_id should raise an error."""
    csv = "id,smiles,subset_id\nA,C,s1\nB,CC,\n"
    p = tmp_path / "bb_mixed.csv"
    p.write_text(csv)

    with pytest.raises(BuildingBlockSetError):
        BuildingBlockSet.load_from_csv(str(p))


def test_duplicate_id_conflicting_smiles_raises(tmp_path):
    """Duplicate IDs with conflicting SMILES should raise a BuildingBlockSetError."""
    csv = "id,smiles\nA,C\nA,CC\n"
    p = tmp_path / "bb_dup_smiles.csv"
    p.write_text(csv)

    with pytest.raises(BuildingBlockSetError):
        BuildingBlockSet.load_from_csv(str(p))


def test_tagged_building_blockset_csv_and_search(tmp_path):
    """Load TaggedBuildingBlockSet from CSV and test tag-based lookup behavior."""
    csv = "id,smiles,tag\nT1,C,AAAA\nT2,CC,TTTT\n"
    p = tmp_path / "tagged.csv"
    p.write_text(csv)

    s = TaggedBuildingBlockSet.load_from_csv(str(p))
    assert len(s) == 2
    # search tags
    assert s.search_tags("AAAA").bb_id == "T1"
    assert s.search_tags("XXXX") is None
    with pytest.raises(KeyError):
        s.search_tags("XXXX", fail_on_missing=True)


def test_tagged_building_blockset_tag_length_mismatch_raises(tmp_path):
    """Tagged CSV with mismatched tag lengths should raise a BuildingBlockSetError."""
    csv = "id,smiles,tag\nT1,C,AAA\nT2,CC,TTTT\n"
    p = tmp_path / "tagged_mismatch.csv"
    p.write_text(csv)

    with pytest.raises(BuildingBlockSetError):
        TaggedBuildingBlockSet.load_from_csv(str(p))
