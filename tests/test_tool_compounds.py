"""Unit tests for deli.dels.tool_compounds

Covers ToolCompound, DopedToolCompound, TaggedToolCompound, ToolCompoundLibrary
and TaggedToolCompoundLibrary behaviors, including loading from JSON and
searching tags.
"""

import json
from pathlib import Path

import pytest

from deli.dels import tool_compounds as tc_mod
from deli.dels.tool_compounds import (
    DopedToolCompound,
    TaggedToolCompound,
    TaggedToolCompoundLibrary,
    ToolCompound,
    ToolCompoundLibrary,
    ToolCompoundParsingError,
)


def test_tool_compound_to_from_dict_roundtrip():
    """Test ToolCompound to_dict and from_dict round-trip behavior."""
    t = ToolCompound(compound_id="C001", smiles=None)
    d = t.to_dict()
    assert d == {"compound_id": "C001"}
    t2 = ToolCompound.from_dict(d)
    assert isinstance(t2, ToolCompound)
    assert t2.compound_id == "C001"

    # with smiles
    t3 = ToolCompound(compound_id="C002", smiles="C(C)C")
    d3 = t3.to_dict()
    assert d3["compound_id"] == "C002" and d3["smiles"] == "C(C)C"
    t4 = ToolCompound.from_dict(d3)
    assert t4._smiles == "C(C)C"


def test_doped_tool_compound_to_from_dict_roundtrip():
    """Test DopedToolCompound to_dict and from_dict round-trip behavior."""
    dtc = DopedToolCompound(compound_id="D001", bb_tags=("A", "B", "C"), smiles="CC")
    d = dtc.to_dict()
    # bb_tags stored as comma separated string
    assert d["bb_tags"] == "A,B,C"
    loaded = DopedToolCompound.from_dict({"compound_id": "D001", "bb_tags": "A,B,C", "smiles": "CC"})
    assert isinstance(loaded, DopedToolCompound)
    assert loaded.bb_tags == ("A", "B", "C")


def test_tagged_tool_compound_to_from_dict_roundtrip():
    """Test TaggedToolCompound to_dict and from_dict round-trip behavior."""
    tagged = TaggedToolCompound(compound_id="T001", tag="TTTT", smiles=None)
    d = tagged.to_dict()
    assert d["compound_id"] == "T001" and d["tag"] == "TTTT"
    t2 = TaggedToolCompound.from_dict(d)
    assert isinstance(t2, TaggedToolCompound)
    assert t2.tag == "TTTT"


def write_temp_json(data, tmp_path: Path) -> Path:
    """Helper to write data as JSON to a temp file and return the path."""
    p = tmp_path / "tmp_tool_compounds.json"
    p.write_text(json.dumps(data))
    return p


def test_tool_compound_library_load_and_get(tmp_path: Path):
    """Test loading ToolCompoundLibrary from JSON and getting compounds."""
    data = {"tool_compounds": [{"compound_id": "A1"}, {"compound_id": "B2", "smiles": "CC"}]}
    p = write_temp_json(data, tmp_path)

    lib = ToolCompoundLibrary.load(str(p))
    assert isinstance(lib, ToolCompoundLibrary)
    # check compounds were created
    assert lib.get_compound("A1").compound_id == "A1"
    assert lib.get_compound("B2")._smiles == "CC"

    # missing id should raise parsing error
    bad_data = {"tool_compounds": [{"smiles": "CC"}]}
    p2 = write_temp_json(bad_data, tmp_path)
    with pytest.raises(ToolCompoundParsingError):
        ToolCompoundLibrary.load(str(p2))


def test_tagged_tool_compound_library_init_and_search():
    """Test TaggedToolCompoundLibrary initialization and tag searching."""
    # create few tagged compounds
    t1 = TaggedToolCompound(compound_id="X1", tag="TAG1")
    t2 = TaggedToolCompound(compound_id="X2", tag="TAG2")

    # fake barcode schema object - only need library_section.section_tag for ctor
    class FakeLib:
        section_tag = "LIBTAG"

    class FakeSchema:
        library_section = FakeLib()

    schema = FakeSchema()

    tagged_lib = TaggedToolCompoundLibrary("libid", [t1, t2], barcode_schema=schema)
    # search existing tag
    found = tagged_lib.search_tags("TAG1")
    assert isinstance(found, TaggedToolCompound) and found.compound_id == "X1"

    # not found returns None by default
    assert tagged_lib.search_tags("NOPE") is None

    # not found with fail_on_missing True raises KeyError and message contains library tag
    with pytest.raises(KeyError) as exc:
        tagged_lib.search_tags("NOPE", fail_on_missing=True)
    assert "LIBTAG" in str(exc.value)


def test_tagged_tool_compound_library_load_uses_schema_from_dict(monkeypatch, tmp_path: Path):
    """Test TaggedToolCompoundLibrary.load uses ToolCompoundBarcodeSchema.from_dict to get schema."""
    # Prepare JSON that TaggedToolCompoundLibrary.load expects
    data = {
        "barcode_schema": {"library": {"tag": "LIBTAG"}},
        "tool_compounds": [
            {"compound_id": "TC1", "tag": "TAG1"},
            {"compound_id": "TC2", "tag": "TAG2"},
        ],
    }
    p = write_temp_json(data, tmp_path)

    # Monkeypatch ToolCompoundBarcodeSchema.from_dict to return a fake schema with library_section
    class FakeLib:
        section_tag = "LIBTAG"

    class FakeSchema:
        library_section = FakeLib()

    monkeypatch.setattr(tc_mod, "ToolCompoundBarcodeSchema", tc_mod.ToolCompoundBarcodeSchema)
    monkeypatch.setattr(tc_mod.ToolCompoundBarcodeSchema, "from_dict", staticmethod(lambda d: FakeSchema()))

    loaded = TaggedToolCompoundLibrary.load(str(p))
    assert isinstance(loaded, TaggedToolCompoundLibrary)
    # ensure tags mapping works
    assert loaded.search_tags("TAG2", fail_on_missing=True).compound_id == "TC2"


def test_get_compound_keyerror():
    """Test ToolCompoundLibrary.get_compound raises KeyError for missing id."""
    lib = ToolCompoundLibrary("id", [])
    with pytest.raises(KeyError):
        lib.get_compound("nope")
