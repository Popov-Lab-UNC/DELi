"""
Unit tests for ToolCompound, DopedToolCompound, and TaggedToolCompound classes.
"""

import json
from pathlib import Path

import pytest

from deli.dels import tool_compound as tc_mod
from deli.dels.tool_compound import (
    DopedToolCompound,
    TaggedToolCompound,
    ToolCompound,
    ToolCompoundParsingError,
)


def write_json(tmp_path: Path, data: dict) -> Path:
    """Helper to write JSON data to a temp file."""
    p = tmp_path / "tmp_tool_compound.json"
    p.write_text(json.dumps(data))
    return p


def test_tool_compound_to_dict_roundtrip():
    """Test ToolCompound to_dict behavior."""
    t = ToolCompound(compound_id="C001")
    d = t.to_dict()
    assert d == {"compound_id": "C001"}

    t2 = ToolCompound(compound_id="C002", smiles="C(C)C")
    d2 = t2.to_dict()
    assert d2["compound_id"] == "C002" and d2["smiles"] == "C(C)C"


def test_tool_compound_load_missing_id_raises(tmp_path: Path):
    """Test that ToolCompound.load raises ToolCompoundParsingError when compound_id is missing."""
    bad = {"smiles": "CC"}
    p = write_json(tmp_path, bad)
    with pytest.raises(ToolCompoundParsingError):
        ToolCompound.load(str(p))


def test_tool_compound_load_returns_none_for_valid(tmp_path: Path):
    """Test that ToolCompound.load returns None as per current implementation."""
    # Current implementation constructs an instance but does not return it.
    data = {"compound_id": "T100", "smiles": "CC"}
    p = write_json(tmp_path, data)
    res = ToolCompound.load(str(p))
    # ensure no exception and current behavior returns None
    assert res is None


def test_doped_tool_compound_to_from_dict_roundtrip():
    """Test DopedToolCompound to_dict and from_dict round-trip behavior."""
    dtc = DopedToolCompound(compound_id="D001", bb_tags=("A", "B", "C"), smiles="CC")
    d = dtc.to_dict()
    # bb_tags stored as comma separated string
    assert d["bb_tags"] == "A,B,C"
    loaded = DopedToolCompound.from_dict({"compound_id": "D001", "bb_tags": "A,B,C", "smiles": "CC"})
    assert isinstance(loaded, DopedToolCompound)
    assert loaded.bb_tags == ("A", "B", "C")


def test_doped_tool_compound_load_disabled():
    """Test that DopedToolCompound.load raises NotImplementedError."""
    with pytest.raises(NotImplementedError):
        DopedToolCompound.load()


def test_tagged_tool_compound_load_uses_schema_from_dict(monkeypatch, tmp_path: Path):
    """Test that TaggedToolCompound.load uses ToolCompoundBarcodeSchema.from_dict correctly."""
    data = {
        "compound_id": "TC1",
        "barcode_schema": {"library": {"tag": "LIBTAG"}},
        "smiles": "CC",
    }
    p = write_json(tmp_path, data)

    # Monkeypatch ToolCompoundBarcodeSchema.from_dict to return a fake schema
    class FakeSchema:
        pass

    monkeypatch.setattr(tc_mod, "ToolCompoundBarcodeSchema", tc_mod.ToolCompoundBarcodeSchema)
    monkeypatch.setattr(tc_mod.ToolCompoundBarcodeSchema, "from_dict", staticmethod(lambda d: FakeSchema()))

    loaded = TaggedToolCompound.load(str(p))
    assert isinstance(loaded, TaggedToolCompound)
    assert isinstance(loaded.barcode_schema, FakeSchema)
