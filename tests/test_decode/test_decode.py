"""tests for deli.decode.decoder module"""

import json
import os
import tempfile
from collections import defaultdict
from copy import deepcopy
from unittest.mock import Mock

import pytest

from deli.decode.barcode_calling import ValidCall
from deli.decode.base import FailedDecodeAttempt
from deli.decode.decoder import (
    DecodedDELCompound,
    DecodedToolCompound,
    DecodeStatistics,
    DecodingSettings,
    SelectionDecoder,
)
from deli.dels.barcode import DELBarcodeSchema, VariableBarcodeSection
from deli.dels.building_block import TaggedBuildingBlock, TaggedBuildingBlockSet
from deli.dels.combinatorial import DELibrary
from deli.dels.compound import DELCompound
from deli.dels.tool_compound import ToolCompound
from deli.enumeration.enumerator import EnumerationRunError


class TestDecodeStatistics:
    """Test cases for the DecodeStatistics class"""

    def test_initialization(self):
        """Test that DecodeStatistics initializes with correct default values"""
        stats = DecodeStatistics()

        assert stats.num_seqs_read == 0
        assert isinstance(stats.num_seqs_decoded_per_lib, defaultdict)
        assert len(stats.num_seqs_decoded_per_lib) == 0

        assert stats.num_failed_too_short == 0
        assert stats.num_failed_too_long == 0
        assert stats.num_failed_alignment == 0
        assert stats.num_failed_library_call == 0
        assert stats.num_failed_building_block_call == 0
        assert stats.num_failed_ambiguous_building_block_call == 0
        assert stats.num_failed_umi == 0

    def test_num_seqs_decoded_property(self):
        """Test the num_seqs_decoded property"""
        stats = DecodeStatistics()

        # Empty dict should return 0
        assert stats.num_seqs_decoded == 0

        # Add some values
        stats.num_seqs_decoded_per_lib["lib1"] = 10
        stats.num_seqs_decoded_per_lib["lib2"] = 20
        stats.num_seqs_decoded_per_lib["lib3"] = 0
        assert stats.num_seqs_decoded == 30

    def test_add_operator(self):
        """Test adding two DecodeStatistics objects"""
        stats1 = DecodeStatistics()
        stats1.num_seqs_read = 100
        stats1.num_seqs_decoded_per_lib["lib1"] = 50
        stats1.num_failed_too_short = 5
        stats1.num_failed_too_long = 3
        stats1.num_failed_alignment = 2
        stats1.num_failed_library_call = 1
        stats1.num_failed_building_block_call = 4
        stats1.num_failed_ambiguous_building_block_call = 6
        stats1.num_failed_umi = 7

        stats2 = DecodeStatistics()
        stats2.num_seqs_read = 200
        stats2.num_seqs_decoded_per_lib["lib1"] = 30
        stats2.num_seqs_decoded_per_lib["lib2"] = 40
        stats2.num_failed_too_short = 10
        stats2.num_failed_too_long = 6
        stats2.num_failed_alignment = 4
        stats2.num_failed_library_call = 2
        stats2.num_failed_building_block_call = 8
        stats2.num_failed_ambiguous_building_block_call = 12
        stats2.num_failed_umi = 14

        result = stats1 + stats2

        assert result.num_seqs_read == 300
        assert result.num_seqs_decoded_per_lib["lib1"] == 80
        assert result.num_seqs_decoded_per_lib["lib2"] == 40
        assert result.num_failed_too_short == 15
        assert result.num_failed_too_long == 9
        assert result.num_failed_alignment == 6
        assert result.num_failed_library_call == 3
        assert result.num_failed_building_block_call == 12
        assert result.num_failed_ambiguous_building_block_call == 18
        assert result.num_failed_umi == 21

    def test_add_operator_type_error(self):
        """Test that adding with non-DecodeStatistics raises TypeError"""
        stats = DecodeStatistics()

        with pytest.raises(TypeError):
            stats + "not a stats object"

    def test_str_representation(self):
        """Test string representation of DecodeStatistics"""
        stats = DecodeStatistics()
        stats.num_seqs_read = 100
        stats.num_failed_too_short = 5

        str_repr = str(stats)
        lines = str_repr.split("\n")

        # Should contain key=value pairs
        assert "num_seqs_read=100" in lines
        assert "num_failed_too_short=5" in lines

    def test_repr_representation(self):
        """Test repr representation of DecodeStatistics"""
        stats = DecodeStatistics()
        stats.num_seqs_read = 100
        stats.num_failed_too_short = 5

        repr_str = repr(stats)
        parts = repr_str.split("; ")

        # Should contain key=value pairs separated by '; '
        assert any("num_seqs_read=100" in part for part in parts)
        assert any("num_failed_too_short=5" in part for part in parts)

    def test_to_file_and_from_file(self):
        """Test saving to file and loading from file"""
        stats = DecodeStatistics()
        stats.num_seqs_read = 1000
        stats.num_seqs_decoded_per_lib["lib1"] = 500
        stats.num_seqs_decoded_per_lib["lib2"] = 300
        stats.num_failed_too_short = 50
        stats.num_failed_too_long = 30
        stats.num_failed_alignment = 20
        stats.num_failed_library_call = 10
        stats.num_failed_building_block_call = 40
        stats.num_failed_ambiguous_building_block_call = 60
        stats.num_failed_umi = 70

        with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".json") as tmp_file:
            tmp_path = tmp_file.name

        try:
            # Save to file
            stats.to_file(tmp_path)

            # Load from file
            loaded_stats = DecodeStatistics.from_file(tmp_path)

            # Verify all attributes match
            assert loaded_stats.num_seqs_read == stats.num_seqs_read
            assert loaded_stats.num_seqs_decoded_per_lib == stats.num_seqs_decoded_per_lib
            assert loaded_stats.num_failed_too_short == stats.num_failed_too_short
            assert loaded_stats.num_failed_too_long == stats.num_failed_too_long
            assert loaded_stats.num_failed_alignment == stats.num_failed_alignment
            assert loaded_stats.num_failed_library_call == stats.num_failed_library_call
            assert loaded_stats.num_failed_building_block_call == stats.num_failed_building_block_call
            assert (
                loaded_stats.num_failed_ambiguous_building_block_call == stats.num_failed_ambiguous_building_block_call
            )
            assert loaded_stats.num_failed_umi == stats.num_failed_umi

        finally:
            os.unlink(tmp_path)

    def test_from_file_with_missing_keys(self):
        """Test loading from file with missing keys (should use defaults)"""
        data = {
            "num_seqs_read": 100,
            "num_seqs_decoded_per_lib": {"lib1": 50},
            # Missing failure counts
        }

        with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".json") as tmp_file:
            json.dump(data, tmp_file)
            tmp_path = tmp_file.name

        try:
            loaded_stats = DecodeStatistics.from_file(tmp_path)

            assert loaded_stats.num_seqs_read == 100
            assert loaded_stats.num_seqs_decoded_per_lib["lib1"] == 50
            # Missing keys should default to 0
            assert loaded_stats.num_failed_too_short == 0
            assert loaded_stats.num_failed_too_long == 0
            assert loaded_stats.num_failed_alignment == 0
            assert loaded_stats.num_failed_library_call == 0
            assert loaded_stats.num_failed_building_block_call == 0
            assert loaded_stats.num_failed_ambiguous_building_block_call == 0
            assert loaded_stats.num_failed_umi == 0

        finally:
            os.unlink(tmp_path)


@pytest.fixture
def mock_del_library():
    """Create a mock DELibrary for testing"""
    library = Mock(spec=DELibrary)
    library.library_id = "TEST_LIB"

    # Mock enumerator that raises error for testing
    enumerator = Mock()
    enumerator.enumerate_by_bbs.side_effect = EnumerationRunError("Enumeration failed")
    library.enumerator = enumerator

    return library


@pytest.fixture
def mock_tagged_building_blocks():
    """Create mock TaggedBuildingBlock objects"""
    bb1 = Mock(spec=TaggedBuildingBlock)
    bb1.bb_id = "BB1"
    bb1.smi = "C1CCCCC1"
    bb1.has_smiles.return_value = True

    bb2 = Mock(spec=TaggedBuildingBlock)
    bb2.bb_id = "BB2"
    bb2.smi = "C1CCNCC1"
    bb2.has_smiles.return_value = True

    return [bb1, bb2]


@pytest.fixture
def mock_building_block_calls(mock_tagged_building_blocks):
    """Create ValidCall objects for building blocks"""
    calls = []
    for i, bb in enumerate(mock_tagged_building_blocks):
        call = ValidCall(obj=bb, score=float(i + 1))
        calls.append(call)
    return calls


@pytest.fixture
def mock_umi():
    """Create a mock UMI"""
    return "ATCGATCG"


@pytest.fixture
def mock_tool_compound():
    """Create a mock ToolCompound"""
    compound = Mock(spec=ToolCompound)
    compound.compound_id = "TOOL001"
    compound.smi = "CCO"
    compound.has_smiles.return_value = True
    return compound


class TestDecodedDELCompound:
    """Test cases for the DecodedDELCompound class"""

    def test_initialization(self, mock_del_library, mock_building_block_calls, mock_umi):
        """Test initialization of DecodedDELCompound"""
        decoded = DecodedDELCompound(
            library=mock_del_library, building_block_calls=mock_building_block_calls, umi=mock_umi
        )

        assert decoded.library == mock_del_library
        assert decoded.building_block_calls == mock_building_block_calls
        assert decoded.umi == mock_umi

    @pytest.fixture
    def mock_decoded_del_compound(self, mock_del_library, mock_building_block_calls, mock_umi):
        """Create a mock DecodedDELCompound for testing"""
        return DecodedDELCompound(
            library=mock_del_library, building_block_calls=mock_building_block_calls, umi=mock_umi
        )

    def test_to_compound(self, mock_decoded_del_compound, mock_del_library):
        """Test converting to DELCompound"""
        compound = mock_decoded_del_compound.to_compound()

        assert isinstance(compound, DELCompound)
        assert compound.library == mock_del_library
        assert len(compound.building_blocks) == 2
        assert compound.building_blocks[0].bb_id == "BB1"
        assert compound.building_blocks[1].bb_id == "BB2"

    def test_to_cube_row_dict(self, mock_decoded_del_compound):
        """Test converting to cube row dict"""
        row_dict = mock_decoded_del_compound.to_cube_row_dict(enumerate_compounds=False)

        assert row_dict["library_id"] == "TEST_LIB"
        assert row_dict["compound_id"] is not None  # compound_id is generated
        assert "smiles" not in row_dict
        assert row_dict["bb1"] == "BB1"
        assert row_dict["bb1_smiles"] == "C1CCCCC1"
        assert row_dict["bb2"] == "BB2"
        assert row_dict["bb2_smiles"] == "C1CCNCC1"

        # failure case with enumeration
        row_dict = mock_decoded_del_compound.to_cube_row_dict(enumerate_compounds=True)

        assert row_dict["library_id"] == "TEST_LIB"
        assert row_dict["compound_id"] is not None
        assert row_dict["smiles"] == "ENUMERATION_FAILED"
        assert row_dict["bb1"] == "BB1"
        assert row_dict["bb1_smiles"] == "C1CCCCC1"
        assert row_dict["bb2"] == "BB2"
        assert row_dict["bb2_smiles"] == "C1CCNCC1"

        # mock enumeration success case
        mocked_enum = deepcopy(mock_decoded_del_compound)
        mocked_enum.library.enumerator.enumerate_by_bbs.side_effect = None
        mocked_enum.library.enumerator.enumerate_by_bbs.return_value = "CCCOCC"
        row_dict = mocked_enum.to_cube_row_dict(enumerate_compounds=True)

        assert row_dict["library_id"] == "TEST_LIB"
        assert row_dict["compound_id"] is not None
        assert row_dict["smiles"] == "CCCOCC"
        assert row_dict["bb1"] == "BB1"
        assert row_dict["bb1_smiles"] == "C1CCCCC1"
        assert row_dict["bb2"] == "BB2"
        assert row_dict["bb2_smiles"] == "C1CCNCC1"

    def test_to_decode_res_row_dict(self, mock_decoded_del_compound):
        """Test converting to decode results row dict"""
        row_dict = mock_decoded_del_compound.get_decode_res_info()

        assert row_dict["bb_ids"] == "BB1,BB2"
        assert row_dict["bb_scores"] == "1.0,2.0"
        assert row_dict["umi"] == "ATCGATCG"

    def test_get_smiles_with_enumeration_failure(self, mock_decoded_del_compound):
        """Test get_smiles method when enumeration fails"""
        with pytest.raises(EnumerationRunError, match="Enumeration failed"):
            mock_decoded_del_compound.get_smiles()

    def test_get_library_id(self, mock_decoded_del_compound):
        """Test get_library_id method"""
        assert mock_decoded_del_compound.get_library_id() == "TEST_LIB"

    def test_get_score(self, mock_decoded_del_compound):
        """Test get_score method"""
        # Score should be sum of building block call scores
        assert mock_decoded_del_compound.get_score() == 3.0


class TestDecodedToolCompound:
    """Test cases for the DecodedToolCompound class"""

    def test_initialization(self, mock_tool_compound, mock_umi):
        """Test initialization of DecodedToolCompound"""
        decoded = DecodedToolCompound(tool_compound=mock_tool_compound, alignment_score=0.95, umi=mock_umi)

        assert decoded.tool_compound == mock_tool_compound
        assert decoded.alignment_score == 0.95

    @pytest.fixture
    def mock_decoded_tool_compound(self, mock_tool_compound, mock_umi):
        """Create a mock DecodedToolCompound for testing"""
        return DecodedToolCompound(tool_compound=mock_tool_compound, alignment_score=0.95, umi=mock_umi)

    def test_to_compound(self, mock_tool_compound, mock_decoded_tool_compound):
        """Test converting to ToolCompound"""
        compound = mock_decoded_tool_compound.to_compound()
        assert compound == mock_tool_compound

    def test_to_cube_row_dict(self, mock_decoded_tool_compound):
        """Test converting to cube row dict when tool compound has SMILES"""
        row_dict = mock_decoded_tool_compound.to_cube_row_dict()

        assert row_dict["library_id"] == "ToolCompound"
        assert row_dict["compound_id"] == "TOOL001"
        assert row_dict["smiles"] == "CCO"

        tool_compound_no_smiles = Mock(spec=ToolCompound)
        tool_compound_no_smiles.compound_id = "TOOL002"
        tool_compound_no_smiles.has_smiles.return_value = False
        no_smiles = deepcopy(mock_decoded_tool_compound)
        no_smiles.tool_compound = tool_compound_no_smiles

        row_dict = no_smiles.to_cube_row_dict()

        assert row_dict["library_id"] == "ToolCompound"
        assert row_dict["compound_id"] == "TOOL002"
        assert "smiles" not in row_dict

    def test_to_decode_res_row_dict(self, mock_decoded_tool_compound):
        """Test converting to decode results row dict"""
        row_dict = mock_decoded_tool_compound.get_decode_res_info()
        assert row_dict["umi"] == "ATCGATCG"

    def test_get_smiles_with_smiles(self, mock_decoded_tool_compound):
        """Test get_smiles method when tool compound has SMILES"""
        smiles = mock_decoded_tool_compound.get_smiles()
        assert smiles == "CCO"

        tool_compound_no_smiles = Mock(spec=ToolCompound)
        tool_compound_no_smiles.compound_id = "TOOL002"
        tool_compound_no_smiles.has_smiles.return_value = False

        mock_decoded_tool_compound.tool_compound = tool_compound_no_smiles

        with pytest.raises(RuntimeError, match="Single compound TOOL002 has no SMILES"):
            mock_decoded_tool_compound.get_smiles()

    def test_get_library_id(self, mock_decoded_tool_compound):
        """Test get_library_id method"""
        assert mock_decoded_tool_compound.get_library_id() == "TOOL001"

    def test_get_score(self, mock_decoded_tool_compound):
        """Test get_score method"""
        assert mock_decoded_tool_compound.get_score() == 0.95


def _make_a_fake_3cycle_lib(lib_id: str, lib_tag: str) -> DELibrary:
    """Create a real DELibrary mock for testing"""
    from deli.dels.barcode import (
        BuildingBlockBarcodeSection,
        LibraryBarcodeSection,
        StaticBarcodeSection,
        UMIBarcodeSection,
    )
    from deli.dels.building_block import (
        TaggedBuildingBlock,
        TaggedNullBuildingBlock,
    )

    sections = [
        StaticBarcodeSection(section_name="primer", section_tag="GTAGATCTTAGATTCGAGAGTC"),
        LibraryBarcodeSection(section_name="library", section_tag=lib_tag),
        StaticBarcodeSection(section_name="static1", section_tag="GTTTAGATCATGC"),
        BuildingBlockBarcodeSection(cycle_number=1, section_name="bb1", section_tag="NNNNNNNN", section_overhang="AG"),
        BuildingBlockBarcodeSection(cycle_number=2, section_name="bb2", section_tag="NNNNNNNN", section_overhang="GC"),
        BuildingBlockBarcodeSection(cycle_number=3, section_name="bb3", section_tag="NNNNNNNN", section_overhang="TC"),
        StaticBarcodeSection(section_name="static2", section_tag="CCTAGATTAGATTC"),
        UMIBarcodeSection(section_name="umi", section_tag="NNNNNNNNNN"),
    ]

    barcode_schema = DELBarcodeSchema(sections)

    # building block sets
    bb_set_1 = TaggedBuildingBlockSet(
        bb_set_id="bb1_TEST",
        building_blocks=[
            TaggedBuildingBlock(bb_id="BB1_A", tag="AAAAAAAA"),
            TaggedBuildingBlock(bb_id="BB1_B", tag="AAAATTTT"),
            TaggedBuildingBlock(bb_id="BB1_C", tag="CCCCCCCC"),
            TaggedBuildingBlock(bb_id="BB1_D", tag="GGGGGGGG"),
            TaggedBuildingBlock(bb_id="BB1_E", tag="AACCGGTT"),
            TaggedBuildingBlock(bb_id="BB1_F", tag="TTTTCCCC"),
            TaggedBuildingBlock(bb_id="BB1_G", tag="AGAGAGAG"),
            TaggedNullBuildingBlock(bb_id="BB1_NULL", tag="CTCTCTCT"),
        ],
    )
    bb_set_2 = TaggedBuildingBlockSet(
        bb_set_id="bb2_TEST",
        building_blocks=[
            TaggedBuildingBlock(bb_id="BB2_A", tag="AAAAAAAA"),
            TaggedBuildingBlock(bb_id="BB2_B", tag="AAAATTTT"),
            TaggedBuildingBlock(bb_id="BB2_C", tag="CCCCCCCC"),
            TaggedBuildingBlock(bb_id="BB2_D", tag="GGGGGGGG"),
            TaggedBuildingBlock(bb_id="BB2_E", tag="AACCGGTT"),
            TaggedBuildingBlock(bb_id="BB2_F", tag="TTTTCCCC"),
            TaggedBuildingBlock(bb_id="BB2_G", tag="AGAGAGAG"),
            TaggedNullBuildingBlock(bb_id="BB2_NULL", tag="CTCTCTCT"),
        ],
    )
    bb_set_3 = TaggedBuildingBlockSet(
        bb_set_id="bb3_TEST",
        building_blocks=[
            TaggedBuildingBlock(bb_id="BB3_A", tag="AAAAAAAA"),
            TaggedBuildingBlock(bb_id="BB3_B", tag="AAAATTTT"),
            TaggedBuildingBlock(bb_id="BB3_C", tag="CCCCCCCC"),
            TaggedBuildingBlock(bb_id="BB3_D", tag="GGGGGGGG"),
            TaggedBuildingBlock(bb_id="BB3_E", tag="AACCGGTT"),
            TaggedBuildingBlock(bb_id="BB3_F", tag="TTTTCCCC"),
            TaggedBuildingBlock(bb_id="BB3_G", tag="AGAGAGAG"),
            TaggedNullBuildingBlock(bb_id="BB3_NULL", tag="CTCTCTCT"),
        ],
    )

    lib = DELibrary(
        library_id=lib_id,
        barcode_schema=barcode_schema,
        bb_sets=[bb_set_1, bb_set_2, bb_set_3],
    )

    return lib


def _make_a_fake_2cycle_lib(lib_id: str, lib_tag: str) -> DELibrary:
    """Create a real DELibrary mock for testing"""
    from deli.dels.barcode import (
        BuildingBlockBarcodeSection,
        LibraryBarcodeSection,
        StaticBarcodeSection,
        UMIBarcodeSection,
    )
    from deli.dels.building_block import (
        TaggedBuildingBlock,
        TaggedNullBuildingBlock,
    )

    sections = [
        StaticBarcodeSection(section_name="primer", section_tag="GTAGATCTTAGATTCGAGAGTC"),
        LibraryBarcodeSection(section_name="library", section_tag=lib_tag),
        StaticBarcodeSection(section_name="static1", section_tag="GTTTAGATCATGC"),
        BuildingBlockBarcodeSection(cycle_number=1, section_name="bb1", section_tag="NNNNNNNN", section_overhang="AG"),
        BuildingBlockBarcodeSection(cycle_number=2, section_name="bb2", section_tag="NNNNNNNN", section_overhang="GC"),
        StaticBarcodeSection(section_name="static2", section_tag="CCTAGATTAGATTC"),
        UMIBarcodeSection(section_name="umi", section_tag="NNNNNNNNNN"),
    ]

    barcode_schema = DELBarcodeSchema(sections)

    # building block sets
    bb_set_1 = TaggedBuildingBlockSet(
        bb_set_id="bb1_TEST",
        building_blocks=[
            TaggedBuildingBlock(bb_id="BB1_A", tag="AAAAAAAA"),
            TaggedBuildingBlock(bb_id="BB1_B", tag="AAAATTTT"),
            TaggedBuildingBlock(bb_id="BB1_C", tag="CCCCCCCC"),
            TaggedBuildingBlock(bb_id="BB1_D", tag="GGGGGGGG"),
            TaggedBuildingBlock(bb_id="BB1_E", tag="AACCGGTT"),
            TaggedBuildingBlock(bb_id="BB1_F", tag="TTTTCCCC"),
            TaggedBuildingBlock(bb_id="BB1_G", tag="AGAGAGAG"),
            TaggedNullBuildingBlock(bb_id="BB1_NULL", tag="CTCTCTCT"),
        ],
    )
    bb_set_2 = TaggedBuildingBlockSet(
        bb_set_id="bb2_TEST",
        building_blocks=[
            TaggedBuildingBlock(bb_id="BB2_A", tag="AAAAAAAA"),
            TaggedBuildingBlock(bb_id="BB2_B", tag="AAAATTTT"),
            TaggedBuildingBlock(bb_id="BB2_C", tag="CCCCCCCC"),
            TaggedBuildingBlock(bb_id="BB2_D", tag="GGGGGGGG"),
            TaggedBuildingBlock(bb_id="BB2_E", tag="AACCGGTT"),
            TaggedBuildingBlock(bb_id="BB2_F", tag="TTTTCCCC"),
            TaggedBuildingBlock(bb_id="BB2_G", tag="AGAGAGAG"),
            TaggedNullBuildingBlock(bb_id="BB2_NULL", tag="CTCTCTCT"),
        ],
    )

    lib = DELibrary(
        library_id=lib_id,
        barcode_schema=barcode_schema,
        bb_sets=[bb_set_1, bb_set_2],
    )

    return lib


def _make_a_fake_tagged_tool_compound(compound_id: str, compound_tag: str):
    """Create a real ToolCompound for testing"""
    from deli.dels.barcode import (
        LibraryBarcodeSection,
        StaticBarcodeSection,
        ToolCompoundBarcodeSchema,
        ToolCompoundRefBarcodeSection,
        UMIBarcodeSection,
    )
    from deli.dels.tool_compound import TaggedToolCompound

    sections = [
        StaticBarcodeSection(section_name="primer", section_tag="GTAGATCTTAGATTCGAGAGTC"),
        LibraryBarcodeSection(section_name="library", section_tag=compound_tag),
        StaticBarcodeSection(section_name="static1", section_tag="GTTTAGATCATGC"),
        ToolCompoundRefBarcodeSection(section_name="tool_compound_ref", section_tag="AGTCACTG", section_overhang="AG"),
        VariableBarcodeSection(section_name="bb2", section_tag="NNNNNNNN", section_overhang="GC"),
        VariableBarcodeSection(section_name="bb3", section_tag="NNNNNNNN", section_overhang="TC"),
        StaticBarcodeSection(section_name="static2", section_tag="CCTAGATTAGATTC"),
        UMIBarcodeSection(section_name="umi", section_tag="NNNNNNNNNN"),
    ]

    schema = ToolCompoundBarcodeSchema(sections)

    tool_compound = TaggedToolCompound(compound_id=compound_id, barcode_schema=schema, smiles="CCO")

    return tool_compound


@pytest.fixture
def test_del_selection():
    """Fixture for a fake DELibrary"""
    import itertools

    from dnaio import SequenceRecord

    from deli.dels.combinatorial import DELibraryCollection
    from deli.dna.io import SequenceReader
    from deli.selection import SequencedSelection

    lib1 = _make_a_fake_3cycle_lib("FAKE_LIB_1", "GTAGCTAG")
    lib2 = _make_a_fake_3cycle_lib("FAKE_LIB_2", "CCTGATAC")
    lib3 = _make_a_fake_3cycle_lib("FAKE_LIB_3", "TTTTGCAA")
    lib_4 = _make_a_fake_2cycle_lib("FAKE_LIB_4", "TTGACTCG")
    tool1 = _make_a_fake_tagged_tool_compound("TOOL_1", "GATCGTAC")
    tool2 = _make_a_fake_tagged_tool_compound("TOOL_2", "CTAGCTAG")

    def bb_triple_to_dna(bb1, bb2, bb3):
        return f"{bb1.tags[0]}AG{bb2.tags[0]}GC{bb3.tags[0]}TC"

    def bb_double_to_dna(bb1, bb2):
        return f"{bb1.tags[0]}AG{bb2.tags[0]}GC"

    example_seqs = []
    for library in [lib1, lib2, lib3]:
        dna_front = f"GTGTGTGTGTAGATCTTAGATTCGAGAGTC{library.barcode_schema.library_section.section_tag}GTTTAGATCATGC"
        bb_chunks = itertools.product(library.bb_sets[0], library.bb_sets[1], library.bb_sets[2])
        dna_back = "CCTAGATTAGATTC" + "A" * 10  # static2 + umi
        example_seqs.extend([f"{dna_front}{bb_triple_to_dna(*bbs)}{dna_back}" for bbs in bb_chunks])

    for library in [lib_4]:
        dna_front = f"GTGTGTGTGTAGATCTTAGATTCGAGAGTC{library.barcode_schema.library_section.section_tag}GTTTAGATCATGC"
        bb_chunks = itertools.product(library.bb_sets[0], library.bb_sets[1])
        dna_back = "CCTAGATTAGATTC" + "A" * 10  # static2 + umi
        example_seqs.extend([f"{dna_front}{bb_double_to_dna(*bbs)}{dna_back}" for bbs in bb_chunks])

    for tool in [tool1, tool2]:
        dna_front = f"GTGTGTGTGTAGATCTTAGATTCGAGAGTC{tool.barcode_schema.library_section.section_tag}GTTTAGATCATGC"
        dna_middle = f"{tool.barcode_schema.tool_compound_ref_section.section_tag}AGNNNNNNNNGCNNNNNNNNTC"
        dna_back = "CCTAGATTAGATTC" + "A" * 10  # static2 + umi
        example_seqs.extend([f"{dna_front}{dna_middle}{dna_back}"])

    mocked_seq_reader = Mock(spec=SequenceReader)
    mocked_seq_reader.iter_seqs.return_value = iter(
        [SequenceRecord(name=f"TEST_SEQ_{i}", sequence=dna) for i, dna in enumerate(example_seqs)]
    )

    collection = DELibraryCollection(libraries=[lib1, lib2, lib3, lib_4])

    selection = SequencedSelection(
        selection_id="TEST_SELECTION",
        library_collection=collection,
        tool_compounds=[tool1, tool2],
        sequence_reader=mocked_seq_reader,
    )

    return selection


@pytest.mark.parametrize(
    "settings",
    [
        DecodingSettings(demultiplexer_mode="single", demultiplexer_algorithm="regex"),
        DecodingSettings(demultiplexer_mode="flanking", demultiplexer_algorithm="regex"),
        DecodingSettings(demultiplexer_mode="library", demultiplexer_algorithm="regex"),
        DecodingSettings(demultiplexer_mode="single", demultiplexer_algorithm="cutadapt"),
        DecodingSettings(demultiplexer_mode="flanking", demultiplexer_algorithm="cutadapt"),
        DecodingSettings(demultiplexer_mode="library", demultiplexer_algorithm="cutadapt"),
        DecodingSettings(demultiplexer_mode="single", demultiplexer_algorithm="regex", wiggle=True),
        DecodingSettings(demultiplexer_mode="flanking", demultiplexer_algorithm="regex", wiggle=True),
        DecodingSettings(demultiplexer_mode="library", demultiplexer_algorithm="regex", wiggle=True),
        DecodingSettings(demultiplexer_algorithm="full"),
    ],
)
def test_fake_decode(settings, test_del_selection):
    """Test a fake decode run using the test DEL selection"""
    decoder = SelectionDecoder(
        selection=test_del_selection,
        decode_settings=settings,
    )

    for seq in test_del_selection.sequence_reader.iter_seqs():
        decoded_read = decoder.decode_read(seq)
        assert not isinstance(decoded_read, FailedDecodeAttempt)
