"""Tests for the barcode schema and parsing functionality

These tests exercise the classes and functions defined in
`src/deli/dels/barcode.py`.
"""

import pytest

from deli.dels.barcode import (
    BarcodeSchemaError,
    BarcodeSection,
    BuildingBlockBarcodeSection,
    DecodeableBarcodeSection,
    DELBarcodeSchema,
    LibraryBarcodeSection,
    MixedBarcodeSection,
    PrimerBarcodeSection,
    StaticBarcodeSection,
    ToolCompoundBarcodeSchema,
    ToolCompoundRefBarcodeSection,
    UMIBarcodeSection,
    VariableBarcodeSection,
    _parse_sections_from_dict,
)


@pytest.fixture
def del_schema():
    """Fixture for a DELBarcodeSchema object."""
    sections = [
        LibraryBarcodeSection("library", "AGCT"),
        StaticBarcodeSection("static1", "AGCT"),
        BuildingBlockBarcodeSection(1, "bb1", "NNNN"),
        VariableBarcodeSection("var1", "NNN"),
        MixedBarcodeSection("mix1", "ANCT"),
        BuildingBlockBarcodeSection(2, "bb2", "NNNN"),
        PrimerBarcodeSection("primer1", "AGCT"),
        UMIBarcodeSection("umi", "NNNNNN"),
    ]
    return DELBarcodeSchema(sections)


class TestBarcodeSections:
    """
    Test suite for BarcodeSection and its subclasses.

    Test cases only exist for barcode sections that build upon the API of BarcodeSection.
    """

    def test_barcode_section_rules_and_api(self):
        """Test BarcodeSection for invalid nucleotide handling and basic API."""
        # valid basic section
        sec = BarcodeSection("barcode", "ATG")
        assert sec.section_name == "barcode"
        assert sec.section_tag == "ATG"
        assert sec.get_dna_sequence() == "ATG"
        assert len(sec) == 3
        assert not sec.has_overhang()
        sec = BarcodeSection("barcode", "ATG", section_overhang="GT")
        assert sec.has_overhang()
        assert sec.get_dna_sequence() == "ATGGT"
        assert len(sec) == 5

        # overhang with invalid nucleotide
        with pytest.raises(BarcodeSchemaError):
            BarcodeSection("bad", "ATG", section_overhang="BX")

        # invalid nucleotide in tag
        with pytest.raises(BarcodeSchemaError):
            BarcodeSection("bad2", "ATX")

        # empty tag not allowed
        with pytest.raises(BarcodeSchemaError):
            BarcodeSection("bad3", "")

        # identical tag and overhang detection
        s1 = BarcodeSection("s1", "AAA", section_overhang="TT")
        s2 = BarcodeSection("s2", "AAA", section_overhang="TT")
        assert s1.has_identical_tag(s2)

        # not equal if name differs
        assert s1 != s2

        # equality requires same name, tag and overhang
        s3 = BarcodeSection("s1", "AAA", section_overhang="TT")
        assert s1 == s3

    def test_variable_section_rules(self):
        """Test VariableBarcodeSection specific rules."""
        VariableBarcodeSection(
            "var",
            "NNN",
        )  # valid
        # must contain only N
        with pytest.raises(BarcodeSchemaError):
            VariableBarcodeSection("var", "ANN")

    def test_static_section_rules(self):
        """Test StaticBarcodeSection rules."""
        StaticBarcodeSection("stat", "AGCT")

        # static cannot contain N
        with pytest.raises(BarcodeSchemaError):
            StaticBarcodeSection("stat", "ANGCT")

    def test_mixed_section_rules(self):
        """Test MixedBarcodeSection specific rules."""
        # mixed allows any combination of A, T, C, G and N
        MixedBarcodeSection("mix1", "ANCT")
        MixedBarcodeSection("mix2", "NNNN")
        MixedBarcodeSection("mix3", "AGCT")

    def test_decodable_barcode_section_rules(self):
        """Test ToolCompoundRefBarcodeSection specific rules."""
        # valid decodable section
        dec = DecodeableBarcodeSection("dec", "NNNN", error_correction_mode="test")
        assert dec.error_correction_mode == "test"
        assert isinstance(dec, VariableBarcodeSection)

    def test_building_block_section(self):
        """Test BuildingBlockBarcodeSection specific rules."""
        # valid building block section
        bb = BuildingBlockBarcodeSection(1, "bb1", "NNNN")
        assert bb.cycle_number == 1

        # cycle number must be positive integer
        with pytest.raises(BarcodeSchemaError):
            BuildingBlockBarcodeSection(0, "bb2", "NNNN")
        with pytest.raises(BarcodeSchemaError):
            BuildingBlockBarcodeSection(-1, "bb3", "NNNN")

        # this will work because it should covert the int
        BuildingBlockBarcodeSection("1", "bb4", "NNNN")

    def test_building_block_and_tool_compound_and_umi_classes(self):
        """Test BuildingBlockBarcodeSection, ToolCompoundBarcodeSection, and UMIBarcodeSection."""
        # Building block needs cycle number first
        bb1 = BuildingBlockBarcodeSection(1, "bb1", "NNNN")
        bb2 = BuildingBlockBarcodeSection(2, "bb2", "NNNN")
        assert bb1.cycle_number == 1
        assert bb2.cycle_number == 2

        # umi
        umi = UMIBarcodeSection("umi", "NNNNNN")
        assert umi.section_name == "umi"


def test_parse_sections_from_dict():
    """Test parsing barcode sections from a dict and constructing DELBarcodeSchema from dict."""
    data = {
        "library": {"tag": "AGCT", "error_correction": "disable"},  # ignores extra parts
        "tool": {"tag": "ACGT", "overhang": "TT"},
        "bb1": {"tag": "NNNN", "error_correction": "disable"},
        "bb2": {"tag": "NNNN"},
        "umi": {"tag": "NNNNNN"},
        "primer_abc": {"tag": "AGCT", "JUNK": "IGNORE_ME"},
        "tool_compound_ref": {"tag": "AGCT"},
        "static1": {"tag": "TTTT"},
        "variable1": {"tag": "NNN"},
        "mixed1": {"tag": "ANCT"},
    }

    section_types = [
        LibraryBarcodeSection,
        LibraryBarcodeSection,
        BuildingBlockBarcodeSection,
        BuildingBlockBarcodeSection,
        UMIBarcodeSection,
        PrimerBarcodeSection,
        ToolCompoundRefBarcodeSection,
        StaticBarcodeSection,
        VariableBarcodeSection,
        MixedBarcodeSection,
    ]

    # parsing is exposed publicly via this class method
    sections = _parse_sections_from_dict(data)

    for section, sec_type in zip(sections, section_types, strict=False):
        assert isinstance(section, sec_type), f"Expected {sec_type}, got {type(section)}"


@pytest.fixture
def library_barcode_section():
    """Fixture for a LibraryBarcodeSection instance."""
    return LibraryBarcodeSection("library", "AGCT")


@pytest.fixture
def static_barcode_section():
    """Fixture for a StaticBarcodeSection instance."""
    return StaticBarcodeSection("static1", "AGCT")


@pytest.fixture
def building_block_barcode_section_cycle1():
    """Fixture for a BuildingBlockBarcodeSection (cycle 1)."""
    return BuildingBlockBarcodeSection(1, "bb1", "NNNN")


@pytest.fixture
def building_block_barcode_section_cycle2():
    """Fixture for a BuildingBlockBarcodeSection (cycle 2)."""
    return BuildingBlockBarcodeSection(2, "bb2", "NNNN")


@pytest.fixture
def variable_barcode_section():
    """Fixture for a VariableBarcodeSection instance."""
    return VariableBarcodeSection("var1", "NNN")


@pytest.fixture
def mixed_barcode_section():
    """Fixture for a MixedBarcodeSection instance."""
    return MixedBarcodeSection("mix1", "ANCT")


@pytest.fixture
def primer_barcode_section():
    """Fixture for a PrimerBarcodeSection instance."""
    return PrimerBarcodeSection("primer1", "AGCT")


@pytest.fixture
def umi_barcode_section():
    """Fixture for a UMIBarcodeSection instance."""
    return UMIBarcodeSection("umi", "NNNNNN")


@pytest.fixture
def tool_compound_ref_barcode_section():
    """Fixture for a ToolCompoundRefBarcodeSection instance."""
    return ToolCompoundRefBarcodeSection("tool_compound_ref", "AGCT")


@pytest.fixture
def barcode_schema(
    library_barcode_section,
    static_barcode_section,
    building_block_barcode_section_cycle1,
    variable_barcode_section,
    building_block_barcode_section_cycle2,
    umi_barcode_section,
):
    """Fixture for a DELBarcodeSchema instance."""
    return DELBarcodeSchema(
        [
            library_barcode_section,
            static_barcode_section,
            building_block_barcode_section_cycle1,
            variable_barcode_section,
            building_block_barcode_section_cycle2,
            umi_barcode_section,
        ]
    )


class TestDELBarcodeSchema:
    """
    Test the DELBarcodeSchema class and its methods.

    Will also include all the test for the abstract BarcodeSchema class
    """

    def test_schema_construction(
        self,
        library_barcode_section,
        building_block_barcode_section_cycle1,
        building_block_barcode_section_cycle2,
        variable_barcode_section,
        static_barcode_section,
        umi_barcode_section,
    ):
        """Test DELBarcodeSchema construction and basic properties."""
        schema = DELBarcodeSchema(
            [
                library_barcode_section,
                static_barcode_section,
                building_block_barcode_section_cycle1,
                variable_barcode_section,
                building_block_barcode_section_cycle2,
                umi_barcode_section,
            ]
        )

        assert len(schema) == 25  # test __len__
        assert len(schema.barcode_sections) == 6
        assert len(schema.building_block_sections) == 2
        assert isinstance(schema.library_section, LibraryBarcodeSection)
        assert isinstance(schema.umi_section, UMIBarcodeSection)
        assert len(schema.static_sections) == 1

        # test __get_item__
        assert schema["library"] == library_barcode_section

        # test __eq__
        schema2 = DELBarcodeSchema(
            [
                library_barcode_section,
                static_barcode_section,
                building_block_barcode_section_cycle1,
                variable_barcode_section,
                building_block_barcode_section_cycle2,
                umi_barcode_section,
            ]
        )
        assert schema == schema2
        schema3 = DELBarcodeSchema(
            [
                library_barcode_section,
                static_barcode_section,
                building_block_barcode_section_cycle1,
                building_block_barcode_section_cycle2,
                umi_barcode_section,  # missing one building block
            ]
        )
        assert schema != schema3

    def test_get_required_sections(self, barcode_schema):
        """Test get_required_sections method."""
        required = barcode_schema.get_required_sections()
        assert len(required) == 4  # library, bb1, bb2, umi
        assert barcode_schema.library_section in required
        assert barcode_schema.building_block_sections[0] in required
        assert barcode_schema.building_block_sections[1] in required
        assert barcode_schema.umi_section in required

    def test_get_section(self, barcode_schema, library_barcode_section):
        """Test get_section method."""
        assert barcode_schema.get_section("library") == library_barcode_section
        with pytest.raises(KeyError):
            barcode_schema.get_section("nonexistent")

    def test_get_full_barcode(self, barcode_schema):
        """Test get_full_barcode method."""
        full = barcode_schema.get_full_barcode()
        expected = "AGCT" + "AGCT" + "NNNN" + "NNN" + "NNNN" + "NNNNNN"
        assert full == expected

    def test_get_section_spans(self, barcode_schema):
        """Test get_section_spans method."""
        spans = barcode_schema.get_section_spans()
        assert spans["library"] == slice(0, 4)
        assert spans["static1"] == slice(4, 8)
        assert spans["bb1"] == slice(8, 12)
        assert spans["var1"] == slice(12, 15)
        assert spans["bb2"] == slice(15, 19)
        assert spans["umi"] == slice(19, 25)

    def test_get_section_span(self, barcode_schema):
        """Test get_section_span method."""
        span = barcode_schema.get_section_span("library")
        assert span == slice(0, 4)
        span_excl = barcode_schema.get_section_span("library", exclude_overhangs=True)
        assert span_excl == slice(0, 4)  # no overhang

    def test_get_required_span(self, barcode_schema):
        """Test get_required_span method."""
        start, end = barcode_schema.get_required_span()
        assert start == 0
        assert end == 25  # full length since all are required

    def test_get_static_sections_before_library(self, barcode_schema):
        """Test get_static_sections_before_library method."""
        before = barcode_schema.get_static_sections_before_library()
        assert len(before) == 0  # static is after

    def test_get_static_sections_after_library(self, barcode_schema):
        """Test get_static_sections_after_library method."""
        after = barcode_schema.get_static_sections_after_library()
        assert len(after) == 1
        assert after[0].section_name == "static1"

    def test_get_length_before_section(self, barcode_schema):
        """Test get_length_before_section method."""
        assert barcode_schema.get_length_before_section("library") == 0
        assert barcode_schema.get_length_before_section("static1") == 4
        assert barcode_schema.get_length_before_section("bb1") == 8

    def test_get_length_after_section(self, barcode_schema):
        """Test get_length_after_section method."""
        assert barcode_schema.get_length_after_section("umi") == 0
        assert barcode_schema.get_length_after_section("bb2") == 6  # umi length
        assert barcode_schema.get_length_after_section("library") == 21  # total - 4

    def test_get_length_between_sections(self, barcode_schema):
        """Test get_length_between_sections method."""
        dist = barcode_schema.get_length_between_sections("library", "bb1")
        assert dist == 4  # static1 length
        dist_rev = barcode_schema.get_length_between_sections("bb1", "library", include_direction=True)
        assert dist_rev == -4

    def test_get_direction_of_sections(self, barcode_schema):
        """Test get_direction_of_sections method."""
        assert barcode_schema.get_direction_of_sections("library", "bb1") == 1
        assert barcode_schema.get_direction_of_sections("bb1", "library") == -1
        assert barcode_schema.get_direction_of_sections("library", "library") == 0

    def test_has_umi(self, barcode_schema):
        """Test has_umi method."""
        assert barcode_schema.has_umi()

    def test_get_section_length(self, barcode_schema):
        """Test get_section_length method."""
        assert barcode_schema.get_section_length("library") == 4
        assert barcode_schema.get_section_length("umi") == 6

    def test_is_schema_align_compatible(self, barcode_schema):
        """Test is_schema_align_compatible method."""
        # Same schema
        assert barcode_schema.is_schema_align_compatible(barcode_schema)
        # Different lengths
        other_sections = [
            LibraryBarcodeSection("library", "AGCT"),
            BuildingBlockBarcodeSection(1, "bb1", "NNNN"),
            UMIBarcodeSection("umi", "NNNNNN"),
        ]
        other_schema = DELBarcodeSchema(other_sections)
        assert not barcode_schema.is_schema_align_compatible(other_schema)
        # Different names
        other_sections2 = [
            LibraryBarcodeSection("lib", "AGCT"),
            StaticBarcodeSection("static1", "AGCT"),
            BuildingBlockBarcodeSection(1, "bb1", "NNNN"),
            VariableBarcodeSection("var1", "NNN"),
            BuildingBlockBarcodeSection(2, "bb2", "NNNN"),
            UMIBarcodeSection("umi", "NNNNNN"),
        ]
        other_schema2 = DELBarcodeSchema(other_sections2)
        assert not barcode_schema.is_schema_align_compatible(other_schema2)

    def test_is_static_library_locate_compatible(self, barcode_schema):
        """Test is_static_library_locate_compatible method."""
        # Same
        assert barcode_schema.is_static_library_locate_compatible(barcode_schema, ["static1"])
        # Different static tag
        other_sections = [
            LibraryBarcodeSection("library", "AGCT"),
            StaticBarcodeSection("static1", "TTTT"),  # different tag
            BuildingBlockBarcodeSection(1, "bb1", "NNNN"),
            VariableBarcodeSection("var1", "NNN"),
            BuildingBlockBarcodeSection(2, "bb2", "NNNN"),
            UMIBarcodeSection("umi", "NNNNNN"),
        ]
        other_schema = DELBarcodeSchema(other_sections)
        assert not barcode_schema.is_static_library_locate_compatible(other_schema, ["static1"])
        # New tag but not between library and static1
        other_sections2 = [
            LibraryBarcodeSection("library", "AGCT"),
            StaticBarcodeSection("static1", "AGCT"),
            StaticBarcodeSection("extra", "AA"),  # not between library and static1
            BuildingBlockBarcodeSection(1, "bb1", "NNNN"),
            VariableBarcodeSection("var1", "NNN"),
            BuildingBlockBarcodeSection(2, "bb2", "NNNN"),
            UMIBarcodeSection("umi", "NNNNNN"),
        ]
        other_schema2 = DELBarcodeSchema(other_sections2)
        assert barcode_schema.is_static_library_locate_compatible(other_schema2, ["static1"])
        # New tag between library and static1
        other_sections3 = [
            LibraryBarcodeSection("library", "AGCT"),
            StaticBarcodeSection("extra", "AA"),  # not between library and static1
            StaticBarcodeSection("static1", "AGCT"),
            BuildingBlockBarcodeSection(1, "bb1", "NNNN"),
            VariableBarcodeSection("var1", "NNN"),
            BuildingBlockBarcodeSection(2, "bb2", "NNNN"),
            UMIBarcodeSection("umi", "NNNNNN"),
        ]
        other_schema3 = DELBarcodeSchema(other_sections3)
        assert barcode_schema.is_static_library_locate_compatible(other_schema3, ["static1"])
        # missing required section
        other_sections4 = [
            LibraryBarcodeSection("library", "AGCT"),
            StaticBarcodeSection("static2", "AGCT"),
            StaticBarcodeSection("extra", "AA"),  # not between library and static1
            BuildingBlockBarcodeSection(1, "bb1", "NNNN"),
            VariableBarcodeSection("var1", "NNN"),
            BuildingBlockBarcodeSection(2, "bb2", "NNNN"),
            UMIBarcodeSection("umi", "NNNNNN"),
        ]
        other_schema4 = DELBarcodeSchema(other_sections4)
        with pytest.raises(KeyError):
            barcode_schema.is_static_library_locate_compatible(other_schema4, ["static1"])
        # not static section
        other_sections5 = [
            LibraryBarcodeSection("library", "AGCT"),
            MixedBarcodeSection("static1", "AGCT"),
            StaticBarcodeSection("extra", "AA"),  # not between library and static1
            BuildingBlockBarcodeSection(1, "bb1", "NNNN"),
            VariableBarcodeSection("var1", "NNN"),
            BuildingBlockBarcodeSection(2, "bb2", "NNNN"),
            UMIBarcodeSection("umi", "NNNNNN"),
        ]
        other_schema5 = DELBarcodeSchema(other_sections5)
        with pytest.raises(BarcodeSchemaError):
            barcode_schema.is_static_library_locate_compatible(other_schema5, ["static1"])

    def test_get_building_block_helpers(self, barcode_schema):
        """Test DELBarcodeSchema building-block helper methods."""
        assert barcode_schema.get_num_building_block_sections() == 2
        assert barcode_schema.get_building_block_section_names() == ["bb1", "bb2"]

    def test_get_length_before_section_with_section_object(self, barcode_schema):
        """Test that get_length_before_section accepts a BarcodeSection object as well as a name."""
        # using the first building block section object
        bb1_obj = barcode_schema.building_block_sections[0]
        assert barcode_schema.get_length_before_section(bb1_obj) == 8

    def test_get_required_span_with_include_sections(self, barcode_schema):
        """Test get_required_span when extra sections are explicitly included."""
        # include a static section (not required by default) and ensure no errors
        start, end = barcode_schema.get_required_span(include_sections=["static1"])
        assert start == 0
        assert end == 25


class TestToolCompoundBarcodeSchema:
    """Tests for ToolCompoundBarcodeSchema construction and behavior."""

    def test_tool_schema_construction_and_required_names(self):
        """Construct a ToolCompoundBarcodeSchema and verify attributes and required names."""
        sections = [
            LibraryBarcodeSection("library", "AGCT"),
            StaticBarcodeSection("prefix", "AA"),
            ToolCompoundRefBarcodeSection("tool_compound_ref", "TTTT"),
            BuildingBlockBarcodeSection(1, "bb1", "NNNN"),
            UMIBarcodeSection("umi", "NNNNNN"),
        ]
        schema = ToolCompoundBarcodeSchema(sections)
        assert isinstance(schema.tool_compound_ref_section, ToolCompoundRefBarcodeSection)
        required = schema.get_required_section_names()
        # required should include library and tool_compound_ref and umi
        assert "library" in required
        assert "tool_compound_ref" in required
        assert "umi" in required

    def test_tool_schema_missing_tool_compound_ref_raises(self):
        """Missing tool_compound_ref section should raise BarcodeSchemaError."""
        sections = [
            LibraryBarcodeSection("library", "AGCT"),
            BuildingBlockBarcodeSection(1, "bb1", "NNNN"),
            UMIBarcodeSection("umi", "NNNNNN"),
        ]
        with pytest.raises(BarcodeSchemaError):
            ToolCompoundBarcodeSchema(sections)

    def test_tool_schema_multiple_tool_compound_ref_raises(self):
        """Multiple tool_compound_ref sections should raise BarcodeSchemaError."""
        sections = [
            LibraryBarcodeSection("library", "AGCT"),
            ToolCompoundRefBarcodeSection("tool_compound_ref", "TTTT"),
            ToolCompoundRefBarcodeSection("tool_compound_ref_2", "CCCC"),
            BuildingBlockBarcodeSection(1, "bb1", "NNNN"),
        ]
        with pytest.raises(BarcodeSchemaError):
            ToolCompoundBarcodeSchema(sections)
