"""define barcode functionality"""

import abc
import math
import re
from typing import Optional
from typing_extensions import Self


class BarcodeSchemaError(Exception):
    """exception raised when a barcode schema is invalid"""

    pass


class BarcodeSection(abc.ABC):
    """base class for all barcode sections"""

    def __init__(
        self, section_name: str, section_tag: str, section_overhang: Optional[str] = None
    ):
        """
        Initialize BarcodeSection

        Parameters
        ----------
        section_name: str
            barcode section name
        section_tag: str
            barcode section DNA tag
            use "N" for variable or regions
        section_overhang: Optional[str]
            DNA of overhang directly after section tag
            leave as `None` if no overhang
        """
        self.section_name = section_name
        self.section_tag = section_tag.upper()
        self.section_overhang = section_overhang.upper() if section_overhang else None

        self._validate()

    def _validate(self):
        """Validate the tha barcode section tag(s) are correct"""
        if isinstance(self.section_overhang, str):
            if len(set(self.section_overhang) - {"A", "C", "G", "T"}) != 0:
                raise BarcodeSchemaError(
                    f"barcode section {self.section_name} overhang contains invalid nucleotides"
                )

    def get_dna_sequence(self) -> str:
        """
        Get the full DNA observed_seq of the barcode section

        Will substitute "N" for any variable nucleotides
        """
        return self.section_tag + (self.section_overhang if self.section_overhang else "")

    def has_overhang(self) -> bool:
        """Returns true if the barcode section has overhang"""
        return self.section_overhang is not None

    def __len__(self) -> int:
        """Gets the full length of the barcode section including the overhang"""
        if self.section_overhang is not None:
            return len(self.section_overhang) + len(self.section_tag)
        else:
            return len(self.section_tag)

    def __eq__(self, other) -> bool:
        """Return `True` if two sections are the same"""
        if isinstance(other, BarcodeSection):
            if self.section_tag == other.section_tag:
                if self.section_overhang == other.section_overhang:
                    if self.section_name == other.section_name:
                        return True
        return False

    @abc.abstractmethod
    def to_regex_pattern(self, error_tolerance: int = 0) -> str:
        """
        Covert the barcode section in a regex pattern for searching

        Variable nucleotides will be wildcards

        Parameters
        ----------
        error_tolerance: int, default = 0
            the among of levenshtein edits (inserts, deletes and substitutions)
            that are allowable for a valid match

        Returns
        -------
        str: A regex match pattern
        """
        raise NotImplementedError()


class VariableBarcodeSection(BarcodeSection):
    """
    Base class for all barcode sections encoding variable regions

    A variable barcode section is one where the nucleotides are
    expected to vary read from read (for example, in building block sections).
    All nucleotides in the section tag must be "N".
    The overhang (if included) must still be static.
    """

    def _validate(self):
        """Makes sure section has correct nucleotides"""
        super()._validate()
        if "N" in self.section_tag:
            if len(set(self.section_tag)) > 1:
                raise BarcodeSchemaError(
                    f"barcode section {self.section_name} contains non-variable "
                    f"nucleotide: {self.section_tag}"
                )

    def to_regex_pattern(self, error_tolerance: int = 0) -> str:
        """
        Covert the barcode section in a regex pattern for searching

        Variable nucleotides will be wildcards

        Parameters
        ----------
        error_tolerance: int, default = 0
            the among of levenshtein edits (inserts, deletes and substitutions)
            that are allowable for a valid match

        Returns
        -------
        str: A regex match pattern
        """
        return (
            f".{{{len(self.section_tag)}}}" + self.section_overhang
            if self.section_overhang is not None
            else ""
        )


class BuildingBlockBarcodeSection(VariableBarcodeSection):
    """Base class for all barcode sections encoding building block regions"""

    def __init__(
        self,
        cycle_number: int,
        section_name: str,
        section_tag: str,
        section_overhang: Optional[str] = None,
        error_correction_mode: Optional[str] = None,
    ):
        """
        Initialize BuildingBlockBarcodeSection

        Parameters
        ----------
        cycle_number: int
            the cycle number for this building block
            the first reactant is always cycle 0 and so on
        section_name: str
            barcode section name
        section_tag: str
            barcode section DNA tag
            use "N" for variable or regions
        section_overhang: Optional[str], default=None
            DNA of overhang directly after section tag
            leave as `None` if no overhang
        error_correction_mode: Optional[str], default=None
            the error correction mode for the building block tags in this section
            See Error Correction docs for more information
        """
        super().__init__(section_name, section_tag, section_overhang)
        self.cycle_number = cycle_number
        self.error_correction_mode = error_correction_mode


class UMIBarcodeSection(VariableBarcodeSection):
    """Class for UMI barcode sections"""

    pass


class StaticBarcodeSection(BarcodeSection):
    """Base class for all static barcode sections"""

    def _validate(self):
        """Validate that all nucleotides are known and valid"""
        super()._validate()
        if isinstance(self.section_overhang, str):
            if len(set(self.section_overhang) - {"A", "C", "G", "T"}) != 0:
                raise BarcodeSchemaError(
                    f"barcode section {self.section_name} contains invalid "
                    f"nucleotides: {self.section_tag}"
                )

    def to_regex_pattern(self, error_tolerance: int = 0) -> str:
        """
        Covert the barcode section in a regex pattern for searching

        Parameters
        ----------
        error_tolerance: int, default = 0
            the among of levenshtein edits (inserts, deletes and substitutions)
            that are allowable for a valid match

        Returns
        -------
        str: A regex match pattern
        """
        if error_tolerance > 0:
            _pattern = f"(?:{self.get_dna_sequence()}){{e<={error_tolerance}}}"
        else:
            _pattern = f"{self.get_dna_sequence()}"
        return _pattern


class HeadpieceBarcodeSection(StaticBarcodeSection):
    """Class for headpiece barcode sections"""

    pass


class LibraryBarcodeSection(StaticBarcodeSection):
    """Class for library barcode sections"""

    pass


class ClosingBarcodeSection(StaticBarcodeSection):
    """Class for closing barcode sections"""

    pass


class BarcodeSchema:
    """
    Barcode schema class

    BarcodeSchemas hold an ordered list of BarcodeSections that make up
    the full barcode definition.
    Only LibraryBarcodeSection and BuildingBlockBarcodeSection (<=1) are required.
    While not required, only up to one (1) UMIBarcodeSection
    and up to one (1) ClosingBarcodeSection are allowed.
    However not including these sections could limit what actions DELi
    can carry out.

    ALl other sections are "custom" (not special to DELi) and are optional and
    allowed in any number. The only limitation is that all BarcodeSections have
    unique section names.

    Attributes
    ----------
    library_section: LibraryBarcodeSection
        the library barcode section
    building_block_sections: List[BuildingBlockBarcodeSection]
        the building block barcode sections
        order of appearance in barcode is preserved
    umi_section: Optional[UMIBarcodeSection]
        the UMI barcode section if there is one, otherwise `None`
    closing_section: Optional[ClosingBarcodeSection]
        the closing barcode section if there is one, otherwise `None`

    """

    def __init__(self, barcode_sections: list[BarcodeSection]):
        """
        Initialize BarcodeSchema

        Parameters
        ----------
        barcode_sections: list[BarcodeSection]
            list of barcode sections *in order* for the barcode schema
        """
        self.barcode_sections = barcode_sections

        # check for uniqueness of barcode section names
        _barcode_section_ids = [section.section_name for section in self.barcode_sections]
        if len(set(_barcode_section_ids)) != len(_barcode_section_ids):
            from collections import Counter

            name_counts = Counter(_barcode_section_ids)
            raise BarcodeSchemaError(
                f"all barcode section names in a barcode schema must be unique;"
                f"section names {[key for key, val in name_counts.items() if val > 1]} "
                f"are not unique"
            )

        # a barcode section name map for easier indexing
        self._barcode_section_map: dict[str, BarcodeSection] = {
            section.section_name: section for section in self.barcode_sections
        }

        # check for one library section
        self.library_section: LibraryBarcodeSection
        _library_sections = [
            (i, section)
            for i, section in enumerate(self.barcode_sections)
            if isinstance(section, LibraryBarcodeSection)
        ]
        if len(_library_sections) == 0:
            raise BarcodeSchemaError("barcode schemas must contain a library barcode section")
        elif len(_library_sections) > 1:
            raise BarcodeSchemaError(
                "barcode schemas must contain only one library barcode section"
            )
        else:
            self.library_section = _library_sections[0][1]
            self._library_section_idx = _library_sections[0][0]

        # check for building block sections
        self.building_block_sections: list[BuildingBlockBarcodeSection] = list()
        self._building_block_section_idxs: list[int] = list()
        _building_block_sections = [
            (i, section)
            for i, section in enumerate(self.barcode_sections)
            if isinstance(section, BuildingBlockBarcodeSection)
        ]
        if len(_building_block_sections) == 0:
            raise BarcodeSchemaError(
                "barcode schemas must contain at least one building block barcode section"
            )
        else:
            for i, (_bb_section_idx, bb_section) in enumerate(_building_block_sections):
                if (i + 1) != bb_section.cycle_number:
                    raise BarcodeSchemaError(
                        f"expected building block section to be for cycle {i + 1},"
                        f"but found cycle number {bb_section.cycle_number} for "
                        f"building block section {bb_section.section_name}"
                    )
                else:
                    self.building_block_sections.append(bb_section)
                    self._building_block_section_idxs.append(_bb_section_idx)

        # check for umi section (optional)
        self.umi_section: Optional[UMIBarcodeSection] = None
        self._umi_section_idx: int = -1
        _umi_section = [
            (i, section)
            for i, section in enumerate(self.barcode_sections)
            if isinstance(section, UMIBarcodeSection)
        ]
        if len(_umi_section) > 1:
            raise BarcodeSchemaError(
                "barcode schemas must contain at most one umi barcode section"
            )
        else:
            if _umi_section:
                self.umi_section = _umi_section[0][1]
                self._umi_section_idx = _umi_section[0][0]

        # check for closing section (optional)
        self.closing_section: Optional[ClosingBarcodeSection] = None
        _closing_section = [
            section
            for section in self.barcode_sections
            if isinstance(section, ClosingBarcodeSection)
        ]
        if len(_closing_section) > 1:
            raise BarcodeSchemaError(
                "barcode schemas must contain at most one closing barcode section"
            )
        else:
            if _closing_section:
                self.closing_section = _closing_section[0]

        # check that library section is not inbetween building block sections
        _found_library = False
        _found_bb_section = False
        self._library_in_front = True
        for section in self.barcode_sections:
            if _found_bb_section and not _found_library:
                _library_in_front = False
            if isinstance(section, LibraryBarcodeSection):
                _found_library = True
            elif isinstance(section, BuildingBlockBarcodeSection):
                if _found_library and _found_bb_section:
                    raise BarcodeSchemaError(
                        "barcode schemas must not have the library section "
                        "between the building block sections"
                    )
                _found_bb_section = True
                _found_library = False

        # get the minimum length needed to extract info from a barcode match
        self.min_length: int
        required_section_indexes = [self._library_section_idx] + self._building_block_section_idxs
        if self.umi_section:
            required_section_indexes.append(self._umi_section_idx)
        start = min(required_section_indexes)
        end = max(required_section_indexes)
        self.min_length = sum([len(section) for section in self.barcode_sections[start : end + 1]])

    def __eq__(self, other):
        """Return `True` if two barcode schemas are the same"""
        if isinstance(other, BarcodeSchema):
            if len(self.barcode_sections) != len(other.barcode_sections):
                return False
            for section in self.barcode_sections:
                if section not in other.barcode_sections:
                    return False
            for section in other.barcode_sections:
                if section not in self.barcode_sections:
                    return False
            return True
        return False

    def __len__(self) -> int:
        """Get the length of all sections of the barcode schema"""
        return sum([len(section) for section in self.barcode_sections])

    @classmethod
    def from_dict(cls, data: dict) -> Self:
        """
        Load a barcode schema from a dictionary defining

        Notes
        -----
        Primarily used to load a barcode schema from a the library json file
        Should not be used outside this context.

        Parameters
        ----------
        data: dict
            a dictionary containing the barcode schema
            each key value pair is section name and its relevant data.
            the data is also a dictionary, with three (3) possible keys
            - 'tag': the section DNA tag (required)
            - 'overhang': the section overhang tag (optional)

        Returns
        -------
        BarcodeSchema
        """
        _sections: list[BarcodeSection] = list()
        for section_name, section_info in data.items():
            if re.match(r"^library$", section_name):
                _sections.append(
                    LibraryBarcodeSection(
                        section_name=section_name,
                        section_tag=section_info["tag"],
                        section_overhang=section_info.get("overhang"),
                    )
                )
            elif re.match(r"^headpiece$", section_name):
                _sections.append(
                    HeadpieceBarcodeSection(
                        section_name=section_name,
                        section_tag=section_info["tag"],
                        section_overhang=section_info.get("overhang"),
                    )
                )
            elif re.match(r"^bb[1-9][0-9]*$", section_name):
                _cycle_number = int(section_name[2:])
                _sections.append(
                    BuildingBlockBarcodeSection(
                        cycle_number=_cycle_number,
                        section_name=section_name,
                        section_tag=section_info["tag"],
                        section_overhang=section_info.get("overhang", None),
                        error_correction_mode=section_info.get("error_correction", None),
                    )
                )
            elif re.match(r"^closing$", section_name):
                _sections.append(
                    ClosingBarcodeSection(
                        section_name=section_name,
                        section_tag=section_info["tag"],
                        section_overhang=section_info.get("overhang"),
                    )
                )
            elif re.match(r"^umi$", section_name):
                _sections.append(
                    UMIBarcodeSection(
                        section_name=section_name,
                        section_tag=section_info["tag"],
                        section_overhang=section_info.get("overhang"),
                    )
                )
            else:
                # a static barcode
                if len(set(section_info["tag"]) - {"A", "C", "G", "T"}) == 0:
                    _sections.append(
                        StaticBarcodeSection(
                            section_name=section_name,
                            section_tag=section_info["tag"],
                            section_overhang=section_info.get("overhang"),
                        )
                    )
                # if not static assume variable
                else:
                    _sections.append(
                        VariableBarcodeSection(
                            section_name=section_name,
                            section_tag=section_info["tag"],
                            section_overhang=section_info.get("overhang"),
                        )
                    )
        return cls(barcode_sections=_sections)

    def __getitem__(self, item: str) -> BarcodeSection:
        """Given a barcode section name, return that BarcodeSection object"""
        return self._barcode_section_map[item]

    def to_regex_pattern(self, error_tolerance: float = 0.1) -> str:
        """
        Generate the full barcode regex pattern for regex searching

        Parameters
        ----------
        error_tolerance: float, default=0.1
            the amount of error to tolerate in each
            barcode section as a percentage of the section size.

        Returns
        -------
        str: regex pattern
        """
        if error_tolerance < 0 or error_tolerance > 1:
            raise ValueError("error_tolerance must be between 0 and 1")

        regex_str: str = ""
        for section in self.barcode_sections:
            regex_str += section.to_regex_pattern(
                error_tolerance=math.ceil(len(section) * error_tolerance)
            )
        return regex_str

    def get_full_barcode(self) -> str:
        """
        Return the full barcode as a string

        Will include variable regions as 'N'

        Returns
        -------
        str
        """
        full_barcode: str = ""
        for section in self.barcode_sections:
            full_barcode += section.get_dna_sequence()
        return full_barcode

    def get_section_spans(self, exclude_overhangs: bool = False) -> dict[str, slice]:
        """
        Get the spans of each barcode section as a dict

        Parameters
        ----------
        exclude_overhangs: bool, default=False
            exclude the overhang region of a tag when
            calculating the spans

        Returns
        -------
        dict[str, slice]
            keys are section ids and values are spans (as a `slice`)
        """
        spans: dict[str, slice] = dict()
        curr_pos: int = 0
        for section in self.barcode_sections:
            _section_len = len(section) - (
                (len(section.section_overhang) if section.section_overhang else 0)
                * exclude_overhangs
            )
            spans[section.section_name] = slice(curr_pos, curr_pos + _section_len)
            curr_pos += len(section)
        return spans

    def get_num_building_block_sections(self) -> int:
        """Get number of building block sections"""
        return len(self.building_block_sections)

    def is_library_tag_in_front(self) -> bool:
        """
        Return `True` if the library tag is at the front of the barcode

        Notes
        -----
        Because library tags must be before or after all barcode sections
        if this is `False`  it means the library tag is at the 'back' of the barcode

        Returns
        -------
        bool
        """
        return self._library_in_front

    def is_library_tag_in_back(self) -> bool:
        """
        Return `True` if the library tag is at the back of the barcode

        Notes
        -----
        Because library tags must be before or after all barcode sections
        if this is `False`  it means the library tag is at the 'front' of the barcode

        Returns
        -------
        bool
        """
        return not self._library_in_front

    def get_length_before_library(self) -> int:
        """
        Get the number of base pairs before the start of the library tag section

        Notes
        -----
        'Before' always referees to the 'left' of the section
        which is towards the 5' end of the barcode

        Returns
        -------
        int
        """
        return sum(
            [len(section) for section in self.barcode_sections[: self._library_section_idx]]
        )

    def get_length_after_library(self) -> int:
        """
        Get the number of base pairs after the end of the library tag section

        Notes
        -----
        'After' always referees to the 'right' of the section
        which is towards the 3' end of the barcode

        Returns
        -------
        int
        """
        return sum(
            [len(section) for section in self.barcode_sections[self._library_section_idx + 1 :]]
        )

    def has_umi(self) -> bool:
        """True is library has UMI barcode section, else False"""
        return self.umi_section is not None

    def get_section_length(self, section_name: str) -> int:
        """
        Get the length of a specific barcode section

        Parameters
        ----------
        section_name: str
            the name of the barcode section

        Returns
        -------
        int: length of the barcode section
        """
        if section_name not in self._barcode_section_map:
            raise KeyError(f"Barcode section '{section_name}' not found in schema")
        return len(self._barcode_section_map[section_name])
