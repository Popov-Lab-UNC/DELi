"""define barcode functionality"""

import enum
import json
import re
from collections import OrderedDict
from typing import List, Optional, Self

from deli.configure import accept_deli_data

from .base import DeliDataLoadableMixin


class BarcodeSectionTrait(enum.Enum):
    """defines traits that barcode sections can have"""

    REQUIRED: str = "REQUIRED"
    STATIC: str = "STATIC"


VALID_BARCODE_SECTIONS = {
    r"^pre-index$": [BarcodeSectionTrait.STATIC],
    r"^index$": [],
    r"^primer$": [BarcodeSectionTrait.STATIC, BarcodeSectionTrait.REQUIRED],
    r"^library$": [],
    r"^bb[1-9][0-9]*$": [],
    r"^pre-umi$": [BarcodeSectionTrait.STATIC],
    r"^umi$": [BarcodeSectionTrait.REQUIRED],
    r"^closing_primer$": [BarcodeSectionTrait.STATIC],
}


def _section_is_valid(section_name: str) -> bool:
    """Check if a section is known to DELi"""
    for valid_section_name, _ in VALID_BARCODE_SECTIONS.items():
        if re.search(valid_section_name, section_name):
            return True
    return False


def _get_valid_section_traits(section_name: str) -> List[BarcodeSectionTrait]:
    """Get the traits from a barcode section"""
    for valid_section_name, valid_section_traits in VALID_BARCODE_SECTIONS.items():
        if re.search(valid_section_name, section_name):
            return valid_section_traits
    raise KeyError("cannot find barcode section '{}'".format(section_name))


class BarcodeSchemaError(Exception):
    """exception raised when a barcode schema is invalid"""

    pass


class BarcodeSchemaFileError(Exception):
    """exception raised when a barcode schema is invalid"""

    pass


class BarcodeSection:
    """object for an individual barcode section"""

    def __init__(self, section_name: str, section_tag: str, overhang_tag: Optional[str] = None):
        """
        Initialize the barcode section

        Parameters
        ----------
        section_name: str
            the name of the barcode section
            must match one of the valid barcode names
        section_tag: str
            DNA tag of the barcode section
        overhang_tag: Optional[str]
            the DNA tag of the overhang associated with this section
            if there is one
        """
        self.section_name = section_name
        self.section_tag = section_tag.upper()
        self.overhang_tag = overhang_tag.upper() if overhang_tag else overhang_tag

        # check that section is recognizable to DELi
        if not _section_is_valid(self.section_name):
            raise BarcodeSchemaError(f"unrecognised barcode section: {self.section_name}")

        # load traits
        self.traits = _get_valid_section_traits(self.section_name)
        self._tag_is_variable = "N" in self.section_tag

        # check for trait conflict with passed tag
        self._check_trait_validity()

        # check that overhang is statis
        if self.overhang_tag and "N" in self.overhang_tag:
            raise BarcodeSchemaError(
                f"overhang_tag cannot contain variable nucleotides: {overhang_tag}"
            )

    def __repr__(self):
        """How the section name and DNA"""
        return (
            f"{self.section_name}: {self.section_tag}"
            f"{self.overhang_tag if self.overhang_tag else ''}"
        )

    def __str__(self):
        """Coverts the section to just the DNA part"""
        return self.section_tag + (self.overhang_tag if self.overhang_tag else "")

    def __len__(self, include_overhang: bool = True):
        """Gets the full length of the barcode section including the overhang"""
        if include_overhang and self.overhang_tag:
            return len(self.overhang_tag) + len(self.section_tag)
        else:
            return len(self.section_tag)

    def __eq__(self, other):
        """Return true if two sections are the same"""
        if isinstance(other, BarcodeSection):
            if self.section_tag == other.section_tag:
                if self.overhang_tag == other.overhang_tag:
                    if self.section_name == other.section_name:
                        return True
        return False

    def _check_trait_validity(self):
        """Helper func to make sure there are no trait conflicts"""
        _has_static_trait = self.has_trait(BarcodeSectionTrait.STATIC)

        if self._tag_is_variable and _has_static_trait:
            raise BarcodeSchemaError(
                f"static barcode section {self.section_name} "
                f"has variable nucleotides: {self.section_tag}"
            )

        if not self._tag_is_variable and not _has_static_trait:
            raise BarcodeSchemaError(
                f"variable barcode section {self.section_name} "
                f"appears static: {self.section_tag}"
            )

    def is_variable(self) -> bool:
        """Returns true if the section have variable nucleotides"""
        return self._tag_is_variable

    def is_static(self) -> bool:
        """Returns true if the section have NO variable nucleotides"""
        return not self._tag_is_variable

    def is_required(self) -> bool:
        """Returns true if the section has required barcode"""
        return self.has_trait(BarcodeSectionTrait.REQUIRED)

    def has_trait(self, trait: BarcodeSectionTrait) -> bool:
        """
        Returns true if the barcode section has the specified trait

        Parameters
        ----------
        trait: BarcodeSectionTrait
            the barcode section trait to check for

        Returns
        -------
        bool
        """
        return trait in self.traits

    def has_overhang(self) -> bool:
        """Returns true if the barcode section has overhang"""
        return self.overhang_tag is not None

    def get_dna_sequence(self) -> str:
        """Get the full DNA sequence of the section"""
        return str(self)

    def to_match_pattern(
        self, wildcard: bool = True, error_tolerance: int = 0, include_overhang: bool = True
    ) -> str:
        """
        Converts this section to a barcode matching regex pattern

        Parameters
        ----------
        wildcard: bool, default= True
            if true, will return a regex of just wild cards
            for the full length (with overhang) of the section
        error_tolerance: int, default=0
            amount of error to tolerate in the pattern
            only used if above 0 and the section has the static trait
        include_overhang: bool, default=True
            even if the section is variable, include the
            explict DNA of the overhang
            overridden by `wildcard`

        Returns
        -------
        pattern: str
        """
        _pattern = ""

        if wildcard:
            return f".{{{len(self)}}}"

        if self.is_static():
            if error_tolerance > 0:
                _pattern = f"(?:{self.section_tag}){{e<={error_tolerance}}}"
            else:
                _pattern = f"{self.section_tag}"
        else:
            f".{{{len(self.section_tag)}}}"

        if include_overhang:
            _pattern += self.overhang_tag if self.overhang_tag is not None else ""
        else:
            _pattern += f".{{{len(self.overhang_tag)}}}" if self.overhang_tag is not None else ""

        return _pattern


class BarcodeSchema(DeliDataLoadableMixin):
    """contains data and metadata about a barcode schema"""

    def __init__(
        self,
        schema_id: str,
        barcode_sections: List[BarcodeSection],
        use_overhang_in_spans: bool = False,
    ):
        """
        Initialize the barcode schema

        Parameters
        ----------
        schema_id: str
            name/id of the barcode schema
        barcode_sections: Dict[str, Optional[str]]
            Dictionary of barcode sections (see barcode schema docs for more info)
        use_overhang_in_spans: bool, default=False
            consider overhang regions in all sections as part of that sections barcode
        """
        self.schema_id = schema_id
        self.barcode_sections = OrderedDict(
            {barcode_section.section_name: barcode_section for barcode_section in barcode_sections}
        )

        self._use_overhang_in_spans = use_overhang_in_spans

        self._check_barcode_sections()

        self.full_barcode = "".join([str(val) for val in self.barcode_sections.values()])
        self.barcode_spans = self._build_barcode_spans()

        self.num_cycles = len(self.get_bb_regions())

        if self.num_cycles < 2:
            raise BarcodeSchemaError(
                f"number of bb cycles must be at least 2; found {self.num_cycles}"
            )

    @classmethod
    @accept_deli_data(sub_dir="barcodes", extension="json")
    def load(cls, path: str) -> Self:
        """
        Load a barcode schema from the DELi data directory

        Notes
        -----
        This is decorated by `accept_deli_data`
        which makes this function actually take
          path_or_name: str
          deli_config: DeliConfig

        `path_or_name` can be the full path to the file
        or it can be the name of the object to load

        See `Storing DEL info` in docs for more details


        Parameters
        ----------
        path: str
            path of the barcode to load

        Returns
        -------
        BarcodeSchema
        """
        return cls.load_from_json(path)

    @classmethod
    def load_from_json(cls, file_path: str) -> Self:
        """
        load a schema from a json file

        Parameters
        ----------
        file_path: str
            path to json file with schema

        Returns
        -------
        BarcodeSchema
            the loaded schema
        """
        data = json.load(open(file_path))

        # load the id
        if "id" not in data.keys():
            raise BarcodeSchemaFileError(f"schema file {file_path} missing 'id' key")
        _id = str(data["id"])

        # load the sections
        if "sections" not in data.keys():
            raise BarcodeSchemaFileError(f"schema file {file_path} missing 'sections' key")
        if not isinstance(data["sections"], dict):
            raise BarcodeSchemaFileError(f"'sections' should be a dictionary got a {type(data)}")

        _sections = []
        _parsed_sections = list(data["sections"].items())
        while _parsed_sections:
            section_name, section_dna_tag = _parsed_sections.pop(0)

            # check for overhang section with no parent section
            if section_name.endswith("_overhang"):
                raise BarcodeSchemaFileError(
                    f"found {section_name} but no " f"{section_name[:-9]} section before it"
                )

            # check for overhang
            _overhang_section_tag = None
            if len(_parsed_sections) > 0:
                _next_section_name, _ = _parsed_sections[0]  # dont pop yet, just peak
                if (
                    _next_section_name.endswith("_overhang")
                    and _next_section_name[:-9] == section_name
                ):
                    _, _overhang_section_tag = _parsed_sections.pop(0)
            _sections.append(BarcodeSection(section_name, section_dna_tag, _overhang_section_tag))

        return cls(_id, _sections)

    def __getitem__(self, item: str) -> Optional[BarcodeSection]:
        """Return the barcode section; None if not present"""
        return self.barcode_sections.get(item)

    def __eq__(self, other):
        """Return true if all sections AND the order of them are equal"""
        if isinstance(other, BarcodeSchema):
            if self.schema_id == other.schema_id:
                if self.barcode_sections == other.barcode_sections:
                    return True
        return False

    def _check_barcode_sections(self):
        """Check barcode sections for errors and standardize text"""
        _seen_sections = set()

        # loop through all passed barcode section and check them all
        for barcode_section_name, _ in self.barcode_sections.items():
            if barcode_section_name not in _seen_sections:
                _seen_sections.add(barcode_section_name)
            else:
                raise BarcodeSchemaError(
                    f"sections must only be defined once; "
                    f"found two definitions for section {barcode_section_name}"
                )

        # check that all required sections are present

        _required_sections = set(
            [
                key
                for key, traits in VALID_BARCODE_SECTIONS.items()
                if BarcodeSectionTrait.REQUIRED in traits
            ]
        )

        for _required_section in _required_sections:
            if not any(
                [re.search(_required_section, _seen_section) for _seen_section in _seen_sections]
            ):
                raise BarcodeSchemaError(f"Missing required barcode section: {_required_section}")

    def _build_barcode_spans(self):
        """Given a barcode section dictionary, build a spans dictionary"""
        _barcode_spans = {}
        prev = 0
        for key, val in self.barcode_sections.items():
            _length = len(val)
            if self._use_overhang_in_spans:
                _section_length = _length
            else:
                _section_length = len(val.section_tag)
            _barcode_spans[key] = (prev, _section_length + prev)
            prev += _length
        return _barcode_spans

    def has_index(self) -> bool:
        """
        Return true if schema include and index region

        Returns
        -------
        bool
        """
        return self.barcode_sections.get("index") is not None

    def has_library(self) -> bool:
        """
        Return true if schema include and library region

        Returns
        -------
        bool
        """
        return self.barcode_sections.get("library") is not None

    def get_bb_regions(self) -> List[str]:
        """
        Return all the region names that are for bb_tags

        Returns
        -------
        List[str]
        """
        _bb_sections = []
        for barcode_section in self.barcode_sections.keys():
            if re.search(r"^bb[1-9][0-9]*$", barcode_section):
                _bb_sections.append(barcode_section)
        return _bb_sections

    def get_all_section_names(self) -> List[str]:
        """
        Return all the barcode section names

        Returns
        -------
        List[str]
        """
        return list(self.barcode_sections.keys())
