"""define barcode functionality"""

import abc
import math
import re
from typing import Optional, Self

from deli.design.hamming import QuaternaryHammingDecoder


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
            if len(set(self.section_overhang) - {"A", "C", "G", "T"}) == 0:
                raise BarcodeSchemaError(
                    f"barcode section {self.section_name} overhang contains invalid nucleotides"
                )

    def get_dna_sequence(self) -> str:
        """
        Get the full DNA sequence of the barcode section

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
        section_name: str,
        section_tag: str,
        section_overhang: Optional[str] = None,
        hamming_decoder: Optional[QuaternaryHammingDecoder] = None,
    ):
        """
        Initialize BuildingBlockBarcodeSection

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
        hamming_decoder: Optional[QuaternaryHammingDecoder]
            if the section tag is encoded with a hamming code,
            the corresponding hamming decoder needed to decode the tags
        """
        super().__init__(section_name, section_tag, section_overhang)
        self.hamming_decoder = hamming_decoder


class UMIBarcodeSection(VariableBarcodeSection):
    """Class for UMI barcode sections"""

    pass


class StaticBarcodeSection(BarcodeSection):
    """Base class for all static barcode sections"""

    def _validate(self):
        """Validate that all nucleotides are known and valid"""
        super()._validate()
        if isinstance(self.section_overhang, str):
            if len(set(self.section_overhang) - {"A", "C", "G", "T"}) == 0:
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
            section
            for section in self.barcode_sections
            if isinstance(section, LibraryBarcodeSection)
        ]
        if len(_library_sections) == 0:
            raise BarcodeSchemaError("barcode schemas must contain a library barcode section")
        elif len(_library_sections) > 1:
            raise BarcodeSchemaError(
                "barcode schemas must contain only one library barcode section"
            )
        else:
            self.library_section = _library_sections[0]

        # check for building block sections
        self.building_block_sections: list[BuildingBlockBarcodeSection]
        _building_block_sections = [
            section
            for section in self.barcode_sections
            if isinstance(section, BuildingBlockBarcodeSection)
        ]
        if len(_building_block_sections) == 0:
            raise BarcodeSchemaError(
                "barcode schemas must contain at least one building block barcode section"
            )
        else:
            self.building_block_sections = _building_block_sections

        # check for umi section (optional)
        self.umi_section: Optional[UMIBarcodeSection] = None
        _umi_section = [
            section for section in self.barcode_sections if isinstance(section, UMIBarcodeSection)
        ]
        if len(_building_block_sections) > 0:
            raise BarcodeSchemaError(
                "barcode schemas must contain at most one umi barcode section"
            )
        else:
            if _umi_section:
                self.umi_section = _umi_section[0]

        # check for closing section (optional)
        self.closing_section: Optional[ClosingBarcodeSection] = None
        _closing_section = [
            section
            for section in self.barcode_sections
            if isinstance(section, ClosingBarcodeSection)
        ]
        if len(_building_block_sections) > 0:
            raise BarcodeSchemaError(
                "barcode schemas must contain at most one closing barcode section"
            )
        else:
            if _closing_section:
                self.closing_section = _closing_section[0]

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
            each key value pair is section name and its relvant data.
            the data is also a dictionary, with three (3) possible keys
            - 'tag': the section DNA tag (required)
            - 'overhang': the section overhang tag (optional)
            - 'hamming_decoder': the hamming decoder format to use
                                 if the section is hamming encoded (optional)

        Returns
        -------
        BarcodeSchema
        """
        _sections: list[BarcodeSection] = list()
        for section_name, section_info in data.items():
            if re.match(r"^library_id$", section_name):
                _sections.append(
                    LibraryBarcodeSection(
                        section_name=section_name,
                        section_tag=section_info["tag"],
                        section_overhang=section_info.get("overhang"),
                    )
                )
            elif re.match(r"^bb[1-9][0-9]*$", section_name):
                decoder = (
                    QuaternaryHammingDecoder.load(section_info["hamming_decoder"])
                    if section_info.get("hamming_decoder") is not None
                    else None
                )
                _sections.append(
                    BuildingBlockBarcodeSection(
                        section_name=section_name,
                        section_tag=section_info["tag"],
                        section_overhang=section_info.get("overhang"),
                        hamming_decoder=decoder,
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

    def __getitem__(self, item: str) -> BarcodeSection:
        """Given a barcode section name, return that BarcodeSection object"""
        return self._barcode_section_map[item]

    def to_regex_pattern(self, error_tolerance: float = 0.1) -> str:
        """
        generate the full barcode regex pattern for regex searching

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


# class BarcodeSectionTrait(enum.Enum):
#     """defines traits that barcode sections can have"""
#
#     REQUIRED: str = "REQUIRED"
#     STATIC: str = "STATIC"
#     VARIABLE: str = "VARIABLE"
#
#
# KNOWN_BARCODE_SECTIONS = {
#     r"^library_id$": [BarcodeSectionTrait.STATIC],
#     r"^bb[1-9][0-9]*$": [BarcodeSectionTrait.VARIABLE],
#     r"^closing$": [BarcodeSectionTrait.STATIC],
#     r"^umi$": [BarcodeSectionTrait.REQUIRED]
# }
#
#
# class BarcodeSchemaError(Exception):
#     """exception raised when a barcode schema is invalid"""
#
#     pass
#
#
# class BarcodeSchemaFileError(Exception):
#     """exception raised when a barcode schema is invalid"""
#
#     pass
#
#
# class BarcodeSection:
#     """object for an individual barcode section"""
#
#     def __init__(
#         self,
#         section_name: str,
#         section_tag: str,
#         overhang_tag: Optional[str] = None,
#         hamming_decoder: Optional[QuaternaryHammingDecoder] = None,
#     ):
#         """
#         Initialize the barcode section
#
#         Parameters
#         ----------
#         section_name: str
#             the name of the barcode section
#         section_tag: str
#             DNA bases of the barcode section
#         overhang_tag: Optional[str], default = None
#             the DNA bases of the overhang associated with this section
#             if there is one
#         hamming_decoder: Optional[QuaternaryHammingDecoder]
#             the decoder associated with this barcode section
#         """
#         self.section_name = section_name
#         self.section_tag = section_tag.upper()
#         self.overhang_tag = overhang_tag.upper() if overhang_tag else overhang_tag
#         self.decoder = hamming_decoder
#
#         self.traits: List[BarcodeSectionTrait] = []
#         for _known_section_name, _traits in KNOWN_BARCODE_SECTIONS.items():
#             if re.search(_known_section_name, section_name):
#                 self.traits = _traits
#                 break
#         if len(self.traits) == 0:
#             if
#
#
#
#         if re.match(self.section_name:
#
#         # check that section is recognizable to DELi
#         if not _section_is_valid(self.section_name):
#             raise BarcodeSchemaError(f"unrecognised barcode section: {self.section_name}")
#
#         # load traits
#         self.traits = _get_valid_section_traits(self.section_name)
#         self._tag_is_variable = "N" in self.section_tag
#
#         # check for trait conflict with passed bases
#         self._check_trait_validity()
#
#         # check that overhang is statis
#         if self.overhang_tag and "N" in self.overhang_tag:
#             raise BarcodeSchemaError(
#                 f"overhang_tag cannot contain variable nucleotides: {overhang_tag}"
#             )
#
#     def _assign_traits(self):
#
#
#
#     def __repr__(self):
#         """How the section name and DNA"""
#         return (
#             f"{self.section_name}: {self.section_tag}"
#             f"{self.overhang_tag if self.overhang_tag else ''}"
#         )
#
#     def __str__(self):
#         """Coverts the section to just the DNA part"""
#         return self.section_tag + (self.overhang_tag if self.overhang_tag else "")
#
#     def __len__(self, include_overhang: bool = True):
#         """Gets the full length of the barcode section including the overhang"""
#         if include_overhang and self.overhang_tag:
#             return len(self.overhang_tag) + len(self.section_tag)
#         else:
#             return len(self.section_tag)
#
#     def __eq__(self, other):
#         """Return true if two sections are the same or None"""
#         if isinstance(other, BarcodeSection):
#             if self.section_tag == other.section_tag:
#                 if self.overhang_tag == other.overhang_tag:
#                     if self.section_name == other.section_name:
#                         return True
#         return False
#
#     def _check_trait_validity(self):
#         """Helper func to make sure there are no trait conflicts"""
#         _has_static_trait = self.has_trait(BarcodeSectionTrait.STATIC)
#
#         if self._tag_is_variable and _has_static_trait:
#             raise BarcodeSchemaError(
#                 f"static barcode section {self.section_name} "
#                 f"has variable nucleotides: {self.section_tag}"
#             )
#
#         if not self._tag_is_variable and not _has_static_trait:
#             raise BarcodeSchemaError(
#                 f"variable barcode section {self.section_name} "
#                 f"appears static: {self.section_tag}"
#             )
#
#     def is_variable(self) -> bool:
#         """Returns true if the section have variable nucleotides"""
#         return self._tag_is_variable
#
#     def is_static(self) -> bool:
#         """Returns true if the section have NO variable nucleotides"""
#         return not self._tag_is_variable
#
#     def is_required(self) -> bool:
#         """Returns true if the section has required barcode"""
#         return self.has_trait(BarcodeSectionTrait.REQUIRED)
#
#     def has_trait(self, trait: BarcodeSectionTrait) -> bool:
#         """
#         Returns true if the barcode section has the specified trait
#
#         Parameters
#         ----------
#         trait: BarcodeSectionTrait
#             the barcode section trait to check for
#
#         Returns
#         -------
#         bool
#         """
#         return trait in self.traits
#
#     def has_overhang(self) -> bool:
#         """Returns true if the barcode section has overhang"""
#         return self.overhang_tag is not None
#
#     def get_dna_sequence(self) -> str:
#         """Get the full DNA sequence of the section"""
#         return str(self)
#
#     def to_match_pattern(
#         self,
#         wildcard: bool = True,
#         error_tolerance: int = 0,
#         include_overhang: bool = True,
#         trim: int = -1,
#     ) -> str:
#         """
#         Converts this section to a barcode matching regex pattern
#
#         Parameters
#         ----------
#         wildcard: bool, default= True
#             if true, will return a regex of just wild cards
#             for the full length (with overhang) of the section
#         error_tolerance: int, default=0
#             amount of error to tolerate in the pattern
#             only used if above 0 and the section has the static trait
#         include_overhang: bool, default=True
#             even if the section is variable, include the
#             explict DNA of the overhang
#             overridden by `wildcard`
#         trim: int, default=-1
#             if -1 match the full static region
#             if not, only match the first <trim> base pairs and treat the rest
#             like wildcards
#
#         Returns
#         -------
#         pattern: str
#         """
#         _pattern = ""
#
#         if wildcard:
#             return f".{{{len(self)}}}"
#
#         if self.is_static():
#             if trim != -1:
#                 _section_tag = self.section_tag[:trim]
#                 _leftover = max(0, len(self.section_tag) - trim)
#             else:
#                 _section_tag = self.section_tag
#                 _leftover = 0
#
#             if error_tolerance > 0:
#                 _pattern = f"(?:{_section_tag}){{e<={error_tolerance}}}"
#             else:
#                 _pattern = f"{_section_tag}"
#
#             if _leftover > 0:
#                 _pattern += f".{{{_leftover}}}"
#
#         else:
#             f".{{{len(self.section_tag)}}}"
#
#         if include_overhang:
#             _pattern += self.overhang_tag if self.overhang_tag is not None else ""
#         else:
#             _pattern += f".{{{len(self.overhang_tag)}}}" if self.overhang_tag is not None else ""
#
#         return _pattern
#
#     def has_hamming_decoder(self) -> bool:
#         """
#         Returns true if the section has hamming decoder (i.e. is hamming encoded)
#
#         Returns
#         -------
#         bool
#         """
#         return self.decoder is not None
#
#
# class BarcodeSchema(DeliDataLoadable):
#     """contains data and metadata about a barcode schema"""
#
#     def __init__(
#         self,
#         schema_id: str,
#         barcode_sections: List[BarcodeSection],
#         use_overhang_in_spans: bool = False,
#         override_required: bool = False,
#     ):
#         """
#         Initialize the barcode schema
#
#         Parameters
#         ----------
#         schema_id: str
#             name/id of the barcode schema
#         barcode_sections: Dict[str, Optional[str]]
#             Dictionary of barcode sections (see barcode schema docs for more info)
#         use_overhang_in_spans: bool, default=False
#             consider overhang regions in all sections as part of that sections barcode
#         override_required: bool, default=False
#
#         """
#         self.schema_id = schema_id
#         self.barcode_sections = OrderedDict(
#             {barcode_section.section_name: barcode_section for
#             barcode_section in barcode_sections}
#         )
#         self.override_required = override_required
#
#         self._use_overhang_in_spans = use_overhang_in_spans
#
#         self._check_barcode_sections()
#
#         self.full_barcode = "".join([str(val) for val in self.barcode_sections.values()])
#         self.barcode_spans = self._build_barcode_spans()
#
#         self.num_cycles = len(self.get_bb_regions())
#
#         if not self.override_required:
#             if self.num_cycles < 2:
#                 raise BarcodeSchemaError(
#                     f"number of bb cycles must be at least 2; found {self.num_cycles}"
#                 )
#
#     @classmethod
#     @accept_deli_data_name(sub_dir="barcodes", extension="json")
#     def load(cls, path: str) -> Self:
#         """
#         Load a barcode schema from the DELi data directory
#
#         Notes
#         -----
#         Call also just pass the name of the file
#         (with or without the file extension)
#         and as long as DELI_DATA_DIR is set
#
#
#         Parameters
#         ----------
#         path: str
#             path of the barcode to load
#
#         Returns
#         -------
#         BarcodeSchema
#         """
#         _cls = cls.load_from_json(path)
#         _cls.loaded_from = path
#         return _cls
#
#     @classmethod
#     def load_from_json(cls, file_path: str) -> Self:
#         """
#         load a schema from a json file
#
#         Parameters
#         ----------
#         file_path: str
#             path to json file with schema
#
#         Returns
#         -------
#         BarcodeSchema
#             the loaded schema
#         """
#         data = json.load(open(file_path))
#
#         # load the id
#         if "id" not in data.keys():
#             raise BarcodeSchemaFileError(f"schema file {file_path} missing 'id' key")
#         _id = str(data["id"])
#
#         # load the sections
#         if "sections" not in data.keys():
#             raise BarcodeSchemaFileError(f"schema file {file_path} missing 'sections' key")
#         if not isinstance(data["sections"], dict):
#             raise BarcodeSchemaFileError(f"'sections' should be a dictionary got a {type(data)}")
#
#         _hamming_encoded_sections: Dict[str, str] = dict()
#         if "hamming_encoded" in data.keys():
#             _hamming_encoded_sections = data["hamming_encoded"]
#
#         _sections = []
#         _parsed_sections = list(data["sections"].items())
#         while _parsed_sections:
#             section_name, section_dna_tag = _parsed_sections.pop(0)
#
#             # check for overhang section with no parent section
#             if section_name.endswith("_overhang"):
#                 raise BarcodeSchemaFileError(
#                     f"found {section_name} but no " f"{section_name[:-9]} section before it"
#                 )
#
#             # check for overhang
#             _overhang_section_tag = None
#             if len(_parsed_sections) > 0:
#                 _next_section_name, _ = _parsed_sections[0]  # dont pop yet, just peak
#                 if (
#                     _next_section_name.endswith("_overhang")
#                     and _next_section_name[:-9] == section_name
#                 ):
#                     _, _overhang_section_tag = _parsed_sections.pop(0)
#
#             _decoder: Optional[QuaternaryHammingDecoder] = None
#             _hamming_encoded = section_name in _hamming_encoded_sections
#             if section_name in _hamming_encoded_sections.keys():
#                 _dec_name = _hamming_encoded_sections[section_name]
#
#                 if _dec_name == "default":
#                     _decoder = QuaternaryHammingDecoder(list(range(len(section_dna_tag))), True)
#                 else:
#                     _decoder = QuaternaryHammingDecoder.load(
#                         _hamming_encoded_sections[section_name]
#                     )
#
#             _sections.append(
#                 BarcodeSection(
#                     section_name=section_name,
#                     section_tag=section_dna_tag,
#                     overhang_tag=_overhang_section_tag,
#                     hamming_decoder=_decoder,
#                 )
#             )
#
#         return cls(_id, _sections)
#
#     def __getitem__(self, item: str) -> BarcodeSection:
#         """Return the barcode section; None if not present"""
#         return self.barcode_sections[item]
#
#     def __eq__(self, other):
#         """Return true if all sections AND the order of them are equal"""
#         if isinstance(other, BarcodeSchema):
#             if self.barcode_sections == other.barcode_sections:
#                 return True
#         return False
#
#     def _check_barcode_sections(self) -> None:
#         """Check barcode sections for errors and standardize text"""
#         _seen_sections = set()
#         _bb_section_order = list()
#
#         # loop through all passed barcode section and check them all
#         _seen_bb_section: bool = False
#         for barcode_section_name, _ in self.barcode_sections.items():
#             if re.match(r"^bb[1-9][0-9]*$", barcode_section_name):
#                 _seen_bb_section = True
#                 _bb_section_order.append(int(barcode_section_name[2:]))  # append only the int
#             if BarcodeSectionTrait.BEFORE_BB in _get_valid_section_traits(barcode_section_name):
#                 if _seen_bb_section:
#                     raise BarcodeSchemaError(
#                         f"section {barcode_section_name} must come before all BB sections"
#                     )
#
#             if barcode_section_name not in _seen_sections:
#                 _seen_sections.add(barcode_section_name)
#             else:
#                 raise BarcodeSchemaError(
#                     f"sections must only be defined once; "
#                     f"found two definitions for section {barcode_section_name}"
#                 )
#
#         # check that the order of the BB integers starts with `1`
#         # and in ascending and full (no missing integers)
#         if _bb_section_order != list(range(1, len(_bb_section_order) + 1)):
#             raise BarcodeSchemaError(
#                 f"the order of BB sections must start with 'bb1' and then "
#                 f"continue to acsend in order; e.g. next is `bb2`, `bb3`... "
#                 f"saw order '{['bb'+str(_) for _ in _bb_section_order]}', "
#                 f"expected order '{list(range(1, len(_bb_section_order)+1))}';"
#                 f"see 'Defining Barcode Schema' docs for more details"
#             )
#
#         # check that all required sections are present
#         _required_sections = set(
#             [
#                 key
#                 for key, traits in VALID_BARCODE_SECTIONS.items()
#                 if BarcodeSectionTrait.REQUIRED in traits
#             ]
#         )
#
#         if not self.override_required:
#             for _required_section in _required_sections:
#                 if not any(
#                     [
#                         re.search(_required_section, _seen_section)
#                         for _seen_section in _seen_sections
#                     ]
#                 ):
#                     raise BarcodeSchemaError(
#                         f"Missing required barcode section: {_required_section}"
#                     )
#
#     def _build_barcode_spans(self) -> Dict[str, Tuple[int, int]]:
#         """Given a barcode section dictionary, build a spans dictionary"""
#         _barcode_spans = {}
#         prev = 0
#         for key, val in self.barcode_sections.items():
#             _length = len(val)
#             if self._use_overhang_in_spans:
#                 _section_length = _length
#             else:
#                 _section_length = len(val.section_tag)
#             _barcode_spans[key] = (prev, _section_length + prev)
#             prev += _length
#         return _barcode_spans
#
#     def has_index(self) -> bool:
#         """
#         Return true if schema include and index region
#
#         Returns
#         -------
#         bool
#         """
#         return self.barcode_sections.get("index") is not None
#
#     def has_library(self) -> bool:
#         """
#         Return true if schema include and library region
#
#         Returns
#         -------
#         bool
#         """
#         return self.barcode_sections.get("library") is not None
#
#     def get_bb_regions(self) -> List[str]:
#         """
#         Return all the region names that are for bb_tags
#
#         Returns
#         -------
#         List[str]
#         """
#         _bb_sections = []
#         for barcode_section in self.barcode_sections.keys():
#             if re.search(r"^bb[1-9][0-9]*$", barcode_section):
#                 _bb_sections.append(barcode_section)
#         return _bb_sections
#
#     def get_all_section_names(self) -> List[str]:
#         """
#         Return all the barcode section names
#
#         Returns
#         -------
#         List[str]
#         """
#         return list(self.barcode_sections.keys())
#
#     def get_position(self, section_name: str) -> int:
#         """
#         Gets the position of the given section in the barcode (starts at 0)
#
#         Notes
#         -----
#         Will return `-1` if the section is not in the barcode
#
#         Parameters
#         ----------
#         section_name: str
#             name of the section to query
#
#         Returns
#         -------
#         int
#         """
#         if section_name in self.barcode_sections.keys():
#             return list(self.barcode_sections.keys()).index(section_name)
#         else:
#             return -1
#
#     def is_library_compatible(self, other: Self) -> bool:
#         """
#         Return True if the schema is library calling compatible with the other
#
#         Notes
#         -----
#         Library calling compatible ("library compatible" for short)
#         means two schema share the same index, primer and library regions
#
#         This is all that is need in order to call the library.
#
#         See "Advanced Barcoding" for more details
#
#         Parameters
#         ----------
#         other: BarcodeSchema
#             BarcodeSchema to compare to
#
#         Returns
#         -------
#         bool
#         """
#         return (
#             (self.barcode_sections.get("index") == other.barcode_sections.get("index"))
#             and (self.barcode_sections.get("primer") == other.barcode_sections.get("primer"))
#             and (self.barcode_sections.get("library") == other.barcode_sections.get("library"))
#             and (self.get_position("index") == other.get_position("index"))
#             and (self.get_position("library") == other.get_position("library"))
#             and (self.get_position("primer") == other.get_position("primer"))
#         )
#
#     def is_experiment_compatible(self, other: Self) -> bool:
#         """
#         Return True if the schema is can be in the same experiment as the other
#
#         Notes
#         -----
#         Experiment compatible means two schemas have the
#         same primer bases (same sequence) and the
#         primer is in the same position as the other
#         in the barcode (if primer is position 1,
#         then the other should have this be true as well
#
#         See "Advanced Barcoding" for more details
#
#         Parameters
#         ----------
#         other: BarcodeSchema
#             BarcodeSchema to compare to
#
#         Returns
#         -------
#         bool
#         """
#         return (self.barcode_sections.get("index") == other.barcode_sections.get("index")) and (
#             self.get_position("index") == other.get_position("index")
#         )
#
#     def to_csv_header(self) -> List[str]:
#         """
#         Returns the header sections post calling required for calls from this schema
#
#         Returns
#         -------
#         List[str]
#         """
#         _headers = ["DEL_ID", "umi"]
#
#         if self.has_index():
#             _headers.append("index")
#         if self.has_library():
#             _headers.append("library")
#         for bb_section in self.get_bb_regions():
#             _headers.append(bb_section)
#
#         return _headers
#
#     def is_hamming_encoded(self) -> bool:
#         """
#         Return True if the schema has 1 or more hamming encoded barcode section
#
#         Returns
#         -------
#         bool
#         """
#         return any([section.has_hamming_decoder() for section in self.barcode_sections.values()])
