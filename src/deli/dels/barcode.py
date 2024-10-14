"""define barcode functionality"""

import json
import re
from collections import OrderedDict
from typing import Dict, Optional, Self, Tuple

from deli.dels.configure import check_file_path


# id: (required, static)
VALID_BARCODE_SECTIONS = {
    "pre-index": (False, True),
    "index": (False, False),
    "primer": (True, True),
    "library_tag": (False, False),
    "library_overhang": (False, True),
    "bb1_tag": (True, False),
    "bb1_overhang": (False, True),
    "bb2_tag": (True, False),
    "bb2_overhang": (False, True),
    r"bb[1-9][0-9]*_tag": (False, False),
    r"bb[1-9][0-9]*_overhang": (False, True),
    "pre-umi": (False, True),
    "umi": (True, False),
    "closing_primer": (False, True),
}


def barcode_section_is_static(section: str) -> bool:
    """
    Return True if barcode section in static

    Parameters
    ----------
    section: str
        name of barcode section

    Returns
    -------
    bool
    """
    for key, (_, is_static) in VALID_BARCODE_SECTIONS.items():
        if re.match(key, section, re.IGNORECASE):
            return is_static
    raise ValueError(f"unrecognized barcode section: '{section}'")


def barcode_section_is_required(section: str) -> bool:
    """
    Return True if barcode section in required

    Parameters
    ----------
    section: str
        name of barcode section

    Returns
    -------
    bool
    """
    for key, (is_required, _) in VALID_BARCODE_SECTIONS.items():
        if re.match(key, section, re.IGNORECASE):
            return is_required
    raise ValueError(f"unrecognized barcode section: '{section}'")


class BarcodeSchemaError(Exception):
    """exception raised when a barcode schema is invalid"""

    pass


class BarcodeSchema:
    """
    contains data and metadata about a barcode schema
    """

    def __init__(
        self,
        schema_id: str,
        barcode_sections: Dict[str, Optional[str]],
        order: Optional[Dict[str, Optional[int]]] = None,
    ):
        """
        Initialize the barcode schema

        Parameters
        ----------
        schema_id: str
            name/id of the barcode schema
        barcode_sections: Dict[str, Optional[str]]
            Dictionary of barcode sections (see barcode schema docs for more info)
        order: Dict[str, Optional[int]]
            Dictionary of barcode section order (see barcode schema docs for more info)
        """
        self.schema_id = schema_id
        self.barcode_sections = barcode_sections
        self.order = order

        self._check_barcode_sections()
        self._sort_barcode_sections()

        self.full_barcode = self._generate_full_barcode()
        self.barcode_lengths = {
            key: len(val) for key, val in self.barcode_sections.items() if val is not None
        }
        self.barcode_spans = self._build_barcode_spans()

        self.num_cycles = sum(
            [
                re.match(r"bb[1-9][0-9]*_tag", section_name) is not None
                for section_name in self.barcode_sections.keys()
            ]
        )

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
        file_path = check_file_path(file_path, "barcodes")
        data = json.load(open(file_path))
        return cls(data["id"], data["sections"], data["order"])

    def has_index(self) -> bool:
        """
        Return true if schema include and index region

        Returns
        -------
        bool
        """
        return self.barcode_sections.get("index") is not None

    def has_library_tag(self) -> bool:
        """
        Return true if schema include and index region

        Returns
        -------
        bool
        """
        return self.barcode_sections.get("library_tag") is not None

    def get_del_tags(self) -> OrderedDict[str, Dict[str, Optional[str]]]:
        """
        get the DEL tags and their respective tag and overhang region

        Notes
        -----
        The tags will be sorted by order of appearance in the schema.

        will be keyed by name of region (ex: "library" or "bb1").
        Val of each key will itself be a dict with "tag" and "overhang" keys.
        Each will be the code encoded in that part (ex: "bb1_tag" or "bb1_overhang").
        If overhang is `null` will be None

        Returns
        -------
        Dict[str, Dict[str, Optional[str]]]
            See Notes
        """
        _return = OrderedDict()
        for key in self.barcode_sections.keys():
            if "_tag" in key:
                _tag_name = key.replace("_tag", "")
                _return[_tag_name] = {
                    "tag": self.barcode_sections[key],
                    "overhang": self.barcode_sections.get(_tag_name + "_overhang"),
                }
        return _return

    def get_del_tag_lengths(self) -> Dict[str, Tuple[int, int]]:
        """
        Get the lengths of each DEL tag region split by tag and overhang

        Notes
        -----
        The order of the List will be the same order the tag regions appear
        Dict keys are the DEL tag names (ex: "library" or "bb1")

        Returns
        -------
        Dict[str, Tuple[int, int]]
            tuple contains length of tag and length of overhang
            if overhang region was `null` will be length 0

        """
        _tags = self.get_del_tags()
        return {
            key: (
                0 if val["tag"] is None else len(val["tag"]),
                0 if val["overhang"] is None else len(val["overhang"]),
            )
            for key, val in _tags.items()
        }

    def __getitem__(self, item: str) -> Optional[str]:
        """Return the barcode section; None if not present"""
        return self.barcode_sections.get(item)

    def _check_barcode_sections(self):
        """Check barcode sections for errors and standardize text"""

        def _is_variable(seq: str) -> bool:
            return all([i == "N" for i in seq])

        _required_sections = set([key for key, (val, _) in VALID_BARCODE_SECTIONS.items() if val])
        _seen_sections = set()
        for key, val in self.barcode_sections.copy().items():
            if not any([re.match(valid_key, key) for valid_key in VALID_BARCODE_SECTIONS.keys()]):
                raise BarcodeSchemaError(f"Unrecognized barcode section: {key}")

            _barcode_section = val.upper()
            _variable = _is_variable(_barcode_section)
            if (not _variable) != barcode_section_is_static(key):
                raise BarcodeSchemaError(
                    f"section {key} should be "
                    f"{'static' if VALID_BARCODE_SECTIONS[key][1] else 'variable'} "
                    f"but found {'variable' if _variable else 'static'}"
                )

            self.barcode_sections[key] = _barcode_section  # DNA is upper case convention
            if key not in _seen_sections:
                _seen_sections.add(key)
            else:
                raise BarcodeSchemaError(
                    f"sections must only be defined once; "
                    f"found two definitions for section {key}"
                )
        if len(_required_sections - _seen_sections) > 0:
            raise BarcodeSchemaError(
                f"Missing required barcode sections: {_required_sections - _seen_sections}"
            )

        # pre-index cannot be present if index is not
        if self.barcode_sections.get("pre-index") is not None:
            if self.barcode_sections.get("index") is None:
                raise BarcodeSchemaError(
                    "pre-index cannot be defined if index is not set;"
                    " set index or use primer section instead"
                )

        # library_overhang cannot be present if library_tag is not:
        if self.barcode_sections.get("library_overhang") is not None:
            if self.barcode_sections.get("library_tag") is None:
                raise BarcodeSchemaError(
                    "library_overhang cannot be defined if library_tag is not set"
                )

        # a bb_overhang cannot be present if the bb_tag is not
        for key in self.barcode_sections.keys():
            if re.match(r"bb[1-9][0-9]*_overhang", key) is not None:
                _bb_tag = key.replace("overhang", "tag")
                if self.barcode_sections.get(_bb_tag) is None:
                    raise BarcodeSchemaError(
                        f"found {key} section but no corresponding {_bb_tag} section"
                    )

    def _sort_barcode_sections(self):
        """Sort barcode sections base on passed order"""
        if self.order is None:  # if no custom order just return barcode sections
            self.barcode_sections = OrderedDict(self.barcode_sections.items())
        order_dict = dict()
        for key, val in self.barcode_sections.items():
            _order = self.order.get(key)
            if _order is None:
                raise BarcodeSchemaError(f"Barcode section {key} not found in order")
            else:
                order_dict[_order] = (key, val)

        _ordered_barcode_sections = OrderedDict()
        for key in sorted(order_dict.keys()):
            section, tag = order_dict[key]
            _ordered_barcode_sections[section] = tag
        self.barcode_sections = _ordered_barcode_sections

    def _generate_full_barcode(self):
        """Given a barcode section dictionary, generate a full barcode"""
        _full_barcode = ""
        for _, dna_tag in self.barcode_sections.items():
            _full_barcode += dna_tag
        return _full_barcode

    def _build_barcode_spans(self):
        """Given a barcode section dictionary, build a spans dictionary"""
        _barcode_spans = {}
        prev = 0
        for key, val in self.barcode_sections.items():
            _length = len(val)
            _barcode_spans[key] = (prev, _length + prev)
            prev += _length
        return _barcode_spans
