"""Dataclasses for various settings contexts in DELi"""

import json
import os
import warnings
from dataclasses import MISSING, asdict, dataclass
from typing import Literal


@dataclass(frozen=True)
class DecodingSettings:
    """
    Define parameters for decoding experiments

    More details about the exact effect of these settings can
    be found in the "Decoding" docs

    Notes
    -----
    Only parameters relating to the algorithm should be here
    Setting relating to IO should be handled outside this context

    Parameters
    ----------
    ignore_tool_compounds: bool, default = False
        if true, will ignore any tool compounds during decoding
    demultiplexer_algorithm: Literal["cutadapt", "regex", "full"], default = "regex"
        The demultiplexing algorithm to use.
        - "cutadapt": use a cutadapt to locate sections
        - "regex": use a regular expression based demultiplexer
        - "full": use a full alignment based demultiplexer
    demultiplexer_mode: Literal["library", "single", "flanking"], default = "flanking"
        The demultiplexing section strategy to use.
        - "library": demultiplex by matching just the library tag
        - "single": demultiplex by matching a single static barcode section
        - "flanking": demultiplex by matching barcode sections that flank the library tag
        (flanking means one before and one after the tag)
    realign: bool, default = False
        if true, will perform a local realignment of the read to the
        libraries barcode schema *after* demultiplexing determine the library.
        This could help recover reads that have complex alignments due multiple indels
    library_error_tolerance: int, default = 1
        The number of errors you are willing to tolerate in any given barcode
        section during library demultiplexing. Will apply to each section
        independently. For example, a flanking demultiplexer will allow for
        1 error in *each* of the flanking sections.
    library_error_correction_mode_str: str, default = "levenshtein_dist:2,asymmetrical"
        The error correction mode string to use for library barcode
        calling during demultiplexing.
    min_library_overlap: int , default = 8
        if using a cutadapt style demultiplexer, this is the minimum number of bases
        that must align to the expected barcode section for a match to be called.
        See the cutadapt documentation for more details on this parameter.
    revcomp: bool, default = True
        If true, search the reverse compliment as well.
        In most cases it is faster to use an external tools
        to align and reverse compliment reads before decoding
    library_wiggle: bool, default = False
        if true, will extend library aligned sections by 1 bp on each sidd.
        Similar to wiggle, but used during library demultiplexing.
    wiggle: bool, default = False
        if true, will extend aligned sections by 1 bp on each side during decoding.
        Can help recover reads with INDELs if not using realignment after demultiplexing
        (and is much faster).
    decode_matching_approach: Literal["greedy", "first_best", "first_perfect", "search_all"] = "first_perfect",
        the approach to use when matching decoded sections to the library during decoding.
        - "greedy": find the best scoring match across all possible windows and return it, even if imperfect
        - "first_best": search all possible windows and return the first occurring best match, even if more
        than one window has the same best score
        - "first_perfect": search all possible windows and return the first perfect match found, and otherwise
        fail if more than one window (with non-perfect score) share the same best score
        - "search_all": search all possible windows and return the best scoring one. Unlike 'first-perfect', if
        multiple windows have a perfect score the read will be rejected as ambiguous rather than just taking
        the first perfect match found.
    max_read_length: int or None, default = None
        maximum length of a read to be considered for decoding
        if above the max, decoding will fail
        if `None` will default to 5x the min_read_length
    min_read_length: int or None, default = None
        minimum length of a read to be considered for decoding
        if below the min, decoding will fail
        if `None` will default to the smallest min match length of
        any library in the collection considered for decoding
        with 10bp of buffer
    default_error_correction_mode_str: str, default = "levenshtein_dist:1,asymmetrical"
        The default error correction mode string to use for decoding.
        If a barcode section lacks a specified error correction mode,
        this mode will be used.
        See the documentation for more details on the format of this string.
    """

    ignore_tool_compounds: bool = False
    demultiplexer_algorithm: Literal["cutadapt", "regex", "full"] = "regex"
    demultiplexer_mode: Literal["library", "single", "flanking"] = "single"
    realign: bool = False
    library_error_tolerance: int = 1
    library_error_correction_mode_str: str = "levenshtein_dist:2,asymmetrical"
    min_library_overlap: int = 8
    revcomp: bool = True
    library_wiggle: bool = False
    wiggle: bool = False
    decode_matching_approach: Literal["greedy", "first_best", "first_perfect", "search_all"] = "first_perfect"
    max_read_length: int | None = None
    min_read_length: int | None = None
    default_error_correction_mode_str: str = "levenshtein_dist:1,asymmetrical"

    def to_json(self, path: os.PathLike[str]) -> None:
        """
        Save settings to a JSON file

        Parameters
        ----------
        path : Pthlike[str]
            path to save settings to
        """
        json.dump(asdict(self), open(path, "w"))

    @classmethod
    def from_json(cls, path: os.PathLike[str]) -> "DecodingSettings":
        """
        Load settings from a JSON file

        Will first check if there is a "decode_settings" key
        and load settings from that sub dict.
        Otherwise, will load from the JSON file keys

        Parameters
        ----------
        path : Pathlike[str]
            Path to JSON file

        Returns
        -------
        DecodingSettings
        """
        return cls(json.load(open(path, "r")))

    @classmethod
    def compare_params(cls, other_params: dict) -> None:
        """
        Check if params provided are missing or include extra keys.

        Will only raise warnings if extra keys are found.
        If a key is missing, will raise a warning only if it has a default value.

        Parameters
        ----------
        **kwargs
            The parameters to check

        Raises
        ------
        ValueError
            If any invalid parameters are provided or if any required parameters are missing
        """
        for key in other_params:
            if key not in cls.__dataclass_fields__:
                warnings.warn(
                    f"Unexpected key '{key}' found in JSON file when loading DecodingSettings", stacklevel=1
                )
        # check no missing keys
        for field_name, field_info in cls.__dataclass_fields__.items():
            if field_name not in other_params and field_info.default is not MISSING: # only warn for keys with defaults
                warnings.warn(
                    f"Expected key '{field_name}' not found in JSON file when loading DecodingSettings", stacklevel=1
                )
