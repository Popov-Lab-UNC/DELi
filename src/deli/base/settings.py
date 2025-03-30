"""base classes for setting objects"""

from os import PathLike
from typing import Literal, Self

import yaml


class BaseSettings:
    """
    Base class for all setting objects

    Setting object are used to organize the possible
    settings that DELi CLI commands might use.

    DELi uses YAML syntax for all setting files and
    only supports safe loading of yaml files.
    This means all variables saved in settings
    must be primitive objects in python
    """

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def to_file(self, path: str | PathLike):
        """
        Save settings to a yaml file

        Parameters
        ----------
        path : str | PathLike
            path to save settings to
        """
        yaml.dump(self.__dict__, open(path, "w"))

    def get(self, key: str, default: object | None = None) -> object | None:
        """
        Get the setting with the matching name

        If setting of that names does not exist, will return `default`

        Parameters
        ----------
        key : str
            setting name to look up
        default : object | None
            the value to return if setting does not exist

        Returns
        -------
        object | None
            the setting value or the default value if setting does not exist
        """
        try:
            return self.__dict__[key]
        except KeyError:
            return default

    @classmethod
    def from_file(cls, path: str | PathLike) -> Self:
        """
        Load settings from a yaml file

        Parameters
        ----------
        path : str | PathLike
            Path to yaml file

        Returns
        -------
        Self
        """
        try:
            return cls(**yaml.safe_load(open(path, "r")))
        except Exception as e:
            raise RuntimeError(f"Failed to load settings from {path}") from e


class DecodingSettings(BaseSettings):
    """
    Define parameters for decoding experiments

    Only parameters relating the algorithm should be here
    Setting relating to IO should be handled outside this context
    (like in the click command)
    """

    def __init__(
        self,
        library_error_tolerance: float = 0.1,
        min_library_overlap: int | None = 10,
        alignment_algorithm: Literal["semi", "hybrid"] = "semi",
        bb_calling_approach: Literal["alignment", "bio"] = "alignment",
        revcomp: bool = False,
        read_type: Literal["single", "paired"] = "single",
        use_hamming: bool = True,
        track_statistics: bool = True,
    ):
        """
        Initialize the decoder settings

        Notes
        -----
        More details about the exact effect of these settings can
        be found in the "Decoding" docs

        Parameters
        ----------
        library_error_tolerance: float, default = 0.2
            the percent error to be tolerated in the library section
            this will be converted to number of errors based on tag size
            and down to the nearest whole number
            for example, a library with 14 nucleotides would tolerate
            1, 2, and 4 errors for an error tolerance of 0.1, 0.2 and 0.3 respectively
        min_library_overlap: int or None, default = 7
            the minimum number of nucleotides required to match
            the library tag
            This is because the demultiplexing will accept truncated matches
            at the front/back of the tag. For example a tag of AGCTGGTTC
            could match a read of GTTC if the min overlap was <=4
            If `None`, will default to the exact length of the tag, meaning
            the whole tag is expected.
            The recommended value is greater than 8, as the odds of a match this strong
            to be accidental are low
        alignment_algorithm: Literal["semi", "hybrid"], default = "semi"
            the algorithm to use for alignment
            only used if bb_calling_approach is "alignment"
        read_type: Literal["single", "paired"], default = "single"
            the type of read
            paired are for paired reads
            all other read types are single
        revcomp: bool, default = False
            If true, search the reverse compliment as well
        bb_calling_approach: Literal["alignment"], default = "alignment"
            the algorithm to use for bb_calling
            right now only "alignment" mode is supported
        use_hamming: bool, default = True
            enable (`True`) or disable (`False`) hamming decoding
            only used if a library specifies tags as hamming encoded
            Note: if hamming encoded libraries are given, and `use_hamming` is
            `False`, the hamming decoding will not occur, even though it is possible
        track_statistics: bool, default = True
            track over statistic during decoding
            statistics include, for example, number or seq that failed
            library demultiplexing, failed bb looked, which cycles failed etc.
            Full details on all statistics can be found in the "Decoding" docs
        """
        super().__init__(
            library_error_tolerance=library_error_tolerance,
            min_library_overlap=min_library_overlap,
            alignment_algorithm=alignment_algorithm,
            read_type=read_type,
            revcomp=revcomp,
            bb_calling_approach=bb_calling_approach,
            use_hamming=use_hamming,
            track_statistics=track_statistics,
        )
