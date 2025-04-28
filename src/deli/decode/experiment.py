"""Decoding Experiment and setting loading and module"""

import random
import time
from os import PathLike
from typing import Any, Literal, Self

import yaml

from deli.dels import DELibrary, DELibraryPool, SequencedSelection


class DecodingExperimentParsingError(Exception):
    """Exception to raise when a decoding experiment file is invalid"""

    pass


class DecodingSettings(dict):
    """
    Define parameters for decoding experiments

    Only parameters relating to the algorithm should be here
    Setting relating to IO should be handled outside this context
    (like in the click command)
    """

    def __init__(
        self,
        library_error_tolerance: float = 0.1,
        min_library_overlap: int | None = 10,
        alignment_algorithm: Literal["semi", "hybrid"] = "semi",
        bb_calling_approach: Literal["alignment", "bio"] = "bio",
        revcomp: bool = False,
        max_read_length: int | None = None,
        min_read_length: int | None = None,
        read_type: Literal["single", "paired"] = "single",
        use_hamming: bool = True,
        umi_clustering: bool = False,
        umi_min_distance: int = 2,
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
        max_read_length: int or None, default = None
            maximum length of a read to be considered for decoding
            if above the max, decoding will fail
            if `None` will default to 5x the min_read_length
        min_read_length: int or None, default = None
            minimum length of a read to be considered for decoding
            if below the min, decoding will fail
            if `None` will default to smallest min match length of
            any library in the pool considered for decoding
        bb_calling_approach: Literal["alignment"], default = "alignment"
            the algorithm to use for bb_calling
            right now only "alignment" mode is supported
        use_hamming: bool, default = True
            enable (`True`) or disable (`False`) hamming decoding
            only used if a library specifies tags as hamming encoded
            Note: if hamming encoded libraries are given, and `use_hamming` is
            `False`, the hamming decoding will not occur, even though it is possible
        umi_clustering: bool, default = False
            when doing degeneration, consider two similar UMIs to be the same
            similarity is based on levenshtein distance and `umi_min_distance`
        umi_min_distance: int, default = 2
            the minimum distance between two UMIs to be considered unique
            only used `umi_clustering` is `True`
        """
        super().__init__(
            library_error_tolerance=library_error_tolerance,
            min_library_overlap=min_library_overlap,
            alignment_algorithm=alignment_algorithm,
            read_type=read_type,
            revcomp=revcomp,
            max_read_length=max_read_length,
            min_read_length=min_read_length,
            bb_calling_approach=bb_calling_approach,
            use_hamming=use_hamming,
            umi_clustering=umi_clustering,
            umi_min_distance=umi_min_distance,
        )

    def to_file(self, path: str | PathLike):
        """
        Save settings to a YAML file

        Parameters
        ----------
        path : str | PathLike
            path to save settings to
        """
        yaml.dump(self.__dict__, open(path, "w"))

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


class DecodingExperiment:
    """
    Defines a decoding experiment for a DEL selection

    Configure using the decoding settings
    """

    def __init__(
        self,
        library_pool: DELibraryPool,
        selections: list[SequencedSelection],
        decode_settings: DecodingSettings | None = None,
        experiment_id: str | None = None,
    ):
        """
        Initialize the experiment with the given settings

        Parameters
        ----------
        library_pool: DELibraryPool
            the library pool used in the section being decoded
        selections: list[SectionCondition]
            the selection conditions for the decoding experiment
            and their corresponding selection file information
        decode_settings: DecodingSettings
            Settings to use for decoding
        experiment_id: str or None, default = None
            the id of the decoding experiment
            if `None` will default to a random number based on
            timestamp the object was created

        """
        self.library_pool = library_pool
        self.selections = selections
        self.decode_settings: DecodingSettings = (
            decode_settings if decode_settings is not None else DecodingSettings()
        )
        self.experiment_id: str = (
            experiment_id
            if experiment_id
            else (str(time.time()).replace(".", "") + f"{random.randint(0, 1000000):06d}")
        )

    def to_file(self, out_path: str | PathLike):
        """
        Write experiment to a human-readable file

        Parameters
        ----------
        out_path: str or PathLike
            path to save experiment to
        """
        data = {
            "experiment_id": self.experiment_id,
            "libraries": [lib.loaded_from for lib in self.library_pool],
            "selections": [
                selection.to_dict(ignore_library_pool=True) for selection in self.selections
            ],
            "decode_settings": self.decode_settings.__dict__,
        }
        yaml.safe_dump(data, open(out_path, "w"))

    @classmethod
    def from_file(cls, file_path: str | PathLike) -> Self:
        """
        Load the experiment from a human-readable file

        Parameters
        ----------
        file_path: str or PathLike
            path to load experiment from

        Returns
        -------
        DecodingExperiment
        """
        data = yaml.safe_load(open(file_path, "r"))

        try:
            _libraries: list[str] = data["libraries"]
        except KeyError as e:
            raise DecodingExperimentParsingError(
                f"{file_path} decoding file does not contain a 'libraries' section"
            ) from e

        _library_pool = DELibraryPool([DELibrary.load(lib_path) for lib_path in _libraries])

        experiment_id = data.get("experiment_id", None)
        selections = [
            SequencedSelection.from_dict({**selection_data, **{"library_pool": _library_pool}})
            for selection_data in data.get("selections")
        ]

        _decode_settings: dict[str, Any] | None = data.get("decode_settings", None)
        if _decode_settings is None:
            _decode_setting_obj = DecodingSettings()
        else:
            try:
                _decode_setting_obj = DecodingSettings(**_decode_settings)
            except TypeError as e:
                _unknown_arg = e.args[0].split()[-1]
                raise DecodingExperimentParsingError(
                    f"unrecognized decoding settings: {_unknown_arg}"
                ) from e

        # check that all selection ids are unique
        selection_ids = [selection.selection_id for selection in selections]
        if len(set(selection_ids)) != len(selection_ids):
            from collections import Counter

            duplicate_ids = [key for key, count in Counter(selection_ids).items() if count > 1]
            raise DecodingExperimentParsingError(
                f"Decoding experiment contains selections with duplicate IDs: {duplicate_ids}\n"
                "ensure all selections have unique IDs"
            )

        return cls(
            experiment_id=experiment_id,
            selections=selections,
            library_pool=_library_pool,
            decode_settings=_decode_setting_obj,
        )
