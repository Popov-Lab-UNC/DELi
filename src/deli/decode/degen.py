"""code for degenerating barcode reads based on UMI tags"""

import abc
import dataclasses
import os
from collections import defaultdict
from typing import Any, Generic, Optional, TypeVar, overload

from typing_extensions import Self

from deli.dels.library import Library
from deli.selection import Selection

from .decoder import DecodedCompound, DecodedDELCompound, DecodedToolCompound
from .umi import UMI


class CompoundUMICounter:
    """
    Keeps track of degenerate compound count for a single compound

    Notes
    -----
    if the UMI is "null" or `None`, it will be counted
    as a *unique* compound, since there is no way to
    determine if it is a duplicate or not

    Attributes
    ----------
    counter: defaultdict[str, int]
        the counter for the compound
    """

    def __init__(self):
        self.counter: dict[str, int] = {"null": 0, "None": 0}

    def __add__(self, other):
        """Add two CompoundUMICounters together by merging their counters"""
        if isinstance(other, CompoundUMICounter):
            new_counter = CompoundUMICounter()
            # merge counters
            for umi, count in self.counter.items():
                new_counter.counter[umi] = new_counter.counter.get(umi, 0) + count
            for umi, count in other.counter.items():
                new_counter.counter[umi] = new_counter.counter.get(umi, 0) + count
            return new_counter
        raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")

    @property
    def raw_count(self) -> int:
        """Get the raw count for the compound (number of times this compound was seen)"""
        return sum(lib_counter for lib_counter in self.counter.values())

    @property
    def degen_count(self) -> int:
        """Get the degenerate count for the compound (number of unique compounds seen)"""
        # exclude null and None from the unique count, but add their counts
        # since we are treating any compound missing a UMI as unique
        return len(self.counter) - 2 + self.counter.get("null", 0) + self.counter.get("None", 0)

    def add_umi(self, umi: str | UMI | None):
        """
        Add a UMI to the compound counter

        Parameters
        ----------
        umi: Any
            The UMI to add to the counter.
            Will cast to string

        Returns
        -------
        bool
            True if UMI was novel, else False
        """
        self.counter[str(umi)] += 1

    def to_json_dict(self) -> dict:
        """
        Convert the compound UMI counter to a JSON serializable dict

        Returns
        -------
        dict
        """
        return {key: val for key, val in self.counter.items() if val > 0}


C = TypeVar("C", bound=DecodedCompound)


class DegenCounter(abc.ABC, Generic[C]):
    """Base class for all degeneration counters"""

    @abc.abstractmethod
    def add_compound(self, compound: C):
        """
        Add a decoded compound to the degenerator

        Parameters
        ----------
        compound: DecodedCompound
            the decoded DEL compound to add
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def get_total_raw_count(self) -> int:
        """
        Get the total raw count

        Returns
        -------
        int
            the total raw count
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def get_total_degen_count(self) -> int:
        """
        Get the total degenerate count

        Returns
        -------
        int
            the total degenerate count
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def to_json_dict(self, *args, **kwargs) -> dict:
        """
        Convert the degenerator to a JSON serializable dict

        Returns
        -------
        dict
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def __add__(self, other: Self) -> Self:
        """
        Add two DegenCounters together by merging their counters

        Parameters
        ----------
        other: DegenCounter
            the other counter to add

        Returns
        -------
        DegenCounter
        """
        raise NotImplementedError()


class CompoundDegenCounter(DegenCounter[DecodedCompound]):
    """
    Counts degenerate compounds for any decoded compound

    Notes
    -----
    unlike `DELDegenCounter` and `ToolDegenCounter`, this will use the
    decoded compound directly as the key, rather than converting to a string
    representation. This is useful when you are degenerating in the same
    process as decoding, and you are not writing to a file after degeneration.
    This way, the information about the compound itself (building block SMILES, any
    enumerated SMILES, etc.) will be preserved. Otherwise, it needs to be reloaded
    later to recover that info.

    This is primarily used for the 'decode' CLI and python API when degeneration
    will occur in the same process as decoding.

    Attributes
    ----------
    counter: defaultdict[DecodedCompound, CompoundUMICounter]
        the counter for the compounds
    """

    def __init__(self):
        self.counter: defaultdict[DecodedCompound, CompoundUMICounter] = defaultdict(CompoundUMICounter)

    def __add__(self, other):
        """Add two CompoundDegenCounters together by merging their counters"""
        if isinstance(other, CompoundDegenCounter):
            new_counter = CompoundDegenCounter()
            # merge counters
            for compound, compound_counter in self.counter.items():
                new_counter.counter[compound] = (
                    new_counter.counter.get(compound, CompoundUMICounter()) + compound_counter
                )
            for compound, compound_counter in other.counter.items():
                new_counter.counter[compound] = (
                    new_counter.counter.get(compound, CompoundUMICounter()) + compound_counter
                )
            return new_counter
        raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")

    def add_compound(self, compound: DecodedCompound):
        """
        Add a decoded compound to the degenerator

        Will adjust the counts based on the compound's UMI

        Parameters
        ----------
        compound: DecodedCompound
        """
        self.counter[compound].add_umi(compound.umi)

    def get_total_raw_count(self) -> int:
        """Get the total raw count for the library"""
        return sum(compound_counter.raw_count for compound_counter in self.counter.values())

    def get_total_degen_count(self) -> int:
        """Get the total degenerate count for the library"""
        return sum(compound_counter.degen_count for compound_counter in self.counter.values())

    def to_json_dict(self, include_bb_smi: bool = False, enumerate_smi: bool = False, *args, **kwargs) -> dict:
        """
        Convert the degenerator to a JSON serializable dict

        Parameters
        ----------
        include_bb_smi: bool, false
            if True, will include building block SMILES in the JSON
        enumerate_smi: bool, false
            if True, will attempt to enumerate SMILES for compounds and include in the output

        Returns
        -------
        dict
        """
        data = {}
        for compound, compound_counter in self.counter.items():
            compound_data: dict[str, Any] = {"umi_counts": compound_counter.to_json_dict()}
            if enumerate_smi:
                compound_data["smi"] = compound.get_smiles()

            if isinstance(compound, DecodedToolCompound):
                _id = compound.tool_compound.compound_id
            elif isinstance(compound, DecodedDELCompound):
                _id = tuple(compound.building_block_ids)
                if include_bb_smi:
                    compound_data["bb_smis"] = [bb.smi for bb in compound.building_blocks]
            else:
                _id = compound.get_library_id()

            data[_id] = compound_data
        return data

    @classmethod
    def from_json_dict(cls, data: dict, library: Library) -> Self:
        """
        Load the degenerator from a JSON serializable dict *and* reconstruct the compound objects

        Parameters
        ----------
        data: dict
            the JSON serializable dict
        library: Library
            the library to use for reconstructing DEL compounds

        Returns
        -------
        CompoundDegenCounter
        """
        new_counter = cls()
        for compound_key, compound_info in data.items():
            try:
                compound = library.get_compound(compound_key)
            except Exception as e:
                raise KeyError(f"Cannot find compound with id {compound_key} in library {library.library_id}") from e

            new_counter.counter[compound].counter.update(compound_info["umi_counts"])
        return new_counter


class DELDegenCounter(DegenCounter[DecodedDELCompound]):
    """
    Counts degenerate DEL compounds for a single DEL

    Notes
    -----
    This assumes that all compounds added to the counter
    belong to the same DEL library. Therefore, the counter
    key is just the building block IDs.

    Attributes
    ----------
    counter: defaultdict[tuple[str, ...], CompoundUMICounter]
        the counter for the compounds
        the key is a tuple of building block IDs
    """

    def __init__(self):
        self.counter: defaultdict[tuple[str, ...], CompoundUMICounter] = defaultdict(CompoundUMICounter)

    def __add__(self, other):
        """Add two DELDegenCounter together by merging their counters"""
        if isinstance(other, DELDegenCounter):
            new_counter = DELDegenCounter()
            # merge counters
            for bb_id_tuple, compound_counter in self.counter.items():
                new_counter.counter[bb_id_tuple] = (
                    new_counter.counter.get(bb_id_tuple, CompoundUMICounter()) + compound_counter
                )
            for bb_id_tuple, compound_counter in other.counter.items():
                new_counter.counter[bb_id_tuple] = (
                    new_counter.counter.get(bb_id_tuple, CompoundUMICounter()) + compound_counter
                )
            return new_counter
        raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")

    def add_compound(self, compound: DecodedDELCompound | tuple[tuple[str, ...], str]):
        """
        Add a decoded DEL compound to the degenerator

        Will also accept the tuple of building block IDs directly
        attached to the UMI string

        Parameters
        ----------
        compound: DecodedDELCompound | tuple[tuple[str, ...], str]
            the decoded DEL compound to add
            can be the building block ID tuple directly with UMI
        """
        if isinstance(compound, tuple):
            self.counter[tuple(compound[0])].add_umi(compound[1])
        else:
            self.counter[tuple(compound.building_block_ids)].add_umi(compound.umi)

    def get_total_raw_count(self) -> int:
        """Get the total raw count for the library"""
        return sum(compound_counter.raw_count for compound_counter in self.counter.values())

    def get_total_degen_count(self) -> int:
        """Get the total degenerate count for the library"""
        return sum(compound_counter.degen_count for compound_counter in self.counter.values())

    def to_json_dict(self, *args, **kwargs) -> dict:
        """
        Convert the degenerator to a JSON serializable dict

        Returns
        -------
        dict
        """
        data = {}
        for compound, compound_counter in self.counter.items():
            data[compound] = {"umi_counts": compound_counter.to_json_dict()}
        return data

    @classmethod
    def from_json_dict(cls, data: dict) -> "DELDegenCounter":
        """
        Load the degenerator from a JSON serializable dict

        Notes
        -----
        Will ignore any non-count information in the JSON file.
        This means any info about SMILES will be lost when
        loaded and cannot be rewritten later

        Parameters
        ----------
        data: dict
            the JSON serializable dict

        Returns
        -------
        DELDegenCounter
        """
        new_counter = cls()
        for compound_key, compound_info in data.items():
            new_counter.counter[compound_key].counter.update(compound_info["umi_counts"])
        return new_counter


class ToolDegenCounter(DegenCounter[DecodedToolCompound]):
    """
    Counts degenerate tool compounds

    Notes
    -----
    Tool compounds are basically libraries with
    1 compound in them (that tool compounds ID)
    so degeneration is just counting unique
    occurrences of each tool compound UMI, rather
    than also tracking building block IDs

    Attributes
    ----------
    counter: CompoundUMICounter
        the counter for the tool compound
    """

    def __init__(self, tool_compound_id: str):
        self.counter = CompoundUMICounter()
        self.tool_compound_id = tool_compound_id

    def __add__(self, other):
        """Add two ToolDegenCounters together by merging their counters"""
        if isinstance(other, ToolDegenCounter):
            if self.tool_compound_id != other.tool_compound_id:
                raise ValueError("Cannot add ToolDegenCounters for different tool compounds")
            new_counter = ToolDegenCounter(self.tool_compound_id)
            new_counter.counter = self.counter + other.counter
            return new_counter
        raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")

    def add_compound(self, compound: DecodedToolCompound | str):
        """
        Add a decoded tool compound to the degenerator

        Can also accept the UMI string directly
        (since there is only one possible compound per in the degen so id not needed)

        Parameters
        ----------
        compound: DecodedToolCompound | str
            the decoded tool compound to add
            or the UMI string directly
        """
        if isinstance(compound, str):
            self.counter.add_umi(compound)
        else:
            self.counter.add_umi(compound.umi)

    def get_total_raw_count(self) -> int:
        """Get the total raw count for the tool compound"""
        return self.counter.raw_count

    def get_total_degen_count(self) -> int:
        """Get the total degenerate count for the tool compound"""
        return self.counter.degen_count

    def to_json_dict(self, *args, **kwargs) -> dict:
        """
        Convert the degenerator to a JSON serializable dict

        Returns
        -------
        dict
        """
        return {self.tool_compound_id: {"umi_counts": self.counter.to_json_dict()}}

    @classmethod
    def from_json_dict(cls, data: dict) -> "ToolDegenCounter":
        """
        Load the degenerator from a JSON serializable dict

        Notes
        -----
        Will ignore any non-count information in the JSON file.
        This means any info about SMILES will be lost when
        loaded and cannot be rewritten later

        Parameters
        ----------
        data: dict
            the JSON serializable dict

        Returns
        -------
        ToolDegenCounter
        """
        tool_compound_id, compound_info = next(iter(data.items()))
        new_counter = cls(tool_compound_id)
        new_counter.counter.counter.update(compound_info["umi_counts"])
        return new_counter


@dataclasses.dataclass(frozen=True)
class DegenSettings:
    """
    Settings for running degeneration

    Degeneration is the process of collapsing reads based on UMI tags
    to get a more accurate count of unique molecules, avoiding noise
    created by PCR amplification and sequencing bias.

    Parameters
    ----------
    umi_clustering: bool, default = False
        if True, will use UMI clustering for UMI degeneration
    ignore_missing_umi: bool, default = False
        if True, will ignore reads missing UMI tags during degeneration;
        they will not be counted at all
    """

    umi_clustering: bool = False
    ignore_missing_umi: bool = False

    def to_file(self, path: str):
        """
        Save settings to a YAML file

        Parameters
        ----------
        path : str
            path to save settings to
        """
        import yaml

        yaml.dump(dataclasses.asdict(self), open(path, "w"))

    @classmethod
    def from_file(cls, path: str) -> "DegenSettings":
        """
        Load settings from a YAML file

        Will first check if there is a "decode_settings" key
        and load settings from that sub dict.
        Otherwise, will load from the YAML file keys

        Parameters
        ----------
        path : str
            Path to YAML file

        Returns
        -------
        DecodingSettings

        Raises
        ------
        RuntimeError
            if valid decode settings cannot be loaded from the passed YAML file
        """
        import yaml

        _data = yaml.safe_load(open(path, "r"))
        if "degen_settings" not in _data:
            try:
                return cls(**yaml.safe_load(open(path, "r")))
            except Exception as e:
                raise RuntimeError(f"Failed to load degen settings from {path}") from e
        else:
            try:
                return cls(**_data["degen_settings"])
            except Exception as e:
                raise KeyError(f"missing degen_settings key in {path}") from e


class Degenerator(abc.ABC):
    """Base class for all degeneration runners"""

    counter: dict[str, DegenCounter]
    settings: DegenSettings

    @abc.abstractmethod
    def degen_decoded_compound(self, compound: DecodedCompound):
        """
        Degenerate a decoded compound

        Parameters
        ----------
        compound: DecodedCompound
            the decoded DEL compound to add
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def combine_degenerator(self, other: Self) -> Self:
        """
        Merge (add) another degenerator into this one

        Parameters
        ----------
        other: Degenerator
            the other degenerator to merge
        """
        raise NotImplementedError()

    def __add__(self, other: Self) -> Self:
        """Add two degenerators together by combining their counters"""
        return self.combine_degenerator(other)

    def write_json(
        self, path: os.PathLike, exclude_umi: bool = False, include_bb_smi: bool = False, enumerate_smi: bool = False
    ):
        """
        Write the degenerator to a JSON file

        Parameters
        ----------
        path: os.PathLike
            the path to write the JSON file to
            will overwrite if file exists
        include_bb_smi: bool, false
            if True, will include building block SMILES in the JSON
            Ignored for tool compounds.
        enumerate_smi: bool, false
            if True, will attempt to enumerate SMILES for decoded compounds
            and include in the output
        """
        import json

        data = {"settings": dataclasses.asdict(self.settings), "counts": {}}
        for lib_id, counter in self.counter.items():
            data["counts"][lib_id] = {
                "type": str(counter.__class__.__name__),
                "count_info": counter.to_json_dict(include_bb_smi=include_bb_smi, enumerate_smi=enumerate_smi),
            }
        with open(path, "w") as f:
            json.dump(data, f, indent=4)


class CompoundSelectionDegenerator(Degenerator):
    """
    Degenerate decoded barcodes based on UMI tags

    Unlike the SelectionDecoder, this is not embarrassingly parallel since degeneration
    requires tracking all observed UMIs for each DEL compound. This means it will not be able
    to operate on the fly and will need to keep everything in memory until all decoded
    compounds have been processed. If it stops before this, the degeneration results will
    have incorrect counts. Counts will never be lower than the actual degenerate counts, but
    they may be higher.

    This also means that if two Degenerator instances are run for the same
    selection (e.g. in parallel seeing different sets of decoded compounds), the results
    cannot just be simply be added together since the UMI sets will not be the same.
    `DELSelectionDegenerator` implement custom adding logic to do this merge correctly.
    They can also be saved and loaded from JSON files to enable merge separate processes.

    Notes
    -----
    This implementation of the degenerator will use the DecodedCompound objects
    directly as keys rather than converting to string representations. This is to
    preserve all information about the compounds during degeneration, rather than
    needing to reload them later. If you do not need this information (for example,
    if you are just writing to a file after degeneration), consider using

    Parameters
    ----------
    degen_settings: DegenSettings
        the settings for degeneration
        see `DegenSettings` for details

    Attributes
    ----------
    counter: dict[str, CompoundDegenCounter]
        the counter for each library ID
        this is a default dict, defaulting to CompoundDegenCounter
        you should not add keys directly, use `degen_decoded_compound` instead
    """

    def __init__(self, degen_settings: Optional[DegenSettings]):
        self.settings = degen_settings if degen_settings is not None else DegenSettings()
        self.counter: defaultdict[str, CompoundDegenCounter] = defaultdict(CompoundDegenCounter)

    def degen_decoded_compound(self, compound: DecodedCompound):
        """
        Degenerate a decoded compound

        Parameters
        ----------
        compound: DecodedCompound
            the decoded DEL compound to add
        """
        library_id = compound.get_library_id()
        self.counter[library_id].add_compound(compound)

    def combine_degenerator(self, other: Self) -> Self:
        """
        Merge (add) another degenerator into this one

        Parameters
        ----------
        other: Degenerator
            the other degenerator to merge

        Returns
        -------
        CompoundSelectionDegenerator
            the merged degenerator
        """
        new_degenerator = CompoundSelectionDegenerator(self.settings)
        new_degenerator.merge_degenerator(self)
        new_degenerator.merge_degenerator(other)
        return new_degenerator

    def merge_degenerator(self, other: Self):
        """
        Merge (add) another degenerator into this one in place

        Parameters
        ----------
        other: Degenerator
            the other degenerator to merge
        """
        for lib_id, counter in other.counter.items():
            self.counter[lib_id] = self.counter[lib_id] + counter


class SelectionDegenerator(Degenerator):
    """
    Degenerate decoded barcodes based on UMI tags

    Unlike the SelectionDecoder, this is not embarrassingly parallel since degeneration
    requires tracking all observed UMIs for each DEL compound. This means it will not be able
    to operate on the fly and will need to keep everything in memory until all decoded
    compounds have been processed. If it stops before this, the degeneration results will
    have incorrect counts. Counts will never be lower than the actual degenerate counts, but
    they may be higher.

    This also means that if two SelectionDegenerator instances are run for the same
    selection (e.g. in parallel seeing different sets of decoded compounds), the results
    cannot just be simply be added together since the UMI sets will not be the same.
    `DELSelectionDegenerator` implement custom adding logic to do this merge correctly.
    They can also be saved and loaded from JSON files to enable merge separate processes.

    Degenerators will group compounds by their library ID automatically.
    This will avoid excessive memory usage encoding the same library multiple times.

    The SelectionDegenerator is also keyed using the tuple of building block id
    strings rather than the full DEL compound ID. This is to avoid issues where
    extracting the building block ids and compound ideas from string representations
    can be error-prone. Also, the compound ID hold the library ID, which is already
    being used to group the compounds in the degenerator making this key redundant.

    Notes
    -----
    This implementation of the degenerator uses string representations of the
    decoded compounds as keys rather than the DecodedCompound objects directly.
    This means it can degenerate from decoded compound files. It will also
    convert any decoded compound it degenerates into a string representation.
    This means that some information about the compounds (building block SMILES,
    any enumerated SMILES, etc.) will be lost during degeneration unless the
    compounds are reloaded later. If you want to preserve this information during
    degeneration, consider using `CompoundSelectionDegenerator`.

    Parameters
    ----------
    degen_settings: DegenSettings
        the settings for degeneration
        see `DegenSettings` for details

    Attributes
    ----------
    counter: dict[str, DELDegenCounter | ToolDegenCounter]
        the counter for each library ID
        will be either a DELDegenCounter or ToolDegenCounter
        based on the compound type.
        This is a default dict, defaulting to CompoundDegenCounter.
        You should not add keys directly, use `degen_decoded_compound` instead.
    """

    def __init__(self, degen_settings: Optional[DegenSettings] = None):
        self.settings = degen_settings if degen_settings is not None else DegenSettings()
        self.counter: dict[str, DELDegenCounter | ToolDegenCounter] = {}

    def degen_decoded_compound(self, compound: DecodedDELCompound | DecodedToolCompound | tuple[tuple[str, ...], str]):
        """
        Degenerate a decoded compound

        Given any decoded compound, will determine its type

        Parameters
        ----------
        compound
        """
        # handle decoded compound object input
        if not isinstance(compound, tuple):
            lib_id = compound.get_library_id()
            if lib_id not in self.counter:
                if isinstance(compound, DecodedDELCompound):
                    self.counter[lib_id] = DELDegenCounter()
                elif isinstance(compound, DecodedToolCompound):
                    self.counter[lib_id] = ToolDegenCounter(lib_id)
                else:
                    raise RuntimeError("This error should be unreachable, please report a bug!")
            self.counter[lib_id].add_compound(compound)
        # handle string based input
        else:
            compound_id, umi = compound
            if len(compound_id) == 1:  # assume tool compound
                if compound_id[0] not in self.counter:
                    self.counter[compound_id[0]] = ToolDegenCounter(compound_id[0])
                self.counter[compound_id[0]].add_compound(umi)
            else:
                lib_id = compound_id[0]
                if lib_id not in self.counter:
                    self.counter[lib_id] = DELDegenCounter()
                self.counter[lib_id].add_compound((compound[1:], umi))

    def combine_degenerator(self, other: Self) -> Self:
        """
        Merge (add) another degenerator into this one

        Parameters
        ----------
        other: Degenerator
            the other degenerator to merge

        Returns
        -------
        SelectionDegenerator
            the merged degenerator
        """
        new_degenerator = SelectionDegenerator(self.settings)
        new_degenerator.merge_degenerator(self)
        new_degenerator.merge_degenerator(other)
        return new_degenerator

    def merge_degenerator(self, other: Self):
        """
        Merge (add) another degenerator into this one in place

        Parameters
        ----------
        other: Degenerator
            the other degenerator to merge
        """
        for lib_id, counter in other.counter.items():
            if lib_id in self.counter:
                self.counter[lib_id] = self.counter[lib_id] + counter
            else:
                self.counter[lib_id] = counter


@overload
def load_degenerator_from_json(path: os.PathLike, selection: Selection) -> CompoundSelectionDegenerator: ...


@overload
def load_degenerator_from_json(path: os.PathLike) -> SelectionDegenerator: ...


@overload
def load_degenerator_from_json(path: os.PathLike, selection: Optional[Selection] = None) -> Degenerator: ...


def load_degenerator_from_json(path: os.PathLike, selection: Optional[Selection] = None) -> Degenerator:
    """
    Load a SelectionDegenerator from a JSON file

    Because of the way the JSON is structured,

    Notes
    -----
    Will ignore non-count information in the JSON file.

    If selection is provided, will load a CompoundSelectionDegenerator
    that reconstructs the Compound objects as associates them
    with their library information (like SMILES).
    This can be expensive if there are many compounds to reconstruct.
    If you are only loading to merge and then write new degenerators,
    it is probably better to not provide a selection to avoid this overhead.

    Parameters
    ----------
    path: os.PathLike
        the path to read the JSON file from
    selection: Selection, optional
        the selection used to decode the compounds.
        If provided, will load a CompoundSelectionDegenerator
        that reconstructs the Compound objects as associates them
        with their library information (like SMILES).
        Otherwise, will load a SelectionDegenerator that only
        contains info on the compounds ids.
    """
    import json

    with open(path, "r") as f:
        data = json.load(f)

    settings = DegenSettings(**data["settings"])
    count_data = data["counts"]

    if selection is not None:
        loaded_counters: defaultdict[str, CompoundDegenCounter] = defaultdict(CompoundDegenCounter)
        for lib_id, lib_counter_info in count_data.items():
            lib_counter = CompoundDegenCounter()
            library = selection.library_collection.get_library(lib_id)
            loaded_counters[lib_id] = lib_counter.from_json_dict(lib_counter_info["count_info"], library)
        degenerator = CompoundSelectionDegenerator(degen_settings=settings)
        degenerator.counter = loaded_counters
        return degenerator

    else:
        loaded_counters: dict[str, DELDegenCounter | ToolDegenCounter] = dict()
        for lib_id, lib_counter_info in count_data.items():
            counter_type = lib_counter_info["type"]
            if counter_type == "DELDegenCounter":
                loaded_counters[lib_id] = DELDegenCounter.from_json_dict(lib_counter_info["count_info"])
            elif counter_type == "ToolDegenCounter":
                loaded_counters[lib_id] = ToolDegenCounter.from_json_dict(lib_counter_info["count_info"])
            else:
                raise RuntimeError(f"Unknown counter type {counter_type} in degenerator JSON file")
        degenerator = SelectionDegenerator(degen_settings=settings)
        degenerator.counter = loaded_counters
        return degenerator
