"""Handle tool compounds for DELi"""

from typing import Literal, Sequence, overload

from deli.configure import DeliDataLoadable, accept_deli_data_name
from deli.utils.mol_utils import SmilesMixin

from .barcode import BarcodedMixin, ToolCompoundBarcodeSchema
from .library import Library
from .compound import Compound


class ToolCompound(Compound, SmilesMixin):
    """
    Contains information about a "tool" compound outside of DEL context.

    In DEL, tool compounds are used often added in as controls/references
    for quality control or experimental purposes. These compounds are not
    part of the combinatorial library, but are still important to track.
    This class provides a way to represent these compounds for DELi.

    As a results, ToolCompound inherits from Compound, but does not
    have any DEL-specific information, like library or building blocks.
    This means it only needs a user provided compound_id and optionally
    a SMILES string if its structure is known.

    See the decoding docs for more info on how to set up tool compounds
    in a decoding context. They are currently not utilized in the
    analysis module, but may be in the future.

    Parameters
    ----------
    compound_id: str
        The unique ID of the compound.
    smiles: str | None
        The SMILES string of the compound, if known.
    """

    def __init__(self, compound_id: str, smiles: str | None = None):
        self._smiles = smiles
        super().__init__(compound_id=compound_id)

    def to_dict(self) -> dict[str, str]:
        """
        Convert this ToolCompound to a dict of info required to recreate it.

        Returns
        -------
        dict
            A dictionary representation of the compound.
        """
        if self._smiles is None:
            return {
                "compound_id": self.compound_id,
            }
        else:
            return {
                "compound_id": self.compound_id,
                "smiles": self._smiles,
            }

    @classmethod
    def from_dict(cls, data: dict) -> "ToolCompound":
        """
        Create a ToolCompound from a dictionary representation.

        Parameters
        ----------
        data: dict
            A dictionary representation of the compound.

        Returns
        -------
        ToolCompound
            The created ToolCompound object.
        """
        return cls(compound_id=data["compound_id"], smiles=data.get("smiles", None))


class DopedToolCompound(ToolCompound):
    """
    A ToolCompound that is doped into a DEL.

    Doped compounds are tool compounds that are intentionally added
    to libraries. They require additional information about building
    block tags inorder to call it from the library.

    See the tool compound docs for more info

    Parameters
    ----------
    compound_id: str
        The unique ID of the compound.
    bb_tags: tuple[str, ...]
        the building block tags associated with the doped compound
    smiles: str | None
        The SMILES string of the compound, if known.
    """

    def __init__(self, compound_id: str, bb_tags: tuple[str, ...], smiles: str | None = None):
        self.bb_tags = bb_tags
        super().__init__(compound_id=compound_id, smiles=smiles)

    def to_dict(self) -> dict[str, str]:
        """
        Convert this DopedToolCompound to a dict of info required to recreate it.

        Returns
        -------
        dict
            A dictionary representation of the compound.
        """
        data = super().to_dict()
        data["bb_tags"] = ",".join(self.bb_tags)
        return data

    @classmethod
    def from_dict(cls, data: dict) -> "DopedToolCompound":
        """
        Create a DopedToolCompound from a dictionary representation.

        Parameters
        ----------
        data: dict
            A dictionary representation of the compound.

        Returns
        -------
        DopedToolCompound
            The created DopedToolCompound object.
        """
        return cls(
            compound_id=data["compound_id"],
            bb_tags=tuple(data["bb_tags"].split(",")),
            smiles=data.get("smiles", None),
        )


class TaggedToolCompound(ToolCompound):
    """
    A ToolCompound that has associated DNA tag.

    Tagged compounds are tool compounds that are not doped
    into a library, but are separate and have a DNA tag associated
    to enable their identification during decoding.

    See the tool compound docs for more info

    Parameters
    ----------
    compound_id: str
        The unique ID of the compound.
    tag: str
        the building block tags associated with the tagged compound
    smiles: str | None
        The SMILES string of the compound, if known.
    """

    def __init__(self, compound_id: str, tag: str, smiles: str | None = None):
        super().__init__(compound_id=compound_id, smiles=smiles)
        self.tag = tag

    def to_dict(self) -> dict[str, str]:
        """
        Convert this ToolCompound to a dict of info required to recreate it.

        Returns
        -------
        dict
            A dictionary representation of the compound.
        """
        data = super().to_dict()
        data["tag"] = self.tag
        return data

    @classmethod
    def from_dict(cls, data: dict) -> "TaggedToolCompound":
        """
        Create a ToolCompound from a dictionary representation.

        Parameters
        ----------
        data: dict
            A dictionary representation of the compound.

        Returns
        -------
        TaggedToolCompound
            The created ToolCompound object.
        """
        return cls(compound_id=data["compound_id"], tag=data["tag"], smiles=data.get("smiles", None))


class ToolCompoundParsingError(Exception):
    """Raised when a tool compound cannot be parsed from input data."""

    pass


class ToolCompoundLibrary(Library[ToolCompound], DeliDataLoadable):
    """
    A collection of ToolCompounds for use in decoding.

    This class manages a set of ToolCompounds, allowing for easy
    loading, saving, and retrieval of tool compounds by their IDs.

    Parameters
    ----------
    library_id: str
        The unique ID of the tool compound library.
    tool_compounds: Sequence[ToolCompound]
        The tool compounds in the tool library.
    """

    def __init__(self, library_id: str, tool_compounds: Sequence[ToolCompound]):
        super().__init__(library_id=library_id)
        self.compounds = tool_compounds

        self._id_to_compound: dict[str, ToolCompound] = {tc.compound_id: tc for tc in tool_compounds}

    @classmethod
    @accept_deli_data_name("tool_compounds", "json", target_param="name_or_path")
    def load(cls, name_or_path: str):
        """
        Load a ToolCompoundLibrary from a JSON file.

        Notes
        -----
        Can load names of tool compounds from the deli data directory subdir "tool_compounds".

        Parameters
        ----------
        name_or_path: str
            The name or path of the JSON file to load.

        Returns
        -------
        ToolCompoundLibrary
            The loaded ToolCompoundLibrary.
        """
        import json
        import os

        with open(name_or_path, "r") as f:
            data = json.load(f)

        tool_compounds: list[ToolCompound] = list()
        for tc_data in data["tool_compounds"]:
            if "compound_id" not in tc_data:
                raise ToolCompoundParsingError(f"tool compound missing 'compound_id' field from file '{name_or_path}'")
            else:
                tool_compounds.append(ToolCompound.from_dict(tc_data))

        library_id = ".".join(os.path.basename(name_or_path).split(".")[:-1])  # filename without extension

        return cls(library_id=library_id, tool_compounds=tool_compounds)

    def get_compound(self, compound_id: str) -> ToolCompound:
        """
        Retrieve a ToolCompound by its ID.

        Parameters
        ----------
        compound_id: str
            The ID of the compound to retrieve.

        Returns
        -------
        ToolCompound
            The ToolCompound.

        Raises
        ------
        KeyError
            If the compound ID is not found in the library.
        """
        try:
            return self._id_to_compound[compound_id]
        except KeyError as e:
            raise KeyError(f"ToolCompound with ID {compound_id} not found in tool library") from e


class TaggedToolCompoundLibrary(ToolCompoundLibrary, BarcodedMixin[ToolCompoundBarcodeSchema]):
    """
    A collection of ToolCompounds for use in decoding.

    This class manages a set of ToolCompounds, allowing for easy
    loading, saving, and retrieval of tool compounds by their IDs.

    Parameters
    ----------
    tool_compounds: Sequence[TaggedToolCompound]
        The tool compounds in the tool library.

    Attributes
    ----------
    library_tag: str
        The library DNA tag for the tool compound library.
    """

    def __init__(
        self, library_id: str, tool_compounds: Sequence[TaggedToolCompound], barcode_schema: ToolCompoundBarcodeSchema
    ):
        super().__init__(library_id, tool_compounds)
        self.barcode_schema = barcode_schema
        self.compounds: Sequence[TaggedToolCompound] = tool_compounds  # type hinting override

        self.library_tag = self.barcode_schema.library_section.section_tag
        self._tag_to_compound: dict[str, TaggedToolCompound] = {tc.tag: tc for tc in tool_compounds}

    @classmethod
    @accept_deli_data_name("tool_compounds", "json", target_param="name_or_path")
    def load(cls, name_or_path: str):
        """
        Load a ToolCompoundLibrary from a JSON file.

        Notes
        -----
        Can load names of tool compounds from the deli data directory subdir "tool_compounds".

        Parameters
        ----------
        name_or_path: str
            The name or path of the JSON file to load.

        Returns
        -------
        ToolCompoundLibrary
            The loaded ToolCompoundLibrary.
        """
        import json
        import os

        with open(name_or_path, "r") as f:
            data = json.load(f)

        barcode_schema = ToolCompoundBarcodeSchema.from_dict(data["barcode_schema"])
        tool_compounds: list[TaggedToolCompound] = list()
        for tc_data in data["tool_compounds"]:
            if "compound_id" not in tc_data:
                raise ToolCompoundParsingError(f"tool compound missing 'compound_id' field from file '{name_or_path}'")
            elif "tag" not in tc_data:
                raise ToolCompoundParsingError(
                    f"tagged tool compound {tc_data['compound_id']} missing 'tag' field from file '{name_or_path}'"
                )
            else:
                tool_compounds.append(TaggedToolCompound.from_dict(tc_data))
        library_id = ".".join(os.path.basename(name_or_path).split(".")[:-1])  # filename without extension

        return cls(library_id=library_id, tool_compounds=tool_compounds, barcode_schema=barcode_schema)

    @overload
    def search_tags(self, query: str, fail_on_missing: Literal[False]) -> TaggedToolCompound | None: ...

    @overload
    def search_tags(self, query: str, fail_on_missing: Literal[True]) -> TaggedToolCompound: ...

    def search_tags(self, query: str, fail_on_missing: bool = False) -> TaggedToolCompound | None:
        """
        Search for a TaggedToolCompound by its tag.

        Parameters
        ----------
        query: str
            The tag to search for.
        fail_on_missing: bool
            Whether to raise an error if the tag is not found.

        Returns
        -------
        TaggedToolCompound | None
            The TaggedToolCompound if found, else None.
            will only return TaggedToolCompound if fail_on_missing is True.

        Raises
        ------
        KeyError
            If the tag is not found and error_on_missing is True.
        """
        try:
            return self._tag_to_compound[query]
        except KeyError as e:
            if fail_on_missing:
                raise KeyError(
                    f"TaggedToolCompound with tag {query} not found in tool library with library tag "
                    f"'{self.barcode_schema.library_section.section_tag}'"
                ) from e
            else:
                return None
