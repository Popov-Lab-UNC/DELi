"""Handle tool compounds for DELi"""

from deli.configure import DeliDataLoadable, accept_deli_data_name

from .barcode import ToolCompoundBarcodeSchema
from .compound import Compound


class ToolCompound(Compound, DeliDataLoadable):
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
    @accept_deli_data_name("tool_compounds", "json", target_param="name_or_path")
    def load(cls, name_or_path: str, load_smiles: bool = True) -> "ToolCompound":
        """
        Load a ToolCompound from a Tool Compound JSON file.

        See the docs for more info on the expected format.

        Will only load information related to the Tool compound from the file.
        Information related to any DNA tag will be ignored even if present in the file.
        If this info is needed, use the `TaggedToolCompound.load()` method instead.

        Notes
        -----
        Can load names of tool compounds from the deli data directory subdir "tool_compounds".

        Parameters
        ----------
        name_or_path: str
            The name or path of the JSON file to load.
        load_smiles: bool, default=True
            Whether to load the SMILES string from the file, if present.

        Returns
        -------
        ToolCompound
            The loaded ToolCompound.
        """
        import json

        with open(name_or_path, "r") as f:
            data = json.load(f)

        # quick validation
        if "compound_id" not in data:
            raise ToolCompoundParsingError(f"tool compound missing 'compound_id' field from file '{name_or_path}'")

        return cls(compound_id=data["compound_id"], smiles=data.get("smiles", None) if load_smiles else None)


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

    @classmethod
    def load(cls, *args, **kwargs):
        """
        Diable loading for DopedToolCompounds.

        Unlike tool compounds, DopedToolCompounds cannot be loaded directly.
        Loading must always be handled by the CombinatorialLibrary they are doped into.
        """
        raise NotImplementedError(
            "DopedToolCompound are loaded via CombinatorialLibrary objects; they cannot be loaded directly"
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
    barcode_schema: ToolCompoundBarcodeSchema
        The barcode schema associated with the tagged tool compound
    smiles: str | None
        The SMILES string of the compound, if known.
    """

    def __init__(self, compound_id: str, barcode_schema: ToolCompoundBarcodeSchema, smiles: str | None = None):
        super().__init__(compound_id=compound_id, smiles=smiles)
        self.barcode_schema = barcode_schema

    @classmethod
    @accept_deli_data_name("tool_compounds", "json", target_param="name_or_path")
    def load(cls, name_or_path: str, load_smiles: bool = True) -> "TaggedToolCompound":
        """
        Load a TaggedToolCompound from a Tool Compound JSON file.

        See the docs for more info on the expected format.

        Will only load information related to the Tagged Tool compound from the file.
        Information related to any doping/building block tags will be ignored even if present in the file.
        If this info is needed, use the `DopedToolCompound.load()` method instead.

        Notes
        -----
        Can load names of tool compounds from the deli data directory subdir "tool_compounds".

        Parameters
        ----------
        name_or_path: str
            The name or path of the JSON file to load.
        load_smiles: bool, default=True
            Whether to load the SMILES string from the file, if present.

        Returns
        -------
        TaggedToolCompound
            The loaded TaggedToolCompound.
        """
        import json

        with open(name_or_path, "r") as f:
            data = json.load(f)

        # quick validation
        if "compound_id" not in data:
            raise ToolCompoundParsingError(
                f"tagged tool compound missing 'compound_id' field from file '{name_or_path}'"
            )
        if "barcode_schema" not in data:
            raise ToolCompoundParsingError(
                f"tagged tool compound missing 'barcode_schema' field from file '{name_or_path}'"
            )

        barcode_schema = ToolCompoundBarcodeSchema.from_dict(data["barcode_schema"])

        return cls(
            compound_id=data["compound_id"],
            barcode_schema=barcode_schema,
            smiles=data.get("smiles", None) if load_smiles else None,
        )


class ToolCompoundParsingError(Exception):
    """Raised when a tool compound cannot be parsed from input data."""

    pass
