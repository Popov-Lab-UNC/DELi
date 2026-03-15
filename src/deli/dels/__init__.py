"""Public exports for DEL classes and functions."""

from .barcode import (
    BarcodedMixin,
    BarcodeSchema,
    BarcodeSchemaError,
    BarcodeSection,
    BuildingBlockBarcodeSection,
    DecodeableBarcodeSection,
    DELBarcodeSchema,
    LibraryBarcodeSection,
    MixedBarcodeSection,
    PrimerBarcodeSection,
    StaticBarcodeSection,
    ToolCompoundBarcodeSchema,
    ToolCompoundRefBarcodeSection,
    UMIBarcodeSection,
    VariableBarcodeSection,
)
from .building_block import (
    BuildingBlock,
    BuildingBlockSet,
    BuildingBlockSetError,
    NullBuildingBlock,
    TaggedBuildingBlock,
    TaggedBuildingBlockSet,
    TaggedFakeBuildingBlock,
    TaggedNullBuildingBlock,
)
from .combinatorial import (
    CombinatorialCollection,
    CombinatorialLibrary,
    CombinatorialLibraryCollection,
    DELibrary,
    DELibraryCollection,
    LibraryBuildError,
)
from .compound import (
    Compound,
    DELCompound,
    DELCompoundException,
    DELCompoundRaw,
    generate_del_compound_id,
)
from .library import Library, LibraryCollection
from .tool_compound import (
    DopedToolCompound,
    TaggedToolCompound,
    ToolCompound,
    ToolCompoundParsingError,
)


__all__ = [
    "BarcodeSchemaError",
    "BarcodeSection",
    "VariableBarcodeSection",
    "DecodeableBarcodeSection",
    "BuildingBlockBarcodeSection",
    "UMIBarcodeSection",
    "StaticBarcodeSection",
    "PrimerBarcodeSection",
    "LibraryBarcodeSection",
    "ToolCompoundRefBarcodeSection",
    "MixedBarcodeSection",
    "BarcodeSchema",
    "DELBarcodeSchema",
    "ToolCompoundBarcodeSchema",
    "BarcodedMixin",
    "BuildingBlockSetError",
    "BaseBuildingBlock",
    "BuildingBlock",
    "NullBuildingBlock",
    "TaggedBuildingBlock",
    "TaggedNullBuildingBlock",
    "TaggedFakeBuildingBlock",
    "BuildingBlockSet",
    "TaggedBuildingBlockSet",
    "LibraryBuildError",
    "CombinatorialLibrary",
    "DELibrary",
    "CombinatorialCollection",
    "CombinatorialLibraryCollection",
    "DELibraryCollection",
    "generate_del_compound_id",
    "DELCompoundException",
    "Compound",
    "DELCompoundRaw",
    "DELCompound",
    "Library",
    "LibraryCollection",
    "ToolCompound",
    "DopedToolCompound",
    "TaggedToolCompound",
    "ToolCompoundParsingError",
]


def load_from_component_id(component_id: str, id_separator: str | None = None) -> DELCompound:
    """
    Given a DEL compound ID in component ID format, load the corresponding compound object.

    Will attempt to load from the configured DELi data directory and fail if the compound cannot be
    found within the directory.

    Parameters
    ----------
    component_id: str
        the DEL compound ID in component ID format, e.g. "DEL004-1-23-45"
    id_separator: str, optional
        the separator used in the component ID format, if not provided will use the value from the DELi config

    Returns
    -------
    DELCompound
        the corresponding DELCompound object
    """
    if id_separator is None:
        from deli.configure import get_deli_config

        id_separator = get_deli_config().comp_id_sep
    comp_ids = component_id.split(id_separator)
    if len(comp_ids) < 2:
        raise ValueError(
            f"Invalid component ID format: {component_id}. "
            f"Expected at least library ID and one building block ID separated by '{id_separator}'."
        )
    library = DELibrary.load(comp_ids[0])

    return library.get_compound(tuple(comp_ids[1:]))
