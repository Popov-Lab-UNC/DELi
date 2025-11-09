"""for classes relating to DEL compounds"""

from typing import TYPE_CHECKING

from .building_block import BuildingBlock


# for mypy to recognize the library type hints
if TYPE_CHECKING:
    from deli.enumeration.enumerator import EnumeratedDELCompound

    from .library import Library, LibraryCollection


def generate_compound_id(library_id: str, building_blocks: list[str]) -> str:
    """
    Generate a compound ID from the library ID and building block IDs.

    Parameters
    ----------
    library_id: str
        The ID of the library.
    building_blocks: list[str]
        The IDs of the building blocks.

    Returns
    -------
    str
        The generated compound ID.
    """
    return f"{library_id}-" + "-".join(building_blocks)


class DELCompoundException(Exception):
    """for errors related to DEL compounds"""

    pass


class DELCompoundRaw:
    """
    A DEL compound that only stores minimal information to identify it.

    In practice, this means it only stores the library ID and the building block IDs.
    This is useful in cases where you lack access to the library files ubt want to load
    in DEL compounds.

    Attributes
    ----------
    self.compound_id: str
        The unique ID of the compound, generated from the library ID and building block IDs.
        Generated using the `generate_compound_id` function.

    Notes
    -----
    Because of the minimal information stored, this class does not
    check that the building blocks are valid or that the library exists.
    """

    def __init__(
        self,
        library_id: str,
        building_blocks_ids: list[str],
    ):
        """
        Initialize the DELCompoundRaw object.

        Parameters
        ----------
        library_id: str
            The ID of the library.
        building_blocks_ids: list[str]
            The IDs of the building blocks.
        """
        self.library_id = library_id
        self.building_block_ids = building_blocks_ids
        self.compound_id = generate_compound_id(self.library_id, self.building_block_ids)

    def __eq__(self, other):
        """Check if two Compound objects are equal (share the same compound_id)"""
        if not isinstance(other, DELCompound):  # any child of this class can be equal
            return False
        return self.compound_id == other.compound_id

    def __hash__(self):
        """Return the hash of the object."""
        return hash(self.compound_id)

    def __repr__(self):
        """Return a string representation of the object."""
        return f"{self.__class__.__name__}({self.compound_id})"

    def load_compound(self, collection: "LibraryCollection") -> "DELCompound":
        """
        Load the full DELCompound information

        Will attempt to load the information about the library and building blocks
        from the provided collection.

        Parameters
        ----------
        collection: DELibraryCollection
            The collection of libraries to load the compound from.

        Returns
        -------
        DELCompound
            The full DELCompound object.

        Raises
        ------
        DELCompoundException
            when the library or building blocks cannot be found to create the compound.
        """
        try:
            library: "Library" = collection.get_library(self.library_id)
        except KeyError as e:
            raise DELCompoundException(f"Library {self.library_id} not found in collection") from e

        _bbs: list[BuildingBlock] = list()
        for i, (bb_set, bb_id) in enumerate(
            zip(library.bb_sets, self.building_block_ids, strict=False)
        ):
            try:
                _bbs.append(bb_set.get_bb_by_id(bb_id, fail_on_missing=True))
            except KeyError as e:
                raise DELCompoundException(
                    f"Building block {bb_id} for cycle {i + 1} "
                    f"not found in library {self.library_id}"
                ) from e
        return DELCompound(library=library, building_blocks=_bbs)

    def to_dict(self) -> dict:
        """
        Convert this DELCompoundRaw to a dict of info required to recreate it.

        Returns
        -------
        dict
            A dictionary representation of the compound.
        """
        return {
            "library_id": self.library_id,
            "building_blocks_ids": self.building_block_ids,
        }


class DELCompound(DELCompoundRaw):
    """
    DEL compound class

    This class is used to represent a compound in a DEL.
    It contains info about the compounds origin in DEL,
    like which DELibrary it belongs to, and the building blocks used to create it.

    Additional info can be generated (or passed directly during initialization),
    like compound_id or enumerated SMILES.

    Attributes
    ----------
    compound_id: str
        The unique ID of the compound, generated from the library ID and building block IDs.
        Generated using the `generate_compound_id` function.
    library_id: str
        The ID of the library this compound belongs to.
    building_block_ids: list[str]
        The IDs of the building blocks used to create this compound.

    See Also
    --------
    generate_compound_id : function to generate compound IDs
    """

    def __init__(
        self,
        library: "Library",
        building_blocks: list[BuildingBlock],
    ):
        """
        Initialize the DELCompound object.

        Parameters
        ----------
        library: DELibrary
            The DELibrary object that this compound belongs to.
        building_blocks: list[BuildingBlock]
            The building blocks of the compound
            should be in cycle order (cycle 1 first, then cycle 2, etc.)
        """
        self.library = library
        self.building_blocks = building_blocks

        super().__init__(
            library_id=library.library_id, building_blocks_ids=[bb.bb_id for bb in building_blocks]
        )

    def to_raw(self) -> "DELCompoundRaw":
        """
        Convert this DELCompound to a LowMemDELCompound.

        Returns
        -------
        DELCompoundRaw
            A low-memory representation of this compound.
        """
        return DELCompoundRaw(
            library_id=self.library.library_id,
            building_blocks_ids=[bb.bb_id for bb in self.building_blocks],
        )

    def enumerate(self) -> "EnumeratedDELCompound":
        """
        Enumerate this DELCompound to get its SMILES representation.

        Will use the enumerator attached to the library to perform the enumeration.

        Returns
        -------
        EnumeratedDELCompound
            The enumerated DEL compound with SMILES.

        Raises
        ------
        DELCompoundException
            when enumeration fails for any reason.
        """
        return self.library.enumerator.enumerate_by_bbs(self.building_blocks)
