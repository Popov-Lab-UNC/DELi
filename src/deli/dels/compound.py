"""for classes relating to DEL compounds"""

from typing import Optional

from deli.dels import BuildingBlock, DELCollection, DELibrary


class DELCompoundException(Exception):
    """for errors related to DEL compounds"""

    pass


class Compound:
    """
    Base compound class

    This class is used to represent a compound in a DEL.
    It contains the compound ID, library ID, and building blocks.
    """

    def __init__(self, compound_id: str):
        """
        Initialize the BaseCompound object.

        Parameters
        ----------
        compound_id: str
            The ID of the compound.
        """
        self.compound_id = compound_id

    def __eq__(self, other):
        """Check if two LowMemDELCompound objects are equal."""
        if not isinstance(other, LowMemDELCompound):
            return False
        return self.compound_id == other.compound_id

    def __hash__(self):
        """Return the hash of the object."""
        return hash(self.compound_id)

    def __repr__(self):
        """Return a string representation of the object."""
        return f"{self.__class__.__name__}: {self.compound_id}"


class DELCompound(Compound):
    """
    DEL compound class

    This class is used to represent a compound in a DEL.
    It contains info about the compounds origin in DEL,
    like which DELibrary it belongs to, and the building blocks used to create it.

    Additional info can be generated (or passed directly during initialization),
    like compound_id or enumerated SMILES.
    """

    def __init__(
        self,
        library: DELibrary,
        building_blocks: list[BuildingBlock],
        compound_id: Optional[str] = None,
    ):
        """
        Initialize the DELCompound object.

        Parameters
        ----------
        compound_id: str
            The ID of the compound.
        library: DELibrary
            The DELibrary object that this compound belongs to.
        building_blocks: dict[str, BuildingBlock]
            The building blocks of the compound
            should be in cycle order (cycle 1 first, then cycle 2, etc.)
        """
        self.library = library
        self.building_blocks = building_blocks

        if compound_id is None:
            _compound_id = f"{self.library.library_id}-" + "-".join(
                [bb.bb_id for bb in self.building_blocks]
            )
        else:
            _compound_id = compound_id

        super().__init__(_compound_id)

    def to_low_mem(self) -> "LowMemDELCompound":
        """
        Convert this DELCompound to a LowMemDELCompound.

        Returns
        -------
        LowMemDELCompound
            A low-memory representation of this compound.
        """
        return LowMemDELCompound(
            compound_id=self.compound_id,
            library_id=self.library.library_id,
            building_blocks_ids=[bb.bb_id for bb in self.building_blocks],
        )


class LowMemDELCompound(Compound):
    """
    Low memory DEL compound class

    Rather than pointing to objects (like DELibrary or BuildingBlocks),
    this class stores info about a compound as the object string ids

    This is useful for large DELs where:
    - memory usage is a concern
    - compounds need to be used as keys in a dictionary or set
    - you need to save info about the compound in a non-binary file
    """

    def __init__(
        self,
        library_id: str,
        building_blocks_ids: list[str],
        compound_id: str,
    ):
        """
        Initialize the LowMemDELCompound object.

        Parameters
        ----------
        library_id: str
            The ID of the library.
        building_blocks_ids: list[str]
            The IDs of the building blocks.
        """
        self.library_id = library_id
        self.building_blocks_ids = building_blocks_ids
        if compound_id is None:
            _compound_id = f"{self.library_id}-" + "-".join(
                [bb_id for bb_id in self.building_blocks_ids]
            )
        else:
            _compound_id = compound_id

        super().__init__(_compound_id)

    def load_compound(self, collection: DELCollection) -> DELCompound:
        """
        Load the full DELCompound object from this low memory representation.

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
            library = collection.get_library(self.library_id)
        except KeyError as e:
            raise DELCompoundException(f"Library {self.library_id} not found in collection") from e

        _bbs: list[BuildingBlock] = list()
        for i, (bb_set, bb_id) in enumerate(zip(library.bb_sets, self.building_blocks_ids)):
            try:
                _bbs.append(bb_set.get_bb_by_id(bb_id, fail_on_missing=True))
            except KeyError as e:
                raise DELCompoundException(
                    f"Building block {bb_id} for cycle {i+1} "
                    f"not found in library {self.library_id}"
                ) from e
        return DELCompound(library=library, building_blocks=_bbs, compound_id=self.compound_id)

    def to_dict(self) -> dict:
        """
        Convert this LowMemDELCompound to a dict of info required to recreate it.

        Returns
        -------
        dict
            A dictionary representation of the compound.
        """
        return {
            "library_id": self.library_id,
            "building_blocks_ids": self.building_blocks_ids,
        }
