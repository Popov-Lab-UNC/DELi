"""for classes relating to DEL compounds"""

import abc
from typing import TYPE_CHECKING, Optional

from rdkit.Chem import Mol

from deli.utils import SmilesMixin

from .building_block import BuildingBlock


# for mypy to recognize the library type hints
if TYPE_CHECKING:
    from .library import Library, LibraryCollection


def _generate_compound_id(library_id: str, building_blocks: list[str]) -> str:
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
        """Check if two Compound objects are equal (share the same compound_id)"""
        if not isinstance(other, Compound):
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
        library: "Library",
        building_blocks: list[BuildingBlock],
    ):
        """
        Initialize the DELCompound object.

        Parameters
        ----------
        library: DELibrary
            The DELibrary object that this compound belongs to.
        building_blocks: dict[str, BuildingBlock]
            The building blocks of the compound
            should be in cycle order (cycle 1 first, then cycle 2, etc.)
        """
        self.library = library
        self.building_blocks = building_blocks

        _compound_id = _generate_compound_id(
            library.library_id, [bb.bb_id for bb in self.building_blocks]
        )

        super().__init__(_compound_id)

    def enumerate(self) -> "EnumeratedDELCompound":
        """
        Attempt to enumerate to get the SMILES for this compound.

        Returns
        -------
        EnumeratedDELCompound
            An EnumeratedDELCompound object with the SMILES string set.

        Raises
        ------
        DELCompoundException
            when the DEL compounds library lacks an enumerator
        EnumerationRunError
            when the enumeration fails for any reason
        """
        return self.library.enumerate_by_bbs(self.building_blocks)

    def to_low_mem(self) -> "LowMemDELCompound":
        """
        Convert this DELCompound to a LowMemDELCompound.

        Returns
        -------
        LowMemDELCompound
            A low-memory representation of this compound.
        """
        return LowMemDELCompound(
            library_id=self.library.library_id,
            building_blocks_ids=[bb.bb_id for bb in self.building_blocks],
        )


class LowMemCompound(Compound, abc.ABC):
    """
    Abstract interface for low memory compound classes

    Low memory compounds are used to avoid excessive pointers/object usage when it is unnecessary.

    Rather than pointing to objects (like DELibrary or BuildingBlocks),
    this class stores info about a compound as the object string ids.
    As long as you have the DELibraryCollection object that originally contained these IDs,
    you can load the full DELCompound object again later on.

    This is useful when you need to save the compounds to a file,
    as they can be saved in an efficient text format.
    They can be converted to the full DELCompound object after loading,
    only requiring the DELibraryCollection object that they originated from.

    An example of where this is useful is when trying to merge several DEL counters
    after parallel decoding runs.
    If the counters are saving the compounds as serialized objects/pickles,
    loading and merging will result in the duplication of all the DELibrary and
    BuildingBlock objects.
    Instead, you can save the compounds as low memory representations,
    merge on the string IDs,
    and reassign pointers to a single object in memory.
    This is after faster and memory efficient, and avoid using pickles,
    which can be a security risk
    """

    @abc.abstractmethod
    def load_compound(self, collection: "LibraryCollection") -> DELCompound:
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
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def to_dict(self) -> dict:
        """
        Convert this LowMemCompound to a dict of info required to recreate it.

        Returns
        -------
        dict
            A dictionary representation of the compound.
        """
        raise NotImplementedError()


class LowMemDELCompound(LowMemCompound):
    """
    Low memory DEL compound class

    Notes
    -----
    Because of the low memory representation, this class does not contain
    checks that the building blocks are valid or that the library exists.
    It is best to create objects of this class by calling "to_low_mem" on the parent
    a DELCompound object.
    """

    def __init__(
        self,
        library_id: str,
        building_blocks_ids: list[str],
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
        _compound_id = _generate_compound_id(self.library_id, self.building_blocks_ids)

        super().__init__(_compound_id)

    def load_compound(self, collection: "LibraryCollection") -> DELCompound:
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
            library: "Library" = collection.get_library(self.library_id)
        except KeyError as e:
            raise DELCompoundException(f"Library {self.library_id} not found in collection") from e

        _bbs: list[BuildingBlock] = list()
        for i, (bb_set, bb_id) in enumerate(zip(library.bb_sets, self.building_blocks_ids)):
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


class LowMemEnumeratedDELCompound(LowMemDELCompound, SmilesMixin):
    """
    Low memory DEL enumerated compound class

    This class is used to represent a DEL compound that has been fully enumerated
    and has a SMILES string mapped to the full compound, represented as string IDs.
    It inherits from LowMemDELCompound and implements the SmilesMixin.

    If Mol is not provided, it will be generated from the SMILES string
    and cached when first accessed.
    """

    def __init__(
        self,
        library_id: str,
        building_blocks_ids: list[str],
        smiles: str,
        mol: Optional[Mol] = None,
    ):
        """
        Initialize the LowMemEnumeratedDELCompound object.

        Parameters
        ----------
        library_id: str
            The ID of the library.
        building_blocks_ids: list[str]
            The IDs of the building blocks.
        smiles: str
            The SMILES string of the compound.
        mol: Optional[Mol], default=None
            The RDKit Mol object for the compound.
            If None, it will be generated from the SMILES string and
            cached when first accessed.
            NOTE: it is not recommended to set this directly unless
                  you are confident the mol object
                  is the results of the provided SMILES
        """
        super().__init__(library_id=library_id, building_blocks_ids=building_blocks_ids)
        self._smiles = smiles
        self._mol = mol


class EnumeratedDELCompound(DELCompound, SmilesMixin):
    """
    DEL enumerated compound class

    This class is used to represent a DEL compound that has been fully enumerated
    and has a SMILES string mapped to the full compound.
    It inherits from DELCompound and implements the SmilesMixin.

    If Mol is not provided, it will be generated from the SMILES string
    and cached when first accessed.
    """

    def __init__(
        self,
        library: "Library",
        building_blocks: list[BuildingBlock],
        smiles: str,
        mol: Optional[Mol] = None,
    ):
        """
        Initialize the EnumeratedDELCompound object.

        Parameters
        ----------
        library: DELibrary
            The DELibrary object that this compound belongs to.
        building_blocks: dict[str, BuildingBlock]
            The building blocks of the compound
            should be in cycle order (cycle 1 first, then cycle 2, etc.)
        smiles: str
            The SMILES string of the compound.
        mol: Optional[Mol], default=None
            The RDKit Mol object for the compound.
            If None, it will be generated from the SMILES string
            and cached when first accessed.
            NOTE: it is not recommended to set this directly
                  unless you are confident the mol object
                  is the results of the provided SMILES
        """
        super().__init__(library=library, building_blocks=building_blocks)
        self._smiles = smiles
        self._mol = mol
