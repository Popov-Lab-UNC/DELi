"""for classes relating to DEL compounds"""

import abc


class ChemicalCompoundException(Exception):
    """for errors related to chemical compounds"""

    pass


class BaseCompound(abc.ABC):
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
        return self.compound_id

    @abc.abstractmethod
    def get_smiles(self) -> str:
        """
        Get the SMILES representation of the full compound

        Returns
        -------
        str

        Raises
        ------
        ChemicalCompoundException
            if the SMILES cannot be generated with provided data
            and is missing
        """
        raise NotImplementedError


class LowMemDELCompound:
    """
    Low memory DEL compound class

    Rather than pointing to objects, this class stores
    info about a compound as strings, based on the ids
    """

    def __init__(
        self,
        compound_id: str,
        library_id: str,
        building_blocks_ids: list[str],
        building_block_smis: list[str | None] | None = None,
        enumerated_smi: str | None = None,
    ):
        """
        Initialize the LowMemDELCompound object.

        Parameters
        ----------
        compound_id: str
            The ID of the compound.
        library_id: str
            The ID of the library.
        building_blocks_ids: list[str]
            The IDs of the building blocks.
        building_block_smis: list[str | None] | None, default None
            The SMILES representations of the building blocks.
            Must be the same length as building_blocks_ids if not None
        enumerated_smi: str | None, default None
            The enumerated SMILES representation of the compound.
            If unknown, set to None.
        """
        self.compound_id = compound_id
        self.library_id = library_id
        self.building_blocks_ids = building_blocks_ids
        self.building_block_smiles = building_block_smis
        self.enumerated_smi = enumerated_smi
