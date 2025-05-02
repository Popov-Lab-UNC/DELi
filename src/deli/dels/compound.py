"""for classes relating to DEL compounds"""


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
            must be same length as building_blocks_ids if not None
        enumerated_smi: str | None, default None
            The enumerated SMILES representation of the compound.
            If unknown, set to None.
        """
        self.compound_id = compound_id
        self.library_id = library_id
        self.building_blocks_ids = building_blocks_ids
        self.building_block_smiles = building_block_smis
        self.enumerated_smi = enumerated_smi

    def __repr__(self):
        """Return a string representation of the object."""
        return self.compound_id
