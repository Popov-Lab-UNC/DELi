"""code for calling the building blocks of DEL reads"""

from typing import TypeAlias, Union

from deli.dels import BuildingBlock, BuildingBlockBarcodeSection, BuildingBlockSet

from .calls import FailedCall, ValidCall


class ValidBuildingBlockCall(ValidCall):
    """base class for all calls"""

    def __init__(self, building_block: BuildingBlock, score: Union[float, int]):
        self.building_block = building_block
        self.score = score


class FailedBuildingBlockCall(FailedCall):
    """call  was found"""

    pass


# type alias for building block calling
BuildingBlockCall: TypeAlias = Union[ValidBuildingBlockCall, FailedBuildingBlockCall]


class BuildingBlockSetTagCaller:
    """Calls building blocks given a tag query and building block set"""

    def __init__(
        self,
        building_block_tag_section: BuildingBlockBarcodeSection,
        building_block_set: BuildingBlockSet,
        use_hamming: bool = True,
    ):
        """
        Initialize a BuildingBlockCaller

        Parameters
        ----------
        building_block_tag_section: BuildingBlockBarcodeSection
            the building block tag section from barcode to decode from
        building_block_set: BuildingBlockSet
            the building block set to look for matches in
        use_hamming: bool
            use hamming decoding if tag is hamming encoded
            if tag is not hamming encoded, will skip hamming decoding
        """
        self.building_block_tag_section = building_block_tag_section
        self.building_block_set = building_block_set
        self.use_hamming = use_hamming and self.building_block_tag_section.is_hamming_encoded()

    def call_building_block(self, tag: str) -> BuildingBlockCall:
        """
        Call a building block given a nucleotide tag query

        Parameters
        ----------
        tag: str
            the nucleotide tag query to search

        Returns
        -------
        BuildingBlockCall
            the building block call
            will be a ValidBuildingBlockCall if match is found
            will be a FailedBuildingBlockCall if match is not found
        """
        if self.use_hamming and self.building_block_tag_section.hamming_decoder is not None:
            _tag = self.building_block_tag_section.hamming_decoder.decode_sequence(tag)
            if _tag is None:
                return FailedBuildingBlockCall()
            else:
                tag = _tag
        _call = self.building_block_set.search_tags(tag, fail_on_missing=False)

        if _call is None:
            return FailedBuildingBlockCall()
        else:
            # score is always 0 for building block matches
            return ValidBuildingBlockCall(_call, 0)
