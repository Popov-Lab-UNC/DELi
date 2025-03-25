"""code for calling DNA barcodes"""

import abc
import warnings
from typing import Literal

from dnaio import SequenceRecord

from deli.dels import DELibrary, DELibraryPool
from deli.dna import Aligner, HybridSemiGlobalAligner, SemiGlobalAligner

from .bb_calling import BuildingBlockCall, BuildingBlockSetTagCaller
from .lib_calling import LibraryCall, LibraryCaller, SingleReadLibraryCaller, ValidLibraryCall
from .settings import DecodingSettings
from .umi import UMI


class DecodedBarcode:
    """
    Holds information about a decoded barcode

    Can handle both successful and failed decodes
    """

    def __init__(
        self,
        library: LibraryCall,
        building_blocks: dict[str, BuildingBlockCall] | None,
        umi: UMI | None = None,
    ):
        """
        initialize the DecodedBarcode object

        Parameters
        ----------
        library : LibraryCall
            the called library
        building_blocks : dict[str, BuildingBlockCall] or None
            the building block calls
            if no calls were attempted (because library failed) pass a `None`
        umi: UMI or `None`
            the UMI for the read
            if not using umi or no umi in the barcode, use a `None`
        """
        self.library = library
        self.building_blocks = building_blocks
        self.umi = umi

    def is_decoded_successfully(self) -> bool:
        """
        True if all required parts of the barcode were successfully decoded else False

        Returns
        -------
        bool
        """
        if self.building_blocks is not None and all(
            [bb.is_valid() for bb in self.building_blocks.values()]
        ):
            return True
        return False


# TODO add a single library decoder class


class DELPoolDecoder:
    """
    A decoder used to decode barcodes from selection using a pool of libraries

    For most use cases, this is the main entry point for decoding.
    Use one of the LibraryDecoders if you only have one (1) DEL in your selection
    As this will be far faster
    """

    def __init__(
        self, library_pool: DELibraryPool, decode_settings: DecodingSettings | None = None
    ):
        """
        Initialize the DELPoolDecoder

        Parameters
        ----------
        library_pool: DELibraryPool
            the library pool used for the selection to decode
        decode_settings: DecodingSettings or None
            the settings used to decode the barcodes
            if none are passed, will use default settings
            see Decoding docs for more details
        """
        self.decode_settings: DecodingSettings
        if decode_settings is None:
            self.decode_settings = DecodingSettings()
        else:
            self.decode_settings = decode_settings

        self.library_caller: LibraryCaller
        if decode_settings.__getattribute__("read_type") == "single":
            self.library_caller = SingleReadLibraryCaller(
                library_pool,
                error_rate=decode_settings.__getattribute__("error_rate"),
                min_overlap=decode_settings.__getattribute__("min_overlap"),
                revcomp=decode_settings.__getattribute__("revcomp"),
            )
        elif decode_settings.__getattribute__("read_type") == "paired":
            warnings.warn(
                "DELi only support single reads currently; if a paired mode would be nice"
                "please raise an issue to request it",
                stacklevel=1,
            )
            self.library_caller = SingleReadLibraryCaller(
                library_pool,
                error_rate=decode_settings.__getattribute__("error_rate"),
                min_overlap=decode_settings.__getattribute__("min_overlap"),
                revcomp=decode_settings.__getattribute__("revcomp"),
            )
        else:
            raise ValueError(
                f"Unknown read_type '{decode_settings.__getattribute__('read_type')}'"
            )

        # determine library caller class
        if decode_settings.__getattribute__("bb_calling_approach") == "alignment":
            self.library_decoders: dict[str, LibraryDecoder] = {
                library.library_id: DynamicAlignmentLibraryDecoder(
                    library,
                    use_hamming=decode_settings.__getattribute__("use_hamming"),
                    alignment_algorithm=decode_settings.__getattribute__("alignment_algorithm"),
                )
                for library in library_pool.libraries
            }
        else:
            raise ValueError(
                f"unrecognized bb_calling_approach "
                f"{decode_settings.__getattribute__('bb_calling_approach')}"
            )

    def decode_read(self, sequence: SequenceRecord) -> DecodedBarcode:
        """Given a sequence read, decode it's barcode"""
        library_call = self.library_caller.call_library(sequence)
        if isinstance(library_call, ValidLibraryCall):
            return self.library_decoders[library_call.library.library_id].call_barcode(
                library_call
            )
        else:
            # if library calling fails cannot continue decoding
            return DecodedBarcode(library_call, None, None)


class LibraryDecoder(abc.ABC):
    """
    Abstract class for all library decoders

    A library decoder assumes the library the read is from is already known.
    These decoders are only concerned with calling or locating all other sections.
    Currently that is the building block sections and UMI.
    If you would like support for calling other section too, raise an issue
    to add that feature
    """

    def __init__(self, library: DELibrary, use_hamming: bool = True):
        """
        Initialize the LibraryDecoder

        Parameters
        ----------
        library: DELibrary
            the library to decode from
        use_hamming: bool, default True
            using hamming correction
            only used for barcodes with hamming encoded parts
        """
        self.use_hamming = use_hamming
        self.library = library

    @abc.abstractmethod
    def call_barcode(self, library_call: ValidLibraryCall) -> DecodedBarcode:
        """Given a library call, decode the barcode"""
        raise NotImplementedError()


class DynamicAlignmentLibraryDecoder(LibraryDecoder):
    """
    Uses a dynamic programing alignment algorithm to call the barcode

    Will align to the library tag (and other known static regions) and use
    that alignment to locate the required sections to call
    """

    def __init__(
        self,
        library: DELibrary,
        alignment_algorithm: Literal["semi", "hybrid"],
        use_hamming: bool = True,
    ):
        """
        Initialize the DynamicAlignmentLibraryDecoder

        Parameters
        ----------
        library: DELibrary
            the library to decode from
        alignment_algorithm: Literal["semi", "hybrid"], default "semi"
            the type alignment algorithm to use
            semi is a semi global and hybrid is a hybrid semi global alignment
            see aligning docs for more details
        use_hamming: bool, default True
            using hamming correction
            only used for barcodes with hamming encoded parts
        """
        super().__init__(library=library, use_hamming=use_hamming)

        # assign the aligner
        self.aligner: Aligner
        if alignment_algorithm == "semi":
            self.aligner = SemiGlobalAligner()
        elif alignment_algorithm == "hybrid":
            self.aligner = HybridSemiGlobalAligner()
        else:
            raise ValueError(
                f"unknown building block alignment algorithm: '{alignment_algorithm}';"
                f"currently supported values are ['semi', 'hybrid']"
            )

        self._bb_callers: dict[str, BuildingBlockSetTagCaller] = {
            bb_section.section_name: BuildingBlockSetTagCaller(
                building_block_tag_section=bb_section,
                building_block_set=bb_set,
                use_hamming=self.use_hamming,
            )
            for bb_section, bb_set in self.library.iter_bb_barcode_sections_and_sets()
        }

        self._barcode_reference = self.library.barcode_schema.get_full_barcode()

        _barcode_section_spans = self.library.barcode_schema.get_section_spans()
        self._bb_sections_spans_to_search_for: dict[str, slice] = {
            section.section_name: _barcode_section_spans[section.section_name]
            for section in self.library.barcode_schema.building_block_sections
        }
        self._umi_span: slice | None = _barcode_section_spans.get("umi", None)

        # check alignment of callers and sections:
        _missing_sections = set(self._bb_callers.keys()) - set(
            self._bb_sections_spans_to_search_for.keys()
        )
        if len(_missing_sections) != 0:
            raise ValueError(
                f"building block caller for section(s) '{_missing_sections}'"
                f"are missing from the library barcode"
            )

    def call_barcode(self, library_call: ValidLibraryCall) -> DecodedBarcode:
        """
        Given a library call, decode the read

        Parameters
        ----------
        library_call: ValidLibraryCall
            the library call to decode

        Returns
        -------
        DecodedBarcode
        """
        _alignment = self.aligner.align(self._barcode_reference, library_call.sequence.sequence)

        _barcode_section_alignment: dict[str, str] = dict()
        for section_name, section_span in self._bb_sections_spans_to_search_for.items():
            _aligned_span = _alignment.map_seq1_span_to_seq2_span(
                section_span.start, section_span.stop
            )
            _barcode_section_alignment[section_name] = library_call.sequence.sequence[
                slice(*_aligned_span)
            ]

        # call the building blocks
        bb_calls: dict[str, BuildingBlockCall] = dict()
        for bb_section_name, bb_caller in self._bb_callers.items():
            bb_calls[bb_section_name] = bb_caller.call_building_block(
                _barcode_section_alignment[bb_section_name]
            )

        # extract the UMI section if there is one
        _umi: UMI | None = None
        if self._umi_span is not None:
            _umi = UMI(
                library_call.sequence.sequence[
                    slice(
                        *_alignment.map_seq1_span_to_seq2_span(
                            self._umi_span.start, self._umi_span.stop
                        )
                    )
                ]
            )

        return DecodedBarcode(library=library_call, building_blocks=bb_calls, umi=_umi)


#####
# Calling by regex alignment has be removed, but these functions
# will stick around incase we want to roll that back
#
#
# def _between(
#         val: int,
#         start: int,
#         stop: int,
#         right_inclusive: bool = True,
#         left_inclusive: bool = False):
#     """
#     Helper function for determining if some value is between two values
#
#     Parameters
#     ----------
#     val: int
#         value to check if in between
#     start: int
#         start of the between range
#     stop: int
#         stop of the between range
#     right_inclusive: bool = True
#         include the `start` value in the range
#     left_inclusive: bool = True
#         include the `stop` value in the range
#
#     Returns
#     -------
#     bool
#         if `val` is between `start` and `stop`
#     """
#     return (start - right_inclusive) < val < (stop + left_inclusive)
#
#
# def _match_call_mode(
#         match: FullBarcodeMatch,
#         barcode_schema: BarcodeSchema
# ) -> BarcodeAlignment:
#     """
#     Call mode function used to generate section alignment in match mode
#
#     Notes
#     -----
#     Needs FullBarcodeMatches (generated with the `full_match=True` flag
#     during the matching)
#
#     Parameters
#     ----------
#     match: FullBarcodeMatch
#         the match to generate the alignment for
#     barcode_schema:
#         the schema defining the barcode to align to
#
#     Returns
#     -------
#     BarcodeAlignment
#     """
#     _tmp_alignment = barcode_schema.barcode_spans.copy()
#
#     _errors: List[Tuple[int, str]] =
#     list(zip(*match.get_sorted_match_sequence_error_idx(exclude_substitute=True)))
#     if len(_errors) == 0:
#         return BarcodeAlignment(_tmp_alignment)
#
#     (_err_idx, _err_type) = _errors.pop(0)
#     _buff = 0
#     _spans = barcode_schema.barcode_spans.copy()
#     for _key, (_start, _stop) in _tmp_alignment.items():
#         _start = _start + _buff
#         _stop = _stop + _buff
#         _end = _stop
#         while _err_idx < _end:
#             if _between(_err_idx, _start, _stop):
#                 if _err_type == "insert":
#                     _stop += 1
#                 elif _err_type == "delete":
#                     _stop -= 1
#                 else:
#                     raise ValueError(f"unrecognized error type {_err_type}")
#
#                 try:
#                     (_err_idx, _err_type) = _errors.pop(0)
#                 except IndexError:
#                     (_err_idx, _err_type) = (1e9, "no-more-errors")
#         _spans[_key] = (_start, _stop)
#         _buff += (_stop - _end)
#     return BarcodeAlignment(_spans)
####
