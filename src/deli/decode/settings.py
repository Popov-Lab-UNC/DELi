"""settings for decoding experiments"""

from typing import Literal

from deli.base import BaseSettings


class DecodingSettings(BaseSettings):
    """
    Define parameters for decoding experiments

    Only parameters relating the algorithm should be here
    Setting relating to IO should be handled outside this context
    (like in the click command)
    """

    def __init__(
        self,
        library_error_tolerance: float = 0.2,
        alignment_algorithm: Literal["semi", "hybrid"] = "semi",
        bb_calling_approach: Literal["alignment"] = "alignment",
        use_hamming: bool = True,
        track_statistics: bool = True,
    ):
        """
        Initialize the decoder settings

        Notes
        -----
        More details about the exact effect of these settings can
        be found in the "Decoding" docs

        Parameters
        ----------
        library_error_tolerance: float, default = 0.2
            the percent error to be tolerated in the library section
            this will be converted to number of errors based on tag size
            and down to the nearest whole number
            for example, a library with 14 nucleotides would tolerate
            1, 2, and 4 errors for an error tolerance of 0.1, 0.2 and 0.3 respectively
        alignment_algorithm: Literal["semi", "hybrid"], default = "semi"
            the algorithm to use for alignment
            only used if bb_calling_approach is "alignment"
        bb_calling_approach: Literal["alignment"], default = "alignment"
            the algorithm to use for bb_calling
            right now only "alignment" mode is supported
        use_hamming: bool, default = True
            enable (`True`) or disable (`False`) hamming decoding
            only used if a library specifies tags as hamming encoded
            Note: if hamming encoded libraries are given, and `use_hamming` is
            `False`, the hamming decoding will not occur, even though it is possible
        track_statistics: bool, default = True
            track over statistic during decoding
            statistics include, for example, number or seq that failed
            library demultiplexing, failed bb looked, which cycles failed etc.
            Full details on all statistics can be found in the "Decoding" docs
        """
        super().__init__(
            library_error_tolerance=library_error_tolerance,
            alignment_algorithm=alignment_algorithm,
            bb_calling_approach=bb_calling_approach,
            use_hamming=use_hamming,
            track_statistics=track_statistics,
        )
