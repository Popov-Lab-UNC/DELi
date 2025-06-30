"""code for library calling"""

import abc
from typing import Optional, TypeAlias, Union

from cutadapt.adapters import BackAdapter, FrontAdapter, MultipleAdapters, SingleMatch
from dnaio import SequenceRecord

from deli.dels import DELibrary, DELibraryCollection

from .calls import FailedCall, ValidCall
from .demultiplex import Demultiplexer


class ValidLibraryCall(ValidCall):
    """base class for all valid library call objects"""

    library: DELibrary
    sequence: SequenceRecord


class FailedLibraryCall(FailedCall):
    """base class for all failed library calls"""

    pass


# type hing for library calls
LibraryCall: TypeAlias = Union[ValidLibraryCall, FailedLibraryCall]


class SingleLibraryCall(ValidLibraryCall):
    """
    Holds information about a single library call

    Single library calls are from single reads (not paired reads)
    """

    def __init__(self, library: DELibrary, sequence: SequenceRecord):
        """
        Initialize the SingleLibraryCall

        Parameters
        ----------
        library: DELibrary
            the called library
        sequence: SequenceRecord
            the observed_seq record trimmed to the just the matching region
        """
        self.library = library
        self.sequence = sequence
        self._success = True


class PairedLibraryCall(ValidLibraryCall):
    """
    Holds information about a paired library call

    Paired library calls are from paired reads
    (meaning the forward and reverse sequences for that read are provided)
    """

    def __init__(self, library: DELibrary, sequence: SequenceRecord, sequence2: SequenceRecord):
        """
        Initialize the PairedLibraryCall

        Parameters
        ----------
        library: DELibrary
            the called library
        sequence: SequenceRecord
            the observed_seq record trimmed to the just the matching region for fwd read
        sequence2: SequenceRecord
            the observed_seq record trimmed to the just the matching region for rev read
        """
        self.library = library
        self.sequence = sequence
        self.sequence2 = sequence2
        self._success = True


class PairedLibraryCallAmbiguous(FailedLibraryCall):
    """
    Library call for when a paired read has an ambiguous library call

    For example, if the fwd read match LIB001 and the rev matched LIB004
    """

    def __init__(
        self,
        library_1: DELibrary,
        library_2: DELibrary,
        sequence_1: SequenceRecord,
        sequence_2: SequenceRecord,
    ):
        """
        Initialize the PairedLibraryCallAmbiguous

        Parameters
        ----------
        library_1: DELibrary
            the called library for the fwd read
        library_2: DELibrary
            the called library for the rev read
        sequence_1: SequenceRecord
            the observed_seq record trimmed to the just the matching region for fwd read
        sequence_2: SequenceRecord
            the observed_seq record trimmed to the just the matching region for rev read
        """
        self.library_1 = library_1
        self.library_2 = library_2
        self.sequence_1 = sequence_1
        self.sequence_2 = sequence_2


class LibraryCallTooShort(FailedLibraryCall):
    """
    Library call for when a trimmed read is not long enough to contain all info

    For example, if a library would need a match with at least 50 nucleotides to
    be able to call the building blocks and umi, and the trimmed match is only 35,
    this is the call that should be made (not FailedLibraryCall)
    """

    def __init__(self, library: DELibrary, sequence: SequenceRecord):
        """
        Initialize the LibraryCallTooShort

        Parameters
        ----------
        library: DELibrary
            the called library for the read
        sequence: SequenceRecord
            the observed_seq record trimmed to the just the matching region of the read
        """
        self.library = library
        self.sequence = sequence


class LibraryCaller(Demultiplexer, abc.ABC):
    """
    Base class for all library callers

    Will build the required demultiplexing adapters from the
    library collection, one for each DEL in it.
    """

    def __init__(
        self,
        library_collection: DELibraryCollection,
        error_rate: Union[float, int],
        min_overlap: Optional[int] = None,
    ):
        """
        Initialize the LibraryCaller

        Parameters
        ----------
        library_collection: DELibraryCollection
            the library collection to demultiplex on
        error_rate: Union[float, int]
            the max error rate for any match
            if int, treat as raw number of errors
            if float between (0,1) treat as error rate
        min_overlap: Optional[int]
            the minimum overlap required for the library call to be valid
            low values will likely cause false positive
            best to keep as high as possible
            cannot be larger than the smallest adapter size
            if `None`, will default to a full match
        """
        self.library_collection = library_collection
        self.error_rate = error_rate
        self.min_overlap = min_overlap

        _adapters = MultipleAdapters(
            [
                FrontAdapter(
                    sequence=lib.barcode_schema.library_section.get_dna_sequence(),
                    max_errors=error_rate,
                    min_overlap=min_overlap
                    if min_overlap is not None
                    else len(lib.barcode_schema.library_section),
                    name=lib.library_id,
                )
                if lib.barcode_schema.is_library_tag_in_front()
                else BackAdapter(
                    sequence=lib.barcode_schema.library_section.get_dna_sequence(),
                    max_errors=error_rate,
                    min_overlap=min_overlap
                    if min_overlap is not None
                    else len(lib.barcode_schema.library_section),
                    name=lib.library_id,
                )
                for lib in library_collection.libraries
            ]
        )

        super().__init__(_adapters)

    def _get_trim_region(
        self, sequence: SequenceRecord, match: SingleMatch
    ) -> tuple[slice, DELibrary]:
        """
        Given a match and observed_seq, figure out the library call and return the trim slice

        Parameters
        ----------
        sequence: SequenceRecord
            the observed_seq record to get trim for
        match: SingleMatch
            the match for the observed_seq

        Returns
        -------
        tuple[slice, DELibrary]
            the trim slice and the called DEL
        """
        called_lib = self.library_collection.get_library(match.adapter.name)
        # a padding of 20 is added to the trim region to allow some leway later alignment
        retain_region = slice(
            max(0, match.rstart - called_lib.barcode_schema.get_length_before_library() - 20),
            min(
                len(sequence),
                match.rstop + called_lib.barcode_schema.get_length_after_library() + 20,
            ),
        )
        return retain_region, called_lib

    @abc.abstractmethod
    def call_library(self, *args, **kwargs) -> LibraryCall:
        """Determine the library call for a given input"""
        raise NotImplementedError()


class SingleReadLibraryCaller(LibraryCaller):
    """
    Library calling for a single reads

    Notes
    -----
    If your reads have the possibility of being duplexed
    (for example using the nanopore)
    you should use the `SingleReadDuplexLibraryCaller` instead
    """

    def __init__(
        self,
        library_collection: DELibraryCollection,
        error_rate: Union[float, int] = 0.1,
        min_overlap: Optional[int] = None,
        revcomp: bool = False,
    ):
        """
        Initialize the SingleReadLibraryCaller

        Parameters
        ----------
        library_collection: DELibraryCollection
            the library collection to demultiplex on
        error_rate: Union[float, int]
            the max error rate for any match
            if int, treat as raw number of errors
            if float between (0,1) treat as error rate
        min_overlap: Optional[int]
            the minimum overlap required for the library call to be valid
            low values will likely cause false positive
            best to keep as high as possible
            cannot be larger than the smallest adapter size
            if `None`, will default to a full match
        revcomp: bool
            search for a match in the reverse compliment as well
            use only if you do not know which direction the
            reads are coming in as.
        """
        super().__init__(library_collection, error_rate, min_overlap)
        self.revcomp = revcomp

    def _pick_seq_and_match(
        self, sequence: SequenceRecord
    ) -> tuple[SequenceRecord, Optional[SingleMatch]]:
        """
        Determines the best observed_seq transformation and its match for a given read

        Notes
        -----
        Right now the only observed_seq modification that exists is `revcomp`

        Parameters
        ----------
        sequence: SequenceRecord
            the read to match to

        Returns
        -------
        tuple[SequenceRecord, Optional[SingleMatch]]
            the best observed_seq transformation and its match
            match is `None` if no library match was found
        """
        match = self._demultiplex(sequence)

        if self.revcomp:
            revcomp_seq: SequenceRecord = sequence.reverse_complement()
            revcomp_match = self._demultiplex(revcomp_seq)
            return self._compare_matches((sequence, match), (revcomp_seq, revcomp_match))
        else:
            return sequence, match

    def call_library(self, sequence: SequenceRecord) -> LibraryCall:
        """
        Call the library for this read

        Will return a `LibraryCall` object.
        This could be a successful call *or* a failed call.
        This function does not remove failures.

        Parameters
        ----------
        sequence: SequenceRecord
            the read to match to

        Returns
        -------
        LibraryCall
            the library call
        """
        sequence, match = self._pick_seq_and_match(sequence)

        if match is None:
            return FailedLibraryCall()

        trim_region, called_lib = self._get_trim_region(sequence, match)
        trimmed_seq = sequence[trim_region]

        # check for too short
        if len(trimmed_seq) < called_lib.barcode_schema.min_length:
            return LibraryCallTooShort(sequence=trimmed_seq, library=called_lib)

        return SingleLibraryCall(library=called_lib, sequence=trimmed_seq)


class SingleReadDuplexLibraryCaller(SingleReadLibraryCaller):
    """
    Library calling for a single reads that might be duplexed

    Notes
    -----
    Duplexing is when a read is not paired, but the rev read is read in
    by the sequencer right after the fwd read, making it appear as a single
    continous read. Nanopore, for example, is prone to doing this.
    In this case, a more accurate match can be made by looking for
    possible duplexes and converting them to paired calls.

    Warning: should probably never use outside of having nanopore data.
    The original authors had only nanopore data, so hence why this niche
    class the really no one should every use exists
    """

    def __init__(
        self,
        library_collection: DELibraryCollection,
        error_rate: Union[float, int] = 0.1,
        min_overlap: int = 3,
        revcomp: bool = False,
        fail_on_ambiguous: bool = True,
    ):
        """
        Initialize the SingleReadDuplexLibraryCaller

        Parameters
        ----------
        library_collection: DELibraryCollection
            the library collection to demultiplex on
        error_rate: Union[float, int]
            the max error rate for any match
            if int, treat as raw number of errors
            if float between (0,1) treat as error rate
        min_overlap: Optional[int]
            the minimum overlap required for the library call to be valid
            low values will likely cause false positive
            best to keep as high as possible
            cannot be larger than the smallest adapter size
            if `None`, will default to a full match
        revcomp: bool
            search for a match in the reverse compliment as well
            use only if you do not know which direction the
            reads are coming in as.
        fail_on_ambiguous: bool = True
            the the duplex rev read matches a different
            library then the fwd read, return a failed call
            instead.
            If false, will just treat return the fwd call as
            a single match if an ambiguous library call is found
        """
        super().__init__(library_collection, error_rate, min_overlap, revcomp=revcomp)
        self.fail_on_ambiguous = fail_on_ambiguous

    def call_library(self, sequence: SequenceRecord) -> LibraryCall:
        """
        Call the library for this read

        Will return a `LibraryCall` object.
        This could be a successful call *or* a failed call.
        This function does not remove failures.

        Parameters
        ----------
        sequence: SequenceRecord
            the read to match to

        Returns
        -------
        LibraryCall
            the library call
        """
        sequence, match = self._pick_seq_and_match(sequence)

        # fail if cannot find any match
        if match is None:
            return FailedLibraryCall()

        trim_region, called_lib = self._get_trim_region(sequence, match)
        trimmed_seq = sequence[trim_region]
        remaining_seq = sequence[trim_region.stop :]

        # check if read has possibility of being duplex
        if len(remaining_seq) >= called_lib.barcode_schema.min_length:
            # duplexes are always the revcomp
            remaining_seq = remaining_seq.reverse_complement()

            # search the duplex
            duplex_match = self.adapters.match_to(remaining_seq.sequence)

            if duplex_match is not None:
                duplex_trim_region, called_duplex_lib = self._get_trim_region(
                    remaining_seq, duplex_match
                )
                trimmed_duplex_seq = remaining_seq[duplex_trim_region]

                if called_duplex_lib.library_id != called_lib.library_id:
                    if self.fail_on_ambiguous:
                        return PairedLibraryCallAmbiguous(
                            library_1=called_lib,
                            library_2=called_duplex_lib,
                            sequence_1=trimmed_seq,
                            sequence_2=trimmed_duplex_seq,
                        )
                    else:
                        return SingleLibraryCall(library=called_duplex_lib, sequence=trimmed_seq)
                else:
                    return PairedLibraryCall(
                        library=called_lib, sequence=trimmed_seq, sequence2=trimmed_duplex_seq
                    )
        # if no duplex found, call like a simplex read
        # check for too short
        if len(trimmed_seq) < called_lib.barcode_schema.min_length:
            return LibraryCallTooShort(sequence=trimmed_seq, library=called_lib)

        return SingleLibraryCall(library=called_lib, sequence=trimmed_seq)
