"""code for program runners"""

import abc
import logging
import os
import warnings

from tqdm import tqdm

from deli.decode import (
    DecodedBarcode,
    DELibraryPoolCounter,
    DELibraryPoolIdCounter,
    DELibraryPoolIdUmiCounter,
    DELPoolDecoder,
)
from deli.dna import SequenceReader

from ._logging import get_dummy_logger, get_logger
from .selection import DecodingExperiment
from .statistics import DecodeStatistics


class BaseRunner(abc.ABC):
    """
    Base class for all runners

    Runners should be implemented for any complex process
    that is going to be given a CLI process.
    This will help keep the CLI script clean(er)
    """

    _runner_name: str

    def __init__(self, disable_logging: bool = False, debug: bool = False):
        """Initialize the runner"""
        _name = f"{self._runner_name}-{os.getpid()}"
        self.logger: logging.Logger = (
            get_dummy_logger() if disable_logging else get_logger(_name, debug=debug)
        )

    @abc.abstractmethod
    def run(self, *args, **kwargs):
        """Run the runner"""
        raise NotImplementedError()


class DecodingExperimentRunner(BaseRunner):
    """
    Main runner for decoding experiments

    For most users, this in the main entry point for doing decoding.

    Notes
    -----
    The runner will write logs tracking progress
    and other warnings.
    Logger will be named after the runner PID
    """

    def __init__(
        self,
        decode_experiment: DecodingExperiment,
        disable_logging: bool = False,
        debug: bool = False,
    ):
        """
        Initialize the DecodingExperimentRunner

        Parameters
        ----------
        decode_experiment: DecodingExperiment
            the experiment to decode
        disable_logging: bool, default = False
            if true, will turn off logging
        debug: bool, default = False
            if true, will enable debug logging
        """
        super(DecodingExperimentRunner, self).__init__(
            disable_logging=disable_logging, debug=debug
        )

        self.experiment = decode_experiment
        self.statistics = DecodeStatistics()  # init and pass to all sub-objects

        # initialize all the decoding object required
        self.decoder = DELPoolDecoder(
            library_pool=self.experiment.library_pool,
            decode_statistics=self.statistics,
            **self.experiment.decode_settings.__dict__,  # unpack the settings
        )

        _has_umi = all(
            [lib.barcode_schema.has_umi() for lib in self.experiment.library_pool.libraries]
        )

        # TODO right now all barcodes must have a UMI to enable, maybe should not be this
        #  will throw warning and ask user to raise issue to see if that ever happens
        if any(
            [lib.barcode_schema.has_umi() for lib in self.experiment.library_pool.libraries]
        ) and (not _has_umi):
            warnings.warn(
                "DELi does not support UMI degeneration for library collections with "
                "only some DELs having a UMI region; if this is an issue for your DEL "
                "please raise an issue and ask for this feature to be added.",
                stacklevel=2,
            )

        self.degen: DELibraryPoolCounter
        if _has_umi:
            self.degen = DELibraryPoolIdUmiCounter(
                umi_clustering=decode_experiment.decode_settings.__getattribute__(
                    "umi_clustering"
                ),
                min_umi_cluster_dist=decode_experiment.decode_settings.__getattribute__(
                    "umi_min_distance"
                ),
            )
        else:
            self.degen = DELibraryPoolIdCounter()

    def run(self, sequence_reader: SequenceReader, use_tqdm: bool = False):
        """
        Run the decoder

        Parameters
        ----------
        sequence_reader: SequenceReader
            the sequence reader with sequences to decode
        use_tqdm: bool, default = False
            turn on a tqdm tracking bar
            only recommended if running a single
            runner
        """
        for seq_record in tqdm(sequence_reader, desc="Running Decoding", disable=not use_tqdm):
            self.statistics.num_seqs_read += 1

            decoded_barcode = self.decoder.decode_read(seq_record)

            # skip failed reads
            if isinstance(decoded_barcode, DecodedBarcode):
                self.statistics.num_seqs_decoded_per_lib[decoded_barcode.library.library_id] += 1
                if self.degen.count_barcode(decoded_barcode):
                    # only up the degen count if not a degenerate read
                    self.statistics.num_seqs_degen_per_lib[decoded_barcode.library.library_id] += 1

    def write_output(
        self,
        decode_out_path: str | os.PathLike = "./decode_out.csv",
        decode_statistics_out_path: str | os.PathLike = "./decode_statistics.json",
    ):
        """
        Write the decoding output files to disk

        Notes
        -----
        Decode output will be a CSV with the columns
        DEL_ID, SMILES, RAW_COUNT, and UMI_CORRECTED_COUNT

        If SMILES are not available for the DEL, the SMILES column will be empty.
        If DELs lack a UMI, the UMI_CORRECTED_COUNT will be the same as the RAW_COUNT

        Parameters
        ----------
        decode_out_path
        decode_statistics_out_path
        """
        self.degen.to_file(decode_out_path, "csv")
        self.statistics.to_file(decode_statistics_out_path)
