"""code for running a decoding experiment"""

import logging
import os
import warnings
from copy import deepcopy
from functools import partial
from typing import Callable

from tqdm import tqdm

from deli._logging import get_dummy_logger, get_logger
from deli.dels import SequencedSelection

from .decoder import DecodedBarcode, DecodeStatistics, DELPoolDecoder
from .degen import DELibraryPoolCounter, DELibraryPoolIdCounter, DELibraryPoolIdUmiCounter
from .experiment import DecodingExperiment
from .report import build_decoding_report


class DecodingExperimentRunner:
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
        _name = f"DELi-DECODER-{os.getpid()}"
        self.logger: logging.Logger = (
            get_dummy_logger() if disable_logging else get_logger(_name, debug=debug)
        )

        self.experiment = decode_experiment

        # initialize all the decoding object required
        self.decoder = DELPoolDecoder(
            library_pool=self.experiment.library_pool,
            decode_statistics=DecodeStatistics(),
            library_error_tolerance=self.experiment.decode_settings.get(
                "library_error_tolerance", 0.1
            ),
            min_library_overlap=self.experiment.decode_settings.get("min_library_overlap", 10),
            alignment_algorithm=self.experiment.decode_settings.get("alignment_algorithm", "semi"),
            bb_calling_approach=self.experiment.decode_settings.get("bb_calling_approach", "bio"),
            revcomp=self.experiment.decode_settings.get("revcomp", False),
            max_read_length=self.experiment.decode_settings.get("max_read_length", None),
            min_read_length=self.experiment.decode_settings.get("min_read_length", None),
            read_type=self.experiment.decode_settings.get("read_type", "single"),
            use_hamming=self.experiment.decode_settings.get("use_hamming", True),
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

        self._degen: Callable
        if _has_umi:
            self._degen = partial(
                DELibraryPoolIdUmiCounter,
                umi_clustering=decode_experiment.decode_settings.get("umi_clustering", False),
                min_umi_cluster_dist=decode_experiment.decode_settings.get("umi_min_distance", 2),
            )
        else:
            self._degen = DELibraryPoolIdCounter

        self.selection_info: list[
            tuple[SequencedSelection, DecodeStatistics, DELibraryPoolCounter]
        ] = list()

    def run(self, save_failed_to: str | os.PathLike | None = None, use_tqdm: bool = False):
        """
        Run the decoder

        Parameters
        ----------
        save_failed_to: str | PathLike | None, default = None
            if provided, will save failed reads to this file
            will include the read_id, the sequence, the quality chain,
            and reason failed
        use_tqdm: bool, default = False
            turn on a tqdm tracking bar
            only recommended if running a single
            runner
        """
        # loop through all selections
        for selection in self.experiment.selections:
            # initalization of selection specific objects
            self.logger.info(f"Running decoding for selection: {selection.selection_id}")
            selection_degen_counter = self._degen()
            selection_decoder = deepcopy(self.decoder)
            selection_statistics = selection_decoder.decode_statistics

            # open failed reads file if specified
            self.logger.info(
                f"Saving failed reads to {save_failed_to} for selection {selection.selection_id}"
            )
            _failed_file = (
                os.path.join(save_failed_to, f"{selection.selection_id}_decode_failed.csv")
                if save_failed_to is not None
                else None
            )
            fail_csv_file = open(_failed_file, "w") if _failed_file is not None else None

            # write header to the failed reads CSV
            if fail_csv_file is not None:
                fail_csv_file.write("read_id,sequence,quality,reason_failed,lib_call\n")

            # look through all sequences in the selection
            for i, seq_record in enumerate(
                tqdm(
                    selection.get_sequence_reader(),
                    desc=f"Running Decoding for selection {selection.selection_id}",
                    disable=not use_tqdm,
                )
            ):
                selection_statistics.num_seqs_read += 1

                # decode the read
                decoded_barcode = selection_decoder.decode_read(seq_record)

                # skip failed reads
                if isinstance(decoded_barcode, DecodedBarcode):
                    selection_statistics.num_seqs_decoded_per_lib[
                        decoded_barcode.library.library_id
                    ] += 1
                    if selection_degen_counter.count_barcode(decoded_barcode):
                        # only up the degen count if not a degenerate read
                        selection_statistics.num_seqs_degen_per_lib[
                            decoded_barcode.library.library_id
                        ] += 1
                elif fail_csv_file is not None:
                    # if we are saving failed reads, save the read
                    fail_csv_file.write(
                        f"{seq_record.id},"
                        f"{seq_record.sequence},"
                        f"{seq_record.qualities},"
                        f"{decoded_barcode.__class__.__name__}\n"
                    )

                if ((i + 1) % 1000) == 0:
                    self.logger.debug(
                        f"Decoded {i+1} reads for selection {selection.selection_id}"
                    )

            if fail_csv_file is not None:
                fail_csv_file.close()
                self.logger.debug(
                    f"Saved failed reads to {save_failed_to} "
                    f"for selection {selection.selection_id}"
                )

            self.selection_info.append((selection, selection_statistics, selection_degen_counter))
            self.logger.info(f"Completed decoding for selection: {selection.selection_id}")

    def write_decode_report(self, out_dir: str | os.PathLike, prefix: str = ""):
        """
        Write the decoding report for this run

        Notes
        -----
        Will write a unique report for each selection in the experiment.
        All files will be `html` and be named "{?prefix}_<selection_id>_decode_report.html"


        Parameters
        ----------
        out_dir: str | PathLike
            path to directory to save results to
            will create the directory if it does not exist
        prefix: str, default = ""
            prefix for output files
        """
        os.makedirs(out_dir, exist_ok=True)
        for selection, selection_stats, _ in self.selection_info:
            _filename = (
                f"{prefix}_{selection.selection_id}_decode_report.html"
                if prefix
                else f"{selection.selection_id}_decode_report.html"
            )
            _out_path = os.path.join(out_dir, _filename)
            build_decoding_report(self.experiment, selection, selection_stats, _out_path)
            self.logger.debug(
                f"Wrote decode report for selection " f"{selection.selection_id} to {_out_path}"
            )

    def write_decode_statistics(self, out_dir: str | os.PathLike, prefix: str = ""):
        """
        Write the decoding statistics for this run

        Notes
        -----
        Will write a unique statistics file for each selection in the experiment.
        All files will be `json` and be named "{?prefix}_<selection_id>_decode_statistics.json"

        Parameters
        ----------
        out_dir: str | PathLike
            path to directory to save results to
            will create the directory if it does not exist
        prefix: str, default = ""
            prefix for output files
        """
        os.makedirs(out_dir, exist_ok=True)
        for selection, selection_stats, _ in self.selection_info:
            _filename = (
                f"{prefix}_{selection.selection_id}_decode_statistics.json"
                if prefix
                else f"{selection.selection_id}_decode_statistics.json"
            )
            _out_path = os.path.join(out_dir, _filename)
            selection_stats.to_file(_out_path)
            self.logger.debug(
                f"Wrote decode statistics for selection "
                f"{selection.selection_id} to {_out_path}"
            )

    def write_decode_results(self, out_dir: str | os.PathLike, prefix: str = ""):
        """
        Write the decoding results for this run

        Notes
        -----
        Will write a unique result file for each selection in the experiment.
        All files will be `csv` and be named "{?prefix}_<selection_id>_decode_results.csv"

        Parameters
        ----------
        out_dir: str | PathLike
            path to directory to save results to
            will create the directory if it does not exist
        prefix: str, default = ""
            prefix for output files
        """
        os.makedirs(out_dir, exist_ok=True)
        for selection, _, selection_degen in self.selection_info:
            _filename = (
                f"{prefix}_{selection.selection_id}_decode_results.csv"
                if prefix
                else f"{selection.selection_id}_decode_results.csv"
            )
            _out_path = os.path.join(out_dir, _filename)
            selection_degen.to_file(_out_path)
            self.logger.debug(
                f"Wrote decode results for selection " f"{selection.selection_id} to {_out_path}"
            )
