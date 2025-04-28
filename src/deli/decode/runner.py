"""code for running a decoding experiment"""

import logging
import warnings
from collections import defaultdict
from copy import deepcopy
from functools import partial
from os import PathLike, getpid, makedirs, path
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
        _name = f"DELi-DECODER-{getpid()}"
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

    def run(self, use_tqdm: bool = False):
        """
        Run the decoder

        Parameters
        ----------
        use_tqdm: bool, default = False
            turn on a tqdm tracking bar
            only recommended if running a single
            runner
        """
        for selection in self.experiment.selections:
            self.logger.info(f"Running decoding for selection: {selection.selection_id}")
            selection_degen_counter = self._degen()
            selection_decoder = deepcopy(self.decoder)
            selection_statistics = selection_decoder.decode_statistics
            for seq_record in tqdm(
                selection.get_sequence_reader(),
                desc=f"Running Decoding for selection {selection.selection_id}",
                disable=not use_tqdm,
            ):
                selection_statistics.num_seqs_read += 1

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
            self.selection_info.append((selection, selection_statistics, selection_degen_counter))
            self.logger.info(f"Completed decoding for selection: {selection.selection_id}")

    def write_output(
        self,
        decode_out_dir: str | PathLike = "./decode_results",
        decode_output_suffix: str = "decode_output",
        decode_statistics_suffix: str = "decode_statistics",
        decode_report_suffix: str = "decode_report",
    ):
        """
        Write the decoding output files to disk

        If two selections in the experiment have the same ID
        (which is bad practice and should be avoided), "_2" will
        be appended to the selection ID to ensure unique file names.
        If more duplicates appear, they will be numbered sequentially
        e.g. "_3" is next, then "_4"

        Notes
        -----
        Decode output will be a CSV with the columns
        DEL_ID, SMILES, RAW_COUNT, and UMI_CORRECTED_COUNT

        If SMILES are not available for the DEL, the SMILES column will be empty.
        If DELs lack a UMI, the UMI_CORRECTED_COUNT will be the same as the RAW_COUNT

        Parameters
        ----------
        decode_out_dir: str | PathLike = "./decode_results",
            path to directory to save results to
        decode_output_suffix: str = "decode_output",
            suffix to attached to selection IDs
        decode_statistics_suffix: str = "decode_statistics",
        decode_report_suffix: str = "decode_report",
        """
        makedirs(decode_out_dir, exist_ok=True)

        existing_ids: defaultdict = defaultdict(int)

        for selection, selection_stats, selection_degen in self.selection_info:
            # check for duplicate selection IDs
            _selection_id: str
            if selection.selection_id in existing_ids:
                _selection_id = (
                    f"{selection.selection_id}_{existing_ids[selection.selection_id] + 1}"
                )
                warnings.warn(
                    f"Duplicate selection ID found: {selection.selection_id}. "
                    f"Saving as duplicate {_selection_id}.",
                    stacklevel=2,
                )
            else:
                _selection_id = selection.selection_id

            existing_ids[selection.selection_id] += 1

            self.logger.info(f"Writing decode results for selection: {_selection_id}")
            _decode_out_path = path.join(
                decode_out_dir, f"{_selection_id}_{decode_output_suffix}.csv"
            )
            _decode_statistics_out_path = path.join(
                decode_out_dir, f"{_selection_id}_{decode_statistics_suffix}.json"
            )
            _decode_report_out_path = path.join(
                decode_out_dir, f"{_selection_id}_{decode_report_suffix}.html"
            )

            selection_degen.to_file(_decode_out_path)
            self.logger.debug(
                f"Wrote decode results for selection " f"{_selection_id} to {_decode_out_path}"
            )
            selection_stats.to_file(_decode_statistics_out_path)
            self.logger.debug(
                f"Wrote decode statistics for selection "
                f"{_selection_id} to {_decode_statistics_out_path}"
            )
            build_decoding_report(
                self.experiment, selection, selection_stats, _decode_report_out_path
            )
            self.logger.debug(
                f"Wrote decode report for selection "
                f"{_selection_id} to {_decode_report_out_path}"
            )
