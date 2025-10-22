"""code for running a decoding experiment"""

import datetime
import logging
import os
import random
import warnings
from itertools import zip_longest
from os import PathLike
from typing import Any, Literal
from typing_extensions import Self

import yaml
from tqdm import tqdm

from deli._logging import get_dummy_logger, get_logger
from deli.dels import DELibrary, DELibraryCollection, Selection, SequencedSelection
from deli.dels.library import EnumerationRunError

from .decoder import DecodedDELCompound, DecodeStatistics, DELCollectionDecoder
from .degen import DELCollectionCounter, DELCollectionIdCounter, DELCollectionIdUmiCounter
from .report import build_decoding_report


class DecodingRunParsingError(Exception):
    """Exception to raise when a decoding experiment file is invalid"""

    pass


class DecodingSettings(dict):
    """
    Define parameters for decoding experiments

    Only parameters relating to the algorithm should be here
    Setting relating to IO should be handled outside this context
    (like in the click command)
    """

    def __init__(
        self,
        library_error_tolerance: float = 0.1,
        min_library_overlap: int | None = 10,
        alignment_algorithm: Literal["semi", "hybrid"] = "semi",
        bb_calling_approach: Literal["alignment", "bio"] = "bio",
        revcomp: bool = False,
        max_read_length: int | None = None,
        min_read_length: int | None = None,
        read_type: Literal["single", "paired"] = "single",
        disable_error_correction: bool = False,
        umi_clustering: bool = False,
        umi_min_distance: int = 2,
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
        min_library_overlap: int or None, default = 7
            the minimum number of nucleotides required to match
            the library tag
            This is because the demultiplexing will accept truncated matches
            at the front/back of the tag. For example a tag of AGCTGGTTC
            could match a read of GTTC if the min overlap was <=4
            If `None`, will default to the exact length of the tag, meaning
            the whole tag is expected.
            The recommended value is greater than 8, as the odds of a match this strong
            to be accidental are low
        alignment_algorithm: Literal["semi", "hybrid"], default = "semi"
            the algorithm to use for alignment
            only used if bb_calling_approach is "alignment"
        read_type: Literal["single", "paired"], default = "single"
            the type of read
            paired are for paired reads
            all other read types are single
        revcomp: bool, default = False
            If true, search the reverse compliment as well
        max_read_length: int or None, default = None
            maximum length of a read to be considered for decoding
            if above the max, decoding will fail
            if `None` will default to 5x the min_read_length
        min_read_length: int or None, default = None
            minimum length of a read to be considered for decoding
            if below the min, decoding will fail
            if `None` will default to smallest min match length of
            any library in the collection considered for decoding
        bb_calling_approach: Literal["alignment"], default = "alignment"
            the algorithm to use for bb_calling
            right now only "alignment" mode is supported
        disable_error_correction: bool, default = False
            disable error correction for barcode sections
            that can be error corrected
        umi_clustering: bool, default = False
            when doing degeneration, consider two similar UMIs to be the same
            similarity is based on levenshtein distance and `umi_min_distance`
        umi_min_distance: int, default = 2
            the minimum distance between two UMIs to be considered unique
            only used `umi_clustering` is `True`
        """
        super().__init__(
            library_error_tolerance=library_error_tolerance,
            min_library_overlap=min_library_overlap,
            alignment_algorithm=alignment_algorithm,
            read_type=read_type,
            revcomp=revcomp,
            max_read_length=max_read_length,
            min_read_length=min_read_length,
            bb_calling_approach=bb_calling_approach,
            disable_error_correction=disable_error_correction,
            umi_clustering=umi_clustering,
            umi_min_distance=umi_min_distance,
        )

    def to_file(self, path: str | PathLike):
        """
        Save settings to a YAML file

        Parameters
        ----------
        path : str | PathLike
            path to save settings to
        """
        yaml.dump(self.__dict__, open(path, "w"))

    @classmethod
    def from_file(cls, path: str | PathLike) -> Self:
        """
        Load settings from a yaml file

        Will first check if there is a "decode_settings" key
        and load settings from that sub dict.
        Otherwise will load from the yaml file keys

        Parameters
        ----------
        path : str | PathLike
            Path to yaml file

        Returns
        -------
        DecodingSettings

        Raises
        ------
        RuntimeError
            if valid decode settings cannot be loaded from the passed yaml file
        """
        _data = yaml.safe_load(open(path, "r"))
        if "decode_settings" not in _data:
            try:
                return cls(**yaml.safe_load(open(path, "r")))
            except Exception as e:
                raise RuntimeError(f"Failed to load decode settings from {path}") from e
        else:
            try:
                return cls(**_data["decode_settings"])
            except Exception as e:
                raise RuntimeError(f"Failed to load decode settings from {path}") from e


class DecodingRunnerResults:
    """Hold the results of a decoding run"""

    def __init__(
        self,
        selection: Selection,
        decode_statistics: DecodeStatistics | None = None,
        degen: DELCollectionCounter | None = None,
    ):
        """
        Initialize the decoding runner results

        Parameters
        ----------
        selection: Selection
            the selection that was decoded
        decode_statistics: DecodeStatistics | None, default = None
            the decode statistics for the run
        degen: DELCollectionCounter | None, default = None
            the degeneration counter for the run
        """
        self.selection = selection
        self.decode_statistics = decode_statistics
        self.degen = degen

    def write_decode_report(self, out_path: str | os.PathLike):
        """
        Write the decoding statistics for this run

        Parameters
        ----------
        out_path: str | PathLike
            path to save decoding statistics to
        """
        if self.decode_statistics is not None:
            if os.path.dirname(out_path) != "":
                os.makedirs(os.path.dirname(out_path), exist_ok=True)
            build_decoding_report(
                selection=self.selection, stats=self.decode_statistics, out_path=out_path
            )
        else:
            raise RuntimeError("Cannot write decode report, no decode statistics found")

    def write_decode_statistics(
        self, out_path: str | os.PathLike, include_read_lengths: bool = False
    ):
        """
        Write the decoding statistics for this run

        Parameters
        ----------
        out_path: str | PathLike
            path to save decoding statistics to
        include_read_lengths: bool, default = False
            if True, will include the read lengths in the file
        """
        if self.decode_statistics is not None:
            if os.path.dirname(out_path) != "":
                os.makedirs(os.path.dirname(out_path), exist_ok=True)
            self.decode_statistics.to_file(out_path, include_read_lengths=include_read_lengths)
        else:
            raise RuntimeError("Cannot write decode statistics, no decode statistics found")

    def write_cube(
        self,
        out_path: str | os.PathLike,
        include_library_id_col: bool = False,
        include_bb_id_cols: bool = False,
        include_bb_smi_cols: bool = False,
        include_raw_count_col: bool = True,
        enumerate_smiles: bool = False,
        file_format: Literal["csv", "tsv"] = "csv",
    ) -> None:
        """
        Write the resulting DEL counts to a human-readable separated value file

        Will write a unique result file for each selection in the experiment.

        Each row will be a unique DEL ID with at least 2 and at most 11 columns.
        If including BB columns, the max cycle number is the number of columns included,
        with DELs from libraries with smaller cycle counts set to 'null'.
        Below are the columns in order of appearance:
        - DEL_ID
        - LIBRARY_ID (if `include_library_id_col` is True)
        - <For [#] of bb cycles in the library>
            - BB[#]_ID (if `include_bb_id_cols` is True)
            - BB[#]_SMILES (if `include_bb_smi_cols` is True)
        - ENUMERATED_SMILES
        - RAW_COUNT (if `include_raw_count_col` is True)
        - UMI_CORRECTED_COUNT (same as raw count if not using UMI degen)

        If `enumerate_smiles` is True, the enumerated SMILES will be generated on the
        fly using the reaction data provided in the library files.
        If this is missing, the SMILES will be set to 'null'.
        Also note that on the fly enumeration will have a dramatic impact on runtime.
        For large DELs, it can be far more efficient to enumerate the SMILES
        once, save in a database, and then map the enumerated SMILES to the DEL ID
        in the CSV file with a join.

        Notes
        -----
        If running decoding in parallel, a call to `to_csv` should be made after
        all decoding is complete and the counters have been merged.
        Merging raw csv files will result is data loss and incorrect counts,
        especially if UMI degen is used.

        Parameters
        ----------
        out_path: str | PathLike
            path to save results to
            will create the directory if it does not exist
        file_format: Literal["csv", "tsv"] = "csv"
            which file format to write to
            either "csv" or "tsv"
        include_library_id_col: bool, default = False
            include the library ID in the output
        include_bb_id_cols: bool, default = False
            include the building block IDs in the output
        include_bb_smi_cols: bool, default = False
            include the building block SMILES in the output
        include_raw_count_col: bool, default = True
            include the raw count in the output
        enumerate_smiles: bool, default = False
            include the enumerated SMILES in the output
            Note: this will have a dramatic impact on runtime

        Warnings
        --------
        UserWarning
            - raised when enumeration is requested but some libraries are
            missing enumerators
            - raised when building block smiles are requested but some libraries
            are missing this info
        """
        delimiter: str
        if file_format == "tsv":
            delimiter = "\t"
        elif file_format == "csv":
            delimiter = ","
        else:
            raise ValueError(f"file format '{file_format}' not recognized")

        # check for enumeration
        if enumerate_smiles:
            if not self.selection.library_collection.all_libs_can_enumerate():
                _libraries_missing_enumerators = [
                    lib.library_id
                    for lib in self.selection.library_collection.libraries
                    if not lib.can_enumerate()
                ]
                warnings.warn(
                    f"Some libraries are missing enumerators: "
                    f"{_libraries_missing_enumerators}. "
                    "On the fly enumeration will be skipped for these libraries. "
                    "Please check the library files to ensure they have the "
                    "correct reaction data.",
                    stacklevel=2,
                )

        # check for building block smiles
        if include_bb_smi_cols:
            if not self.selection.library_collection.all_libs_have_building_block_smiles():
                _libraries_missing_bb_smi = [
                    lib.library_id
                    for lib in self.selection.library_collection.libraries
                    if not lib.building_blocks_have_smi()
                ]
                warnings.warn(
                    f"Some libraries are missing building block smiles: "
                    f"{_libraries_missing_bb_smi}. "
                    "Building block smiles will be skipped for these libraries. "
                    "Please check the building block files to ensure they have "
                    "the correct smiles.",
                    stacklevel=2,
                )

        # get max_cycle_size
        _max_cycle_size = (
            self.selection.library_collection.max_cycle_size()
            if (include_bb_id_cols or include_bb_smi_cols)
            else 0
        )

        # set header
        _header = "DEL_ID"
        if include_library_id_col:
            _header += f"{delimiter}LIBRARY_ID"
        for i in range(_max_cycle_size):
            if include_bb_id_cols:
                _header += f"{delimiter}BB{i + 1}_ID"
            if include_bb_smi_cols:
                _header += f"{delimiter}BB{i + 1}_SMILES"
        if enumerate_smiles:
            _header += f"{delimiter}SMILES"
        if include_raw_count_col:
            _header += f"{delimiter}RAW_COUNT"
        _header += f"{delimiter}UMI_CORRECTED_COUNT\n"

        if self.degen is not None:
            if os.path.dirname(out_path) != "":
                os.makedirs(os.path.dirname(out_path), exist_ok=True)

            # write the file
            with open(out_path, "w") as f:
                # write header
                f.write(_header)
                for library_count_set in self.degen.del_counter.values():
                    for decode_barcode, counts in library_count_set.items():
                        _row = f"{decode_barcode.compound_id}"
                        if include_library_id_col:
                            _row += f"{delimiter}{decode_barcode.library.library_id}"
                        for _, bb in zip_longest(
                            range(_max_cycle_size), decode_barcode.building_blocks
                        ):
                            if include_bb_id_cols:
                                _row += f"{delimiter}{bb.bb_id if bb is not None else 'null'}"
                            if include_bb_smi_cols:
                                _row += (
                                    f"{delimiter}"
                                    f"{bb.smi if (bb is not None and bb.has_smiles()) else 'null'}"
                                )
                        if enumerate_smiles:
                            try:
                                _smi = decode_barcode.enumerate().smi
                            except EnumerationRunError:
                                _smi = "null"
                            _row += f"{delimiter}{_smi}"
                        if include_raw_count_col:
                            _row += f"{delimiter}{counts.get_raw_count()}"
                        _row += f"{delimiter}{counts.get_degen_count()}\n"
        else:
            raise RuntimeError("Cannot write cube, no degeneration counter found")


class DecodingRunner:
    """
    Main runner for decoding selections

    For most users, this in the main entry point for doing decoding.

    Notes
    -----
    The runner will write logs tracking progress
    and other warnings.
    Logger will be named after the runner PID
    """

    def __init__(
        self,
        selection: SequencedSelection,
        decode_settings: DecodingSettings | None = None,
        disable_logging: bool = False,
        debug: bool = False,
        run_id: str | None = None,
    ):
        """
        Initialize the DecodingRunner

        Parameters
        ----------
        selection: SequencedSelection
            the selection to decode
        decode_settings: DecodingSettings | None, default = None
            the settings to use for decoding
            if `None`, will use the default settings
        disable_logging: bool, default = False
            if true, will turn off logging
        debug: bool, default = False
            if true, will enable debug logging
        run_id: str, default = None
            a unique id for this run
            if `None`, will default to the PID and random number
        """
        if run_id is None:
            self.run_id = f"DELi-DECODER-{os.getpid()}-" + f"{random.randint(0, 100000):06d}"
        else:
            self.run_id = run_id

        self.logger: logging.Logger = (
            get_dummy_logger() if disable_logging else get_logger(self.run_id, debug=debug)
        )

        self.selection: SequencedSelection = selection
        self.decode_settings = (
            decode_settings if decode_settings is not None else DecodingSettings()
        )

        # initialize all the decoding object required
        self.decoder = DELCollectionDecoder(
            library_collection=self.selection.library_collection,
            decode_statistics=DecodeStatistics(),
            library_error_tolerance=self.decode_settings.get("library_error_tolerance", 0.1),
            min_library_overlap=self.decode_settings.get("min_library_overlap", 10),
            alignment_algorithm=self.decode_settings.get("alignment_algorithm", "semi"),
            bb_calling_approach=self.decode_settings.get("bb_calling_approach", "bio"),
            revcomp=self.decode_settings.get("revcomp", False),
            max_read_length=self.decode_settings.get("max_read_length", None),
            min_read_length=self.decode_settings.get("min_read_length", None),
            read_type=self.decode_settings.get("read_type", "single"),
            disable_error_correction=self.decode_settings.get("disable_error_correction", False),
        )

        _has_umi = all(
            [lib.barcode_schema.has_umi() for lib in self.selection.library_collection.libraries]
        )

        # TODO right now all barcodes must have a UMI to enable, maybe should not be this
        #  will throw warning and ask user to raise issue to see if that ever happens
        if any(
            [lib.barcode_schema.has_umi() for lib in self.selection.library_collection.libraries]
        ) and (not _has_umi):
            warnings.warn(
                "DELi does not support UMI degeneration for library collections with "
                "only some DELs having a UMI region; if this is an issue for your DEL "
                "please raise an issue and ask for this feature to be added.",
                stacklevel=2,
            )

        self.degen: DELCollectionCounter
        if _has_umi:
            self.degen = DELCollectionIdUmiCounter(
                umi_clustering=self.decode_settings.get("umi_clustering", False),
                min_umi_cluster_dist=self.decode_settings.get("umi_min_distance", 2),
            )
        else:
            self.degen = DELCollectionIdCounter()

    def run(
        self, save_failed_to: str | os.PathLike | None = None, use_tqdm: bool = False
    ) -> DecodingRunnerResults:
        """
        Run the decoder

        Parameters
        ----------
        save_failed_to: str | PathLike | None, default = None
            if provided, will save failed reads to this directory
            file will be named <selection_id>_decode_failed.tsv
            will include the read_id, the observed_seq, the quality chain,
            and reason failed
        use_tqdm: bool, default = False
            turn on a tqdm tracking bar
            only recommended if running a single
            runner

        Returns
        -------
        DecodingRunnerResults
            the results of the decode run
        """
        self.logger.info(f"Running decoding for selection: {self.selection.selection_id}")

        # open failed reads file if specified
        self.logger.info(
            f"Saving failed reads to {save_failed_to} for selection {self.selection.selection_id}"
        )
        _failed_file = (
            os.path.join(save_failed_to, f"{self.selection.selection_id}_decode_failed.tsv")
            if save_failed_to is not None
            else None
        )
        fail_csv_file = open(_failed_file, "w") if _failed_file is not None else None

        # write header to the failed reads CSV
        if fail_csv_file is not None:
            fail_csv_file.write("read_id\tobserved_seq\tquality\treason_failed\tlib_call\n")

        # look through all sequences in the selection
        for i, seq_record in enumerate(
            tqdm(
                self.selection.get_sequence_reader(),
                desc=f"Running Decoding for selection {self.selection.selection_id}",
                disable=not use_tqdm,
            )
        ):
            self.decoder.decode_statistics.num_seqs_read += 1

            # decode the read
            decoded_barcode = self.decoder.decode_read(seq_record)
            # skip failed reads
            if isinstance(decoded_barcode, DecodedDELCompound):
                self.decoder.decode_statistics.num_seqs_decoded_per_lib[
                    decoded_barcode.library.library_id
                ] += 1
                if self.degen.count_barcode(decoded_barcode):
                    # only up the degen count if not a degenerate read
                    self.decoder.decode_statistics.num_seqs_degen_per_lib[
                        decoded_barcode.library.library_id
                    ] += 1
            elif fail_csv_file is not None:
                # if we are saving failed reads, save the read
                fail_csv_file.write(
                    f"{seq_record.id}\t"
                    f"{seq_record.sequence}\t"
                    f"{seq_record.qualities}\t"
                    f"{decoded_barcode.__class__.__name__}\n"
                )

            if ((i + 1) % 1000) == 0:
                self.logger.debug(
                    f"Decoded {i + 1} reads for selection {self.selection.selection_id}"
                )

        if fail_csv_file is not None:
            fail_csv_file.close()
            self.logger.debug(
                f"Saved failed reads to {save_failed_to} "
                f"for selection {self.selection.selection_id}"
            )

        self.logger.info(f"Completed decoding for selection: {self.selection.selection_id}")
        return DecodingRunnerResults(
            selection=self.selection,
            decode_statistics=self.decoder.decode_statistics,
            degen=self.degen,
        )

    def to_file(self, out_path: str | PathLike):
        """
        Write decode runner config to a human-readable file

        Parameters
        ----------
        out_path: str or PathLike
            path to save experiment to
        """
        data = {
            "selection_id": self.selection.selection_id,
            "target_id": self.selection.selection_condition.target_id,
            "selection_condition": self.selection.selection_condition.selection_condition,
            "additional_info": self.selection.selection_condition.additional_info
            if self.selection.selection_condition.additional_info
            else "NA",
            "data_ran": self.selection.get_run_date_as_str(),
            "sequence_files": self.selection.sequence_files,
            "libraries": [
                str(lib.loaded_from) for lib in self.selection.library_collection.libraries
            ],
            "decode_settings": self.decode_settings.__dict__,
        }
        yaml.safe_dump(data, open(out_path, "w"))

    @classmethod
    def from_file(
        cls,
        decode_file: str | PathLike,
        fastq_files: list[str | PathLike] | None = None,
        ignore_decode_seqs: bool = False,
        debug: bool = False,
        disable_logging: bool = False,
    ) -> "DecodingRunner":
        """
        Load the decode run from a human-readable file

        Notes
        -----
        This will raise an exception if both `fastq_files` and the "sequence_files" key
        in the decode file are not provided OR if both are provided.
        This is to avoid confusion about which files to use for decoding.
        Only one of these two parameters should be used to provide fastq files for the
        runner.

        The exception to this is when 'ignore_decode_seqs' is set to `True`.
        In this case, the observed_seq files used will always be the one provided
        to the function and the decode file sequences will be ignored.

        NOTE: it is best practice to add the sequences to the decode file.
        This way the full settings needed for the run are in a single file.
        The ability to also pass fastq files is for convenience and mainly
        meant to support parallelization efforts when splitting fastq files
        into smaller chunks on the fly.

        Parameters
        ----------
        decode_file: str or PathLike
            path to load experiment from
        fastq_files: list[str | PathLike], default = None
            list of paths to fastq files to decode
            if `None`, will use the observed_seq files from the decode file
        ignore_decode_seqs: bool, default = False
            if true, will ignore the observed_seq files in the decode file
            in this case, the `fastq_files` parameter must be provided
        debug: bool, default = False
            if true, will enable debug logging
        disable_logging: bool, default = False
            if true, will disable logging
            this is useful for running the experiment without logging
            or when running in a script where logging is not needed

        Raises
        ------
        DecodingRunParsingError
            raised if both `fastq_files` and the "sequence_files" key in the
            passed `decode_file` are not provided OR if both are provided and
            `ignore_decode_seqs` is `False`.

        Returns
        -------
        DecodingExperiment
        """
        # just to make it list for typing constancy
        _fastq_files: list[str | PathLike]
        if fastq_files is None:
            _fastq_files = []
        else:
            _fastq_files = fastq_files

        data = yaml.safe_load(open(decode_file, "r"))

        # check this early so we do not waste time loading libraries if format error
        _seq_files: list[str | PathLike]
        _decode_seq_files = data.get("sequence_files", [])
        if (len(_fastq_files) > 0) and (len(_decode_seq_files) > 0):
            # if we are ignoring the decode seqs, just use the fastq files
            if ignore_decode_seqs:
                _seq_files = _fastq_files
            else:
                raise DecodingRunParsingError(
                    f"`fastq_files` cannot be provided when "
                    f"{decode_file} has a 'sequence_files' section; "
                    f"did you mean to set `ignore_decode_seqs` to True?"
                )
        elif (len(_fastq_files) == 0) and (len(_decode_seq_files) == 0):
            raise DecodingRunParsingError(
                f"either {decode_file} should have a 'sequence_files' section OR "
                f"`fastq_files` is a non-None list of files; found neither"
            )
        elif (len(_fastq_files) == 0) and ignore_decode_seqs:
            raise DecodingRunParsingError(
                "``ignore_decode_seqs`` is True, but no `fastq_files` "
                "provided; pass `fastq_files` OR set `ignore_decode_seqs` to True"
            )
        elif len(_fastq_files) > 0:
            _seq_files = _fastq_files
        else:
            _seq_files = _decode_seq_files

        # load in the libraries
        try:
            _libraries: list[str] = data["libraries"]
        except KeyError as e:
            raise DecodingRunParsingError(
                f"{decode_file} decoding run config file does not contain a 'libraries' section"
            ) from e
        _library_collection = DELibraryCollection(
            [DELibrary.load(lib_path) for lib_path in _libraries]
        )

        _selection_id = data.get("selection_id", None)
        _target_id = data.get("target_id", None)
        _selection_condition = data.get("selection_condition", None)
        _additional_info = data.get("additional_info", None)
        _date_ran = data.get("data_ran", None)

        # parse date ran if given
        _date_ran_timestamp: datetime.datetime | None
        try:
            _date_ran_timestamp = datetime.datetime.fromisoformat(_date_ran) if _date_ran else None
        except ValueError:
            _date_ran_timestamp = None

        # make selection object
        _selection = SequencedSelection(
            library_collection=_library_collection,
            sequence_files=_seq_files,
            date_ran=_date_ran_timestamp,
            target_id=_target_id,
            selection_condition=_selection_condition,
            selection_id=_selection_id,
            additional_info=_additional_info,
        )

        # if decode setting provided, load otherwise use defaults
        _decode_settings: dict[str, Any] | None = data.get("decode_settings", None)
        if _decode_settings is None:
            _decode_setting_obj = DecodingSettings()
        else:
            try:
                _decode_setting_obj = DecodingSettings(**_decode_settings)
            except TypeError as e:
                _unknown_arg = e.args[0].split()[-1]
                raise DecodingRunParsingError(
                    f"unrecognized decoding settings: {_unknown_arg}"
                ) from e

        return cls(
            selection=_selection,
            decode_settings=_decode_setting_obj,
            disable_logging=disable_logging,
            debug=debug,
        )
