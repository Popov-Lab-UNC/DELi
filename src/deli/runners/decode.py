import datetime
import logging
import os
import random
import warnings
from pathlib import Path
from typing import Literal, Optional, Any, Iterator

import yaml
from tqdm import tqdm

from deli._logging import get_dummy_logger, get_logger

from deli.decode.base import FailedDecodeAttempt
from deli.decode.decoder import DecodeStatistics, DELCollectionDecoder, DecodedDELCompound
from deli.decode.degen import DELCollectionCounter
from deli.decode.library_demultiplex import LibraryDemultiplexer, FullSeqAlignmentLibraryDemultiplexer, \
    LibraryTagCutadaptLibraryDemultiplexer, SinglePrimerCutadaptLibraryDemultiplexer, \
    FlankingPrimersCutadaptLibraryDemultiplexer, LibraryTagRegexLibraryDemultiplexer, \
    SinglePrimerRegexLibraryDemultiplexer, FlankingPrimersRegexLibraryDemultiplexer, get_library_demultiplexer_type
from deli.decode.report import build_decoding_report

from deli.dels.combinatorial import DELibraryCollection, DELibrary
from deli.dna.io import get_reader
from deli.enumeration.enumerator import EnumerationRunError
from deli.selection import Selection, SequencedSelection


class DecodingRunParsingError(Exception):
    """Exception to raise when a decoding experiment file is invalid"""

    pass


# TODO this should be a dataclass
class DecodingSettings(dict):
    """
    Define parameters for decoding experiments

    More details about the exact effect of these settings can
    be found in the "Decoding" docs

    Notes
    -----
    Only parameters relating to the algorithm should be here
    Setting relating to IO should be handled outside this context

    Parameters
    ----------
    ignore_tool_compounds: bool, default = False
        if true, will ignore any tool compounds during decoding
    demultiplexer_algorithm: Literal["cutadapt", "regex", "full"], default = "regex"
        The demultiplexing algorithm to use.
        - "cutadapt": use a cutadapt to locate sections
        - "regex": use a regular expression based demultiplexer
        - "full": use a full alignment based demultiplexer
    demultiplexer_mode: Literal["library", "single", "flanking"], default = "flanking"
        The demultiplexing section strategy to use.
        - "library": demultiplex by matching just the library tag
        - "single": demultiplex by matching a single static barcode section
        - "flanking": demultiplex by matching barcode sections that flank the library tag
        (flanking means one before and one after the tag)
    realign: bool, default = False
        if true, will perform a local realignment of the read to the
        libraries barcode schema *after* demultiplexing determine the library.
        This could help recover reads that have complex alignments due multiple indels
    library_error_tolerance: int, default = 1
        The number of errors you are willing to tolerate in any given barcode
        section during library demultiplexing. Will apply to each section
        independently. For example, a flanking demultiplexer will allow for
        1 error in *each* of the flanking sections.
    library_error_correction_mode_str: str, default = "levenshtein_dist:2,asymmetrical"
        The error correction mode string to use for library barcode
        calling during demultiplexing.
    min_library_overlap: int , default = 8
        if using a cutadapt style demultiplexer, this is the minimum number of bases
        that must align to the expected barcode section for a match to be called.
        See the cutadapt documentation for more details on this parameter.
    wiggle: bool, default = False
        if true, will extend aligned sections by 1 bp on each side
        and then and scan all possible chunks of the expected barcode length,
        1 smaller and 1 larger than expected length (in that order). If not
        using a local realignment post demultiplexing this can help recover
        reads lost to indels in the barcode region.
    revcomp: bool, default = False
        If true, search the reverse compliment as well.
        In most cases it is faster to use an external tools
        to align and reverse compliment reads before decoding
    max_read_length: int or None, default = None
        maximum length of a read to be considered for decoding
        if above the max, decoding will fail
        if `None` will default to 5x the min_read_length
    min_read_length: int or None, default = None
        minimum length of a read to be considered for decoding
        if below the min, decoding will fail
        if `None` will default to the smallest min match length of
        any library in the collection considered for decoding
        with 10bp of buffer
    default_error_correction_mode_str: str, default = "levenshtein_dist:1,asymmetrical"
        The default error correction mode string to use for decoding.
        If a barcode section lacks a specified error correction mode,
        this mode will be used.
        See the documentation for more details on the format of this string.
    """

    def __init__(
        self,
        ignore_tool_compounds: bool = False,
        demultiplexer_algorithm: Literal["cutadapt", "regex", "full"] = "regex",
        demultiplexer_mode: Literal["library", "single", "flanking"] = "single",
        realign: bool = False,
        library_error_tolerance: int = 1,
        library_error_correction_mode_str: str = "levenshtein_dist:2,asymmetrical",
        min_library_overlap: int = 8,
        revcomp: bool = False,
        wiggle: bool = False,
        max_read_length: Optional[int] = None,
        min_read_length: Optional[int] = None,
        default_error_correction_mode_str: str = "levenshtein_dist:1,asymmetrical",
    ):
        super().__init__(
            ignore_tool_compounds=ignore_tool_compounds,
            demultiplexer_algorithm=demultiplexer_algorithm,
            demultiplexer_mode=demultiplexer_mode,
            realign=realign,
            library_error_tolerance=library_error_tolerance,
            library_error_correction_mode_str=library_error_correction_mode_str,
            min_library_overlap=min_library_overlap,
            revcomp=revcomp,
            wiggle=wiggle,
            max_read_length=max_read_length,
            min_read_length=min_read_length,
            default_error_correction_mode_str=default_error_correction_mode_str,
        )

    def to_file(self, path: str):
        """
        Save settings to a YAML file

        Parameters
        ----------
        path : str
            path to save settings to
        """
        yaml.dump(dict(self), open(path, "w"))

    @classmethod
    def from_file(cls, path: str) -> "DecodingSettings":
        """
        Load settings from a yaml file

        Will first check if there is a "decode_settings" key
        and load settings from that sub dict.
        Otherwise will load from the yaml file keys

        Parameters
        ----------
        path : str
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
    """
    Holds the results of a decoding run

    Parameters
    ----------
    selection: Selection
        the selection that was decoded
    decode_statistics: DecodeStatistics | None, default = None
        the decode statistics for the run
    degen: DELCollectionCounter | None, default = None
        the degeneration counter for the run
    """

    def __init__(
        self,
        selection: Selection,
        decode_statistics: DecodeStatistics | None = None,
        degen: DELCollectionCounter | None = None,
    ):
        self.selection = selection
        self.decode_statistics = decode_statistics
        self.degen = degen

        self._max_cycle_size: int = max([library.num_cycles for library in self.selection.library_collection])

    def write_decode_report(self, out_path: str):
        """
        Write the decoding statistics for this run

        Parameters
        ----------
        out_path: str
            path to save decoding statistics to
        """
        if self.decode_statistics is not None:
            if os.path.dirname(out_path) != "":
                os.makedirs(os.path.dirname(out_path), exist_ok=True)
            build_decoding_report(selection=self.selection, stats=self.decode_statistics, out_path=out_path)
        else:
            raise RuntimeError("Cannot write decode report, no decode statistics found")

    def write_decode_statistics(self, out_path: str, include_read_lengths: bool = False):
        """
        Write the decoding statistics for this run

        Parameters
        ----------
        out_path: str
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
        out_path: str,
        include_library_id_col: bool = True,
        include_bb_id_cols: bool = True,
        include_bb_smi_cols: bool = False,
        include_raw_count_col: bool = True,
        enumerate_smiles: bool = False,
        fail_on_failed_numeration: bool = False,
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
        out_path: str
            path to save results to
            will create the directory if it does not exist
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
        fail_on_failed_numeration: bool, default = False
            if `True`, will raise an exception if enumeration fails for any DEL
            if `False`, will set the SMILES to 'null' if enumeration fails
            only used if `enumerate_smiles` is `True`
        file_format: Literal["csv", "tsv"] = "csv"
            which file format to write to
            either "csv" or "tsv"

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
                    lib.library_id for lib in self.selection.library_collection.libraries if not lib.can_enumerate()
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
            self.selection.library_collection.max_cycle_size() if (include_bb_id_cols or include_bb_smi_cols) else 0
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
                    for decoded_compound, counts in library_count_set.items():
                        _row_dict = decoded_compound.to_cube_row_dict()
                        _row = _row_dict["compound_id"]
                        if include_library_id_col:
                            _row += f"{delimiter}{_row_dict.get('library_id', 'null')}"
                        for i in range(1, self._max_cycle_size + 1):
                            if include_bb_id_cols:
                                _row += f"{delimiter}{_row_dict.get(f'BB{i:02d}_id', 'null')}"
                            if include_bb_smi_cols:
                                _row += f"{delimiter}{_row_dict.get(f'BB{i:02d}_smiles', 'null')}"

                        if enumerate_smiles:
                            try:
                                _smi = decoded_compound.get_smiles()
                            except (EnumerationRunError, RuntimeError, ValueError) as e:
                                if fail_on_failed_numeration:
                                    raise RuntimeError(
                                        f"Failed to enumerate SMILES for decoded DEL ID {decoded_compound.compound_id}"
                                    ) from e
                                _smi = "null"
                            _row += f"{delimiter}{_smi}"
                        if include_raw_count_col:
                            _row += f"{delimiter}{counts.get_raw_count()}"
                        _row += f"{delimiter}{counts.get_degen_count()}\n"
                        f.write(_row)
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

        def __init__(
                self,
                selection: SequencedSelection,
                decode_settings: DecodingSettings | None = None,
                disable_logging: bool = False,
                debug: bool = False,
                run_id: str | None = None,
        ):
            if run_id is None:
                self.run_id = f"DELi-DECODER-{os.getpid()}-" + f"{random.randint(0, 100000):06d}"
            else:
                self.run_id = run_id

            self.logger: logging.Logger = get_dummy_logger() if disable_logging else get_logger(self.run_id,
                                                                                                debug=debug)

            self.selection: SequencedSelection = selection
            self.decode_settings: DecodingSettings = decode_settings if decode_settings is not None else DecodingSettings()
            self.decode_stats: DecodeStatistics = DecodeStatistics()

            # parse the demultiplexer settings
            demultiplexer: LibraryDemultiplexer

            demultiplex_algorithm = self.decode_settings["demultiplexer_algorithm"]
            demultiplex_mode = self.decode_settings["demultiplexer_mode"]


            if demultiplex_algorithm == "full":
                demultiplexer = FullSeqAlignmentLibraryDemultiplexer(
                    libraries=self.selection.library_collection,
                    tool_compounds=list(self.selection.tool_compounds),
                    revcomp=self.decode_settings["revcomp"],
                )
            elif demultiplex_algorithm == "cutadapt":
                if demultiplex_mode == "library":
                    demultiplexer = LibraryTagCutadaptLibraryDemultiplexer(
                        libraries=self.selection.library_collection,
                        tool_compounds=list(self.selection.tool_compounds),
                        min_overlap=self.decode_settings["min_library_overlap"],
                        error_tolerance=self.decode_settings["library_error_tolerance"],
                        revcomp=self.decode_settings["revcomp"],
                        realign=self.decode_settings["realign"],
                    )
                elif demultiplex_mode == "single":
                    demultiplexer = SinglePrimerCutadaptLibraryDemultiplexer(
                        libraries=self.selection.library_collection,
                        tool_compounds=list(self.selection.tool_compounds),
                        min_overlap=self.decode_settings["min_library_overlap"],
                        error_tolerance=self.decode_settings["library_error_tolerance"],
                        revcomp=self.decode_settings["revcomp"],
                        error_correction_mode_str=self.decode_settings["default_error_correction_mode_str"],
                        realign=self.decode_settings["realign"],
                    )
                elif demultiplex_mode == "flanking":
                    demultiplexer = FlankingPrimersCutadaptLibraryDemultiplexer(
                        libraries=self.selection.library_collection,
                        tool_compounds=list(self.selection.tool_compounds),
                        min_overlap=self.decode_settings["min_library_overlap"],
                        error_tolerance=self.decode_settings["library_error_tolerance"],
                        revcomp=self.decode_settings["revcomp"],
                        error_correction_mode_str=self.decode_settings["default_error_correction_mode_str"],
                        realign=self.decode_settings["realign"],
                    )
                else:
                    raise DecodingRunParsingError(
                        f"demultiplexer_sections '{demultiplex_mode}' not recognized for approach 'cutadapt'"
                    )
            elif demultiplex_algorithm == "regex":
                if demultiplex_mode == "library":
                    demultiplexer = LibraryTagRegexLibraryDemultiplexer(
                        libraries=self.selection.library_collection,
                        tool_compounds=list(self.selection.tool_compounds),
                        error_tolerance=self.decode_settings["library_error_tolerance"],
                        revcomp=self.decode_settings["revcomp"],
                        realign=self.decode_settings["realign"],
                    )
                elif demultiplex_mode == "single":
                    demultiplexer = SinglePrimerRegexLibraryDemultiplexer(
                        libraries=self.selection.library_collection,
                        tool_compounds=list(self.selection.tool_compounds),
                        error_tolerance=self.decode_settings["library_error_tolerance"],
                        revcomp=self.decode_settings["revcomp"],
                        error_correction_mode_str=self.decode_settings["default_error_correction_mode_str"],
                        realign=self.decode_settings["realign"],
                    )
                elif demultiplex_mode == "flanking":
                    demultiplexer = FlankingPrimersRegexLibraryDemultiplexer(
                        libraries=self.selection.library_collection,
                        tool_compounds=list(self.selection.tool_compounds),
                        error_tolerance=self.decode_settings["library_error_tolerance"],
                        revcomp=self.decode_settings["revcomp"],
                        error_correction_mode_str=self.decode_settings["default_error_correction_mode_str"],
                        realign=self.decode_settings["realign"],
                    )
                else:
                    raise DecodingRunParsingError(
                        f"demultiplexer_sections '{demultiplex_mode}' not recognized for approach 'regex'"
                    )
            else:
                raise DecodingRunParsingError(f"demultiplexer_mode '{demultiplex_algorithm}' not recognized")

            # initialize all the decoding object required
            self.decoder = DELCollectionDecoder(
                library_demultiplexer=demultiplexer,
                decode_statistics=self.decode_stats,
                wiggle=self.decode_settings["wiggle"],
                max_read_length=self.decode_settings["max_read_length"],
                min_read_length=self.decode_settings["min_read_length"],
                default_error_correction_mode_str=self.decode_settings["default_error_correction_mode_str"],
            )

        def run(self, use_tqdm: bool = False) -> Iterator[DecodedDELCompound | FailedDecodeAttempt]:


    def run_to_file(
        self,
        save_to: os.PathLike = "./decodes.tsv",
        split_by_library: bool = False,
        include_fastq_info: bool = False,
        save_failed_to: Optional[os.PathLike] = None,
        use_tqdm: bool = False,
    ) -> list[Path]:
        """
        Run the decoder

        Notes
        -----
        The output decode files are written to on the fly to both save memory
        and to allow for stopping and restarting.
        The file format is always a TSV with the following columns:
        - FASTQ-source (if `include_fastq_info` is True)
        - READ_ID (if `include_fastq_info` is True)
        - DEL_ID
        - LIB_ID
        - BB<#>_ID (for each cycle #)
        - UMI (if present, otherwise 'null')

        The purpose of the TSV format is to avoid issues with commas in various ids

        Parameters
        ----------
        save_to: PathLike, default = "decodes.tsv"
            save all decodes to this location.
            Can be a file (ends in file extension) or a directory.
            If directory, will save a "decodes.tsv" file in that directory.
            The file will always be in TSV format, even if a different
            file extension is provided.
            If `split_by_library` is True, will save to
            <save_to>/<library_id>_decodes.tsv (ff `save_to` is a file,
            will drop the extension and make
            a directory with that name to save the files to).
        split_by_library: bool, default = False
            if true, will save decodes to separate files
            for each library in the collection.
            files will be named <library_id>_decodes.tsv
        include_fastq_info: bool, default = False
            if true, will include the fastq source file
            and read id in the decode output file(s)
        save_failed_to: PathLike, default = None
            if provided, will save failed reads to this directory
            file will be named <selection_id>_decode_failed.tsv
            will include the read_id, the barcode, the quality chain,
            and reason failed
        use_tqdm: bool, default = False
            turn on a tqdm tracking bar
            only recommended if running a single
            runner

        Returns
        -------
        list[Path]
            list of paths to the decode output file(s)
            can be anywhere from 1 to N files, where
            N is the number of libraries in the selection
        """
        self.logger.info(f"Running decoding for selection: {self.selection.selection_id}")

        # parse save_to path
        save_to_path = Path(save_to)

        if save_to_path.is_dir() or save_to_path.suffix:
            save_to_dir_path = save_to_path.absolute()
            save_to_file_path = save_to_dir_path / "decodes.tsv"
        else:
            save_to_dir_path = save_to_path.parent.absolute()
            save_to_file_path = save_to_path.with_suffix(".tsv").absolute()

        os.makedirs(save_to_dir_path, exist_ok=True)

        # deal with logging save path
        if split_by_library:
            self.logger.info(
                f"Saving library specific decode files to directory {str(save_to_dir_path)}"
            )
        else:
            self.logger.info(
                f"Saving all decodes to {str(save_to_file_path)}"
            )

        # to avoid opening and closing files repeatedly if there are multiple
        output_files = {"all": open(save_to_file_path, "w")} if not split_by_library else {}

        # open failed reads file if specified
        failed_file_path: Path | None = None
        if save_failed_to is not None:
            failed_file_path = Path(save_failed_to) / f"{self.selection.selection_id}_decode_failed.tsv"
            self.logger.info(f"Saving failed reads to {str(failed_file_path)}")

        fail_tsv_file = open(failed_file_path, "w") if failed_file_path is not None else None
        # write header to the failed reads TSV
        if fail_tsv_file is not None:
            fail_tsv_file.write("read_id\tbarcode\tquality\tfail_type\tfail_desc\n")

        # look through all sequences in the selection
        for i, (fastq_file, seq_record) in enumerate(
            tqdm(
                self.selection.sequence_reader.iter_seqs_with_filenames(),
                desc=f"Running Decoding for selection {self.selection.selection_id}",
                disable=not use_tqdm,
            )
        ):
            self.decode_stats.num_seqs_read += 1
            # decode the read
            decoded_barcode = self.decoder.decode_read(seq_record)

            # process decoded reads only
            if isinstance(decoded_barcode, DecodedDELCompound):
                self.decode_stats.num_seqs_decoded_per_lib[decoded_barcode.library.library_id] += 1
                # determine which output file to use
                if split_by_library:
                    if decoded_barcode.library.library_id not in output_files:
                        lib_file_path = save_to_dir_path / f"{decoded_barcode.library.library_id}_decodes.tsv"
                        output_files[decoded_barcode.library.library_id] = open(lib_file_path, "w")
                    curr_decode_out_file = output_files[decoded_barcode.library.library_id]
                else:
                    curr_decode_out_file = output_files["all"]
                # write fastq info if specified
                output_str = f"{seq_record.id}\t{decoded_barcode.compound_id}\t" if include_fastq_info else ""
                # write decoded barcode info
                output_str += (
                    f"{decoded_barcode.compound_id}\t"
                    f"{decoded_barcode.library.library_id}\t"
                    
                    f"{decoded_barcode.umi if decoded_barcode.has_umi() else 'null'}\n"
                )
                curr_decode_out_file.write(output_str)

            # write failed reads if specified
            elif (isinstance(decoded_barcode, FailedDecodeAttempt)) and (fail_tsv_file is not None):
                # if we are saving failed reads, save the read
                fail_tsv_file.write(
                    f"{seq_record.id}\t"
                    f"{seq_record.sequence}\t"
                    f"{seq_record.qualities}\t"
                    f"{decoded_barcode.__class__.__name__}\t"
                    f"{decoded_barcode.reason}\n"
                )

            if ((i + 1) % 100000) == 0:
                self.logger.debug(f"Decoded {i + 1} reads for selection {self.selection.selection_id}")

        # close open files
        for f in output_files.values():
            f.close()
        if fail_tsv_file is not None:
            fail_tsv_file.close()

        self.logger.info(f"Completed decoding for selection: {self.selection.selection_id}")

        return list([Path(f.name).absolute() for f in output_files.values()])

    def to_file(self, out_path: str):
        """
        Write decode runner config to a human-readable file

        Parameters
        ----------
        out_path: str
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
            "sequence_files": [str(p) for p in self.selection.sequence_files],
            "libraries": [str(lib.loaded_from) for lib in self.selection.library_collection.libraries],
            "decode_settings": self.decode_settings.__dict__,
        }
        yaml.safe_dump(data, open(out_path, "w"))

    @classmethod
    def from_file(
        cls,
        decode_file: str,
        fastq_files: list[os.PathLike] | None = None,
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
        In this case, the barcode files used will always be the one provided
        to the function and the decode file sequences will be ignored.

        NOTE: it is best practice to add the sequences to the decode file.
        This way the full settings needed for the run are in a single file.
        The ability to also pass fastq files is for convenience and mainly
        meant to support parallelization efforts when splitting fastq files
        into smaller chunks on the fly.

        Parameters
        ----------
        decode_file: str
            path to load experiment from
        fastq_files: list[str], default = None
            list of paths to fastq files to decode
            if `None`, will use the barcode files from the decode file
        ignore_decode_seqs: bool, default = False
            if true, will ignore the barcode files in the decode file
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
            raised if both `fastq_files` and the "sequence_files" key in the `decode_file`
            are not provided OR if both are provided and `ignore_decode_seqs` is `False`.

        Returns
        -------
        DecodingExperiment
        """
        # just to make it list for typing constancy
        _fastq_files: list[os.PathLike]
        if fastq_files is None:
            _fastq_files = []
        else:
            _fastq_files = fastq_files

        data = yaml.safe_load(open(decode_file, "r"))

        # check this early so we do not waste time loading libraries if format error
        _seq_files: list[os.PathLike] = []
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

        _seq_reader = get_reader(_seq_files)

        # load in the libraries
        try:
            _libraries: list[str] = data["libraries"]
        except KeyError as e:
            raise DecodingRunParsingError(
                f"{decode_file} decoding run config file does not contain a 'libraries' section"
            ) from e
        _library_collection = DELibraryCollection([DELibrary.load(lib_path) for lib_path in _libraries])

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
            sequence_reader=_seq_reader,
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
                raise DecodingRunParsingError(f"unrecognized decoding settings: {_unknown_arg}") from e

        return cls(
            selection=_selection,
            decode_settings=_decode_setting_obj,
            disable_logging=disable_logging,
            debug=debug,
        )
