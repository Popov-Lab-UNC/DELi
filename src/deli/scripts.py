import os
import json
import argparse
import logging
import pickle
import subprocess
import sys
import time
import warnings
from typing import Literal

from tqdm import tqdm

from deli.constants import MODULE_DIR, FILE_HEADER, INDEX, LIBS, MAKES, DeliConfigError
from deli.utils import lib_to_make
from deli.html_report import HTMLReport

logging.captureWarnings(True)  # hook warning to logger


def _wrapper_handle_dup(func, *args, **kwargs):
    """
    A  wrapper script to handle argparse receiving the same argument twice without crashing

    Parameters
    ----------
    func: Callable
        some function (parser.add_argument in this case)
    args:
        positional args to pass to function
    kwargs
        keyword args to pass to function
    """
    try:
        func(*args, **kwargs)
    except argparse.ArgumentError:
        pass


def _add_general_args(parser: argparse.ArgumentParser):
    """
    Adds general DELi arguments to argparser

    Parameters
    ----------
    parser: ArgumentParser
        argument parser to add args to

    Notes
    -----
    edits argparser inplace
    """
    _wrapper_handle_dup(parser.add_argument, '-o', '--outdir', type=str, default=os.getcwd(),
                        help='directory to write the output files to')
    
    _wrapper_handle_dup(parser.add_argument, '-p', '--prefix', type=str, default="",
                        help='prefix to add to the filename of output called csv file(s)')

    _wrapper_handle_dup(parser.add_argument, '--debug', action='store_true', default=False,
                        help='enable debug logging')

    _wrapper_handle_dup(parser.add_argument, '--print', action='store_true',
                        help="print updating status to console")
    
    
def _add_matching_args(parser: argparse.ArgumentParser):
    """
    Adds the DELi matching arguments to argparser

    Parameters
    ----------
    parser: ArgumentParser
        argument parser to add args to

    Notes
    -----
    edits argparser inplace
    """
    _wrapper_handle_dup(parser.add_argument, '--make', required=True, type=str,
                        help="the 'make' ID defining the makeup of the library")

    _wrapper_handle_dup(parser.add_argument, '-s', '--save-match', action='store_true',
                        help='whether to save the matches as a file')

    _wrapper_handle_dup(parser.add_argument, '-r', "--recursive", action="store_true",
                        help='recursively loop through all file in a directory')

    _wrapper_handle_dup(parser.add_argument, '-m', "--mode", choices=['full', 'minimum', 'double', 'single'],
                        default='single',
                        help="which barcode query to use when searching for matches, see readme.md for more details")

    _wrapper_handle_dup(parser.add_argument, '--strict', action='store_true',
                        help='use a more strict alignment and pattern match that uses the closing primer post UMI')

    _wrapper_handle_dup(parser.add_argument, '-e', '--error', type=int, required=False, default=0,
                        help='how many edit distance errors you are willing to allow when attempting a match')

    _wrapper_handle_dup(parser.add_argument, '-c', '--compliment', action='store_true',
                        help='also search for the reverse compliment of the DEL')

    _wrapper_handle_dup(parser.add_argument, '-j', '--n-workers', type=int, default=-1,
                        help='number of worked to use for calling')


def _add_calling_args(parser: argparse.ArgumentParser):
    """
    Adds the DELi calling arguments to argparser

    Parameters
    ----------
    parser: ArgumentParser
        argument parser to add args to

    Notes
    -----
    edits argparser inplace
    """
    _wrapper_handle_dup(parser.add_argument, '--make', required=True, type=str,
                        help="the 'make' ID defining the makeup of the library")

    _wrapper_handle_dup(parser.add_argument, '-s', '--save-match', action='store_true',
                        help='whether to save the matches as a file')

    _wrapper_handle_dup(parser.add_argument, '--strict', action='store_true',
                        help='use a more strict alignment and pattern match that uses the closing primer post UMI')
    
    _wrapper_handle_dup(parser.add_argument, '-i', "--index", type=str, nargs="+", required=True,
                        help='index(s) included in the DEL selection being called; pass "all" for all indexes')

    _wrapper_handle_dup(parser.add_argument, '-l', "--library", type=str, nargs="+", required=True,
                        help='library(s) included in the DEL selection being called; pass "all" for all library(s)')
    
    _wrapper_handle_dup(parser.add_argument, '-j', '--n-workers', type=int, default=-1,
                        help='number of worked to use for calling')
    
    
def _add_cube_gen_args(parser: argparse.ArgumentParser):
    """
    Adds the DELi cube_gen arguments to argparser

    Parameters
    ----------
    parser: ArgumentParser
        argument parser to add args to

    Notes
    -----
    edits argparser inplace
    """
    _wrapper_handle_dup(parser.add_argument, '-i', "--index", type=str, nargs="+", required=True,
                        help='index(s) included in the DEL selection being called; pass "all" for all indexes')

    _wrapper_handle_dup(parser.add_argument, "--index_name", type=str, required=False, default=None,
                        help='path to json file containing a dictionary mapping passed index to selection names')
    
    _wrapper_handle_dup(parser.add_argument, '-l', "--library", type=str, nargs="+", required=True,
                        help='library(s) included in the DEL selection being called; pass "all" for all library(s)')

    _wrapper_handle_dup(parser.add_argument, '-t', "--control", type=str, required=False, default=None,
                        help='index of selection that is the negative control')
    
    _wrapper_handle_dup(parser.add_argument, '-n', '--normalize', action='store_true',
                        help='normalize the DEL_ID sequence counts by the control (requires -t/--control is set)')

    _wrapper_handle_dup(parser.add_argument, '--monosynthon', action='store_true',
                        help='conduct a monosynthon analysis on all generated cubes')

    _wrapper_handle_dup(parser.add_argument, '--disynthon', action='store_true',
                        help='conduct a monosynthon analysis on all generated cubes')
    
    _wrapper_handle_dup(parser.add_argument, '-u', '--umi-cluster', action='store_true',
                        help='conduct a UMI clustering to further denoise results')
    
    _wrapper_handle_dup(parser.add_argument, '-j', '--n-workers', type=int, default=-1,
                        help='number of worked to use for calling')


def get_args(parser: argparse.ArgumentParser, mode: Literal["matching", "calling", "cube_gen", "all"]):
    """
    Adds the correct argument to the passed argparser based on DELi step

    Parameters
    ----------
    parser: ArgumentParser
        parser to add args to
    mode: Literal["matching", "calling", "cube_gen", "all"]
        which steps need to be added (all will add all steps)

    Notes
    -----
    edits argparser inplace
    """
    _add_general_args(parser)
    if mode in ["matching", "all"]:
        _add_matching_args(parser)
    if mode in ["calling", "all"]:
        _add_calling_args(parser)
    if mode in ["cube_gen", "all"]:
        _add_cube_gen_args(parser)


def _check_prefix(args: argparse.Namespace):
    """
    Cleans up the prefix argument to have nice underscoring when used in a script.

    Notes
    -----
    Edits the argument inplace (will NOT return a new argument Namespace)

    Parameters
    ----------
    args: Namespace
        the parsed Namespace or arguments passes to the script
    """
    if not args.prefix.endswith("_"):
        args.prefix += "_"


def _check_index(args: argparse.Namespace):
    """
    Handles checking the index arguments to make sure:
     1. The index(s) passed exist
     2. The "all" keyword is converted to a list of all the index names

    Notes
    -----
    Edits the argument inplace (will NOT return a new argument Namespace)

    Parameters
    ----------
    args: Namespace
        the parsed Namespace or arguments passes to the script

    Raises
    ------
    DeliConfigError
        "passed index {idx} does not exist in DELi configs"
    """
    if args.index[0] == "all":
        args.index = list(INDEX.keys())
    else:
        for idx in args.index:
            if idx not in INDEX.keys():
                raise DeliConfigError(f"passed index {idx} does not exist in DELi configs")


def _check_libs(args: argparse.Namespace):
    """
    Handles checking the library arguments to make sure:
     1. The library(s) passed exist
     2. The "all" keyword is converted to a list of all the library names

    Notes
    -----
    Edits the argument inplace (will NOT return a new argument Namespace)

    Parameters
    ----------
    args: Namespace
        the parsed Namespace or arguments passes to the script

    Raises
    ------
    DeliConfigError
        "passed library {lib} does not exist in DELi configs"
    """
    if args.library[0] == "all":
        args.library = list(LIBS.keys())
    else:
        for lib in args.library:
            if lib not in LIBS.keys():
                raise DeliConfigError(f"passed library {lib} does not exist in DELi configs")


def _check_make(args: argparse.Namespace, check_lib: bool = True):
    """
    Handles checking the make arguments to make sure:
     1. The 'make' passed exist
     2. All passed library(s) are built using the passed barcode make

    Notes
    -----
    If checking library(s) will raise a warning if a passed library is does not have the correct 'make'
    and then remove it form the list of library(s) used during downstream execution

    Parameters
    ----------
    args: Namespace
        the parsed Namespace or arguments passes to the script
    check_lib: bool
        whether to also make sure that the passed library(s) are built using the passed barcode make

    Raises
    ------
    DeliConfigError
        Raised when passed make is not found in config files
        "passed library {lib} does not exist in DELi configs"

    Warnings
    --------
    Warning
        raised when `check_lib` is `True` and a library(s) is found using the wrong make
        "library {lib} does not have make {args.make}; removing this library"
    """
    if args.make not in MAKES.keys():
        raise DeliConfigError(f"passed 'make' {args.make} does not exist in DELi configs")
    if check_lib:
        _checked_libs = []
        for lib in args.library:
            if lib_to_make(lib) != args.make:
                warnings.warn(f"library {lib} does not have make {args.make}; removing this library")
            else:
                _checked_libs.append(lib)
        args.library = _checked_libs


def check_args(args: argparse.Namespace, check_prefix: bool = True, check_lib: bool = True,
               check_index: bool = True, check_make: bool = True):
    """
    Check the passed DEL arguments to make sure they are compatible with DELi configs

    Parameters
    ----------
    args: Namespace
        the parsed arguments to check
    check_prefix: bool
        check the `prefix` argument
    check_lib: bool
        check the `library` argument
    check_index: bool
        check the `index` argument
    check_make: bool
        check the `make` argument

    Notes
    -----
    Will alter the Namespace inplace if changes are needed/occur

    Raises
    ------
    DeliConfigError
        raised if there is any issue with the argument check finding argument values that are not found in the DELi
        configs (located in '<PACKAGE_DIR>/data')
    Warning
        raised if the check finds incompatible settings that it can fix.
        for example, a library with the wrong make (which can be removed and the execution can still progress)
    """
    if check_prefix:
        _check_prefix(args)
    if check_lib:
        _check_libs(args)
    if check_index:
        _check_index(args)
    if check_make:
        _check_make(args)


def run_matching(args: argparse.Namespace, logger: logging.Logger, report: HTMLReport) -> list[str]:
    """
    This script runs the "matching" job, the first of the three steps DELi needs to do for barcode calling

    Parameters
    ----------
    args: Namespace
        the arguments defined run settings
    logger: Logger
        a logger object to handle logging
    report: HTMLReport
        a report object to handle report generation

    Notes
    -----
    will save matches to match file (see readme for information on how match files are formated)

    Returns
    -------
    match_files: list[str]
        a list of all the absolute paths to the "match files" generated during the matching run
    """
    import fastq
    from deli.utils import reverse_compliment
    from deli.del_make import DELMake

    # look up make settings
    make = DELMake(args.make)

    # for legacy reasons, code needs fastq files to be a list
    if args.recursive:
        logger.info("searching for fastq files in {}".format(args.input))
        fastq_files = [os.path.join(args.input, f) for f in os.listdir(args.input) if f.endswith('.fastq')]
        logger.info("found {} fastq files".format(len(fastq_files)))
    else:
        fastq_files = [args.input]

    # determine the number of sub-jobs to launch
    _pool_size = args.n_workers if args.n_workers != -1 else os.cpu_count() - 1

    # determine the matching pattern to use (defined in the `barcode_makes.json` config file)
    _pattern = None
    if args.error <= 0:
        if args.mode == "full":
            _pattern = make["full"]
        if args.mode == "minimum":
            _pattern = make["minimal"]
        if args.mode == "double":
            _pattern = make["double_barcode"]
        if args.mode == "single":
            # only the single matching mode has a strict matching format
            if args.strict:
                _pattern = make["strict_single_barcode"]
            else:
                _pattern = make["single_barcode"]
    else:
        # if using error tolerance, add that into the regex string
        logger.debug("error mode {} active".format(args.error))
        if args.mode == "full":
            _pattern = make["full_error"].replace("ERR", str(args.error))
        if args.mode == "minimum":
            _pattern = make["minimal_error"].replace("ERR", str(args.error))
        if args.mode == "double":
            _pattern = make["double_barcode_error"].replace("ERR", str(args.error))
        if args.mode == "single":
            # only the single matching mode has a strict matching format
            if args.strict:
                _pattern = make["strict_single_barcode_error"].replace("ERR", str(args.error))
            else:
                _pattern = make["single_barcode_error"].replace("ERR", str(args.error))

    # raise an error if no pattern was found
    if _pattern is None:
        raise ValueError(f"invalid mode and error settings; mode:{args.mode}, error:{args.error}")
    logger.debug("searching for pattern {}".format(_pattern))

    # log some reporting data
    report.from_file = fastq_files[0] if len(fastq_files) == 1 else str(fastq_files)
    report.error_tolerance = args.error
    report.pattern = _pattern

    # read in the all the sequences from the fastq file into a list
    seqs = []
    lengths = []  # keep track of lengths for the reporter to use later
    print("reading sequence files...")
    for i, fastq_file in enumerate(fastq_files):  # loop through all the fastq files (deprecated; only ever 1 file now)
        _seqs = []
        _j = 0
        logger.info("reading sequences from {}".format(fastq_file))
        for j, read in enumerate(fastq.read(fastq_file)):
            _seqs.append((f"{i}_{_j}", read.getSeq()))
            _j += 1

        lengths.extend([len(_[1]) for _ in _seqs])  # for the reporter

        # if compliment mode is on also get the revcomp sequences are all read in sequences
        if args.compliment:
            logger.info("reading rev comp sequences from {}".format(fastq_file))
            _seqs = _seqs + [(f"{i}_{_j + ii}", reverse_compliment(s[1])) for ii, s in enumerate(_seqs)]

        seqs += _seqs
        logger.info(f"read in {len(seqs) if not args.compliment else len(seqs) // 2} sequences from {fastq_file}")

    _total_reads = len(seqs) if not args.compliment else len(seqs) // 2  # total reads is doubled if compliment is on
    report.total_reads = _total_reads
    report.seq_lengths = lengths

    _too_long = (len(make["BARCODE_FULL_SEQUENCE"]) * 2) + 120  # too long is a sequence over twice the expected length

    seqs = [seq for seq in seqs if len(seq[1]) < _too_long]  # remove the sequences that are too long

    # do some logging and reporting
    report.right_size = len(seqs) if not args.compliment else len(seqs) // 2
    logger.info(f"found {len(seqs) if not args.compliment else len(seqs) // 2} sequences <{_too_long} bp")

    # write N sequence files (1 per sub-job)
    _seq_files = []
    _batch_size = (len(seqs) // _pool_size) + 1  # batch size is the number of seqs per file
    for i in range(_pool_size):
        _seqs = seqs[i * _batch_size:(i + 1) * _batch_size]
        _seq_files.append(os.path.join(args.outdir, args.prefix + f"seqs_{i}.tsv"))
        with open(_seq_files[-1], "w") as f:
            for _s in _seqs:
                f.write(f"{_s[0]}\t{_s[1]}\n")

    # launch sub-jobs; 1 job per sequence file created
    _processes = []
    _match_files = []
    for i, _seq_file in enumerate(_seq_files):
        _match_file_loc = os.path.join(args.outdir,
                                       args.prefix + f"_match_{_seq_file.split('_')[-1].split('.')[0]}.pkl")
        _match_files.append(str(_match_file_loc))
        p = subprocess.Popen([sys.executable, str(os.path.join(MODULE_DIR, "matching_subprocess.py")), _pattern,
                              _seq_file, _match_files[-1]])
        logger.debug(f"launched matching subprocess {p.pid}")
        _processes.append(p)

    # wait for all sub-jobs to terminate
    _p_bar = tqdm(total=len(seqs), desc=f"matching sequences", disable=not args.print)
    _pbar_count = 0
    while _processes:
        time.sleep(0.5)
        if _processes[-1].poll() is None:  # check if the last job in the list is done
            _tmp_count = 0
            # loop through the tmp `.count` files the sub-jobs make and use them to determine the progress
            for file in [f for f in os.listdir(args.outdir) if f.endswith(".count")]:
                _file = os.path.join(args.outdir, file)
                try:
                    _tmp_count += int(open(_file, "r").readline().strip())
                except ValueError:
                    continue
            _p_bar.update(_tmp_count - _pbar_count)  # update the progress bar
            _pbar_count = _tmp_count
            continue
        else:
            _processes.pop()  # if a job is done remove it from the list of active jobs
    _p_bar.update(len(seqs)-_pbar_count)  # finish the tqdm bar for aesthetics
    _p_bar.close()

    # remove count files created by the sub-jobs
    for file in [f for f in os.listdir(args.outdir) if f.endswith(".count")]:
        os.remove(os.path.join(args.outdir, file))

    # remove seq files
    for seq_file in _seq_files:
        os.remove(seq_file)

    logger.debug("all matching subprocesses complete")
    return _match_files


def run_calling(args: argparse.Namespace, logger: logging.Logger, report: HTMLReport, input_: list[str] = None) -> str:
    """
    This script runs the "calling" job, the second of the three steps DELi needs to do for barcode calling

    Parameters
    ----------
    args: Namespace
        the arguments used to define the run settings
    logger: Logger
        the logger object for logging
    report: HTMLReport
        the report object to track reporting
    input_: list[str], default = None
        the list of match files to use
        if not passed will read in match files from `args.input`

    Notes
    -----
    will save a single called csv file will all calls for this run

    Returns
    -------
    call_file: str
        the absolute path to the merged calling file csv
    """

    # if an input was given use that to determine all the match_files, others read from the args
    if input_ is None:
        _match_files = args.input
    else:
        _match_files = input_

    # loop through each match fill and launch a calling sub-job for it
    _call_files = []
    _processes = []
    for i, _match_file in enumerate(_match_files):
        # create calling file name
        _call_file_loc = os.path.join(args.outdir, args.prefix + f"call_{_match_file.split('_')[-1].split('.')[0]}.csv")
        _call_files.append(_call_file_loc)
        # log some matching info
        try:
            _match_info = json.load(open(_match_file.replace(".pkl", ".json"), "r"))
            report.num_matches_found += _match_info['_num_single_matches'] + _match_info['_num_multiple_matches']
            logger.debug(f"subprocess {i} found "
                         f"{_match_info['_num_single_matches'] + _match_info['_num_multiple_matches']} matches")
        except FileNotFoundError:
            logger.warning("cannot find matching info.json for corresponding match.pkl")
        # launch sub job
        # goes match_file, library(s) as comma sep list, index(s) as comma sep list, strict as 0/1, outfile, make_id
        p = subprocess.Popen([sys.executable, str(os.path.join(MODULE_DIR, "calling_subprocess.py")), _match_file,
                              ",".join(args.library), ",".join(args.index), str(int(args.strict)),
                              _call_files[-1], args.make])
        _processes.append(p)
        logger.debug(f"launched calling subprocess {p.pid}")

    logger.info(f"found {report.num_matches_found} matches")

    # wait for all sub-jobs to terminate
    _p_bar = tqdm(total=report.num_matches_found, desc=f"calling sequences", disable=not args.print)
    _pbar_count = 0
    while _processes:
        time.sleep(0.5)
        if _processes[-1].poll() is None:  # check if last job is list of done
            _tmp_count = 0
            # collect current progress info for all jobs
            for file in [f for f in os.listdir(args.outdir) if f.endswith(".count")]:
                _file = os.path.join(args.outdir, file)
                try:
                    _tmp_count += int(open(_file, "r").readline().strip())
                except ValueError:
                    continue
            # update progress
            _p_bar.update(_tmp_count - _pbar_count)
            _pbar_count = _tmp_count
            continue
        else:
            _processes.pop()  # remove job from list if it is complete
    _p_bar.update(report.num_matches_found)  # finish the tqdm bar for aesthetics
    _p_bar.close()

    # remove count files
    for file in [f for f in os.listdir(args.outdir) if f.endswith(".count")]:
        os.remove(os.path.join(args.outdir, file))

    logger.debug("all calling subprocesses complete")

    # save merged match file if wanted
    _match_output_file = os.path.join(args.outdir, args.prefix + "_matches.pkl")
    if args.save_match:
        data = []
        logger.info("saving matches to {}".format(_match_output_file))
        # merge all the match files to 1 file for simplicity
        for _match_file in _match_files:
            data.extend(pickle.load(open(_match_file, "rb")))
        pickle.dump(data, open(_match_output_file, "wb"))

    # remove the tmp match files
    for _match_file in _match_files:
        # os.remove(_match_file)
        try:
            os.remove(_match_file.replace(".pkl", ".json"))
        except FileNotFoundError:
            pass

    # variables for reporting
    _num_called = 0
    _num_error = 0
    _num_failed = 0
    _num_passed = 0

    # merge call files
    _call_out_file = os.path.join(args.outdir, args.prefix + "_calls.csv")
    _call_out_file_object = open(_call_out_file, "w")
    _call_out_file_object.write(FILE_HEADER + "\n")
    # write all the calls to single call csv
    for i, _call_file in enumerate(_call_files):
        with open(_call_file, "r") as f:
            _header = f.readline()  # ignore the header
            for line in f:
                _call_out_file_object.write(line)
        # log some info for reporting
        _info = json.load(open(_call_file.replace(".csv", ".json"), "r"))
        logger.debug(f"subprocess {i} attempted calling {_info['num_called']}; "
                     f"{_info['num_passed']} passed, {_info['num_failed']} failed calling")

        _num_called += _info['num_called']
        _num_passed += _info['num_passed']
        _num_error += _info['num_error']
        _num_failed += _info['num_failed']
        # clean up
        os.remove(_call_file)  # remove file
        os.remove(_call_file.replace(".csv", ".json"))  # remove info file

    # log some more info for reporting
    logger.info(f"attempted {_num_called} calls; num called successfully: {_num_passed}; calls failed: {_num_failed}")
    report.num_matches_called = _num_passed
    report.called_csv_path = _call_out_file

    return str(_call_out_file)


def run_cube_gen(args: argparse.Namespace, logger: logging.Logger, report: HTMLReport, input_: str = None):
    """
    This script is used to convert a generated called.csv file into cube files that can be used for data analysis

    args: Namespace
        the arguments used to define the run settings
    logger: Logger
        the logger object for logging
    report: HTMLReport
        the report object to track reporting
    input_: str, default = None
        the list of match files to use
        if not passed will read in called csv file from `args.input`

    Notes
    -----
    will save a cube file for each library ID containing every index and a cube file for each index containing every
    library.
    """
    import polars as pl
    import polars.selectors as cs

    from deli.umi import cluster_umis

    # read in the input
    if input_ is None:
        df = pl.read_csv(args.input).filter("OVERALL_PASS")
    else:
        df = pl.read_csv(input_).filter("OVERALL_PASS")

    report.num_matches_called = len(df)

    # determine which libraries and index are in the file
    found_libs = df.unique(subset="CALLED_LIB").select(pl.col("CALLED_LIB")).to_series().to_list()
    logger.debug(f"found libs: {len(found_libs)}")
    found_idxs = df.unique(subset="CALLED_INDEX").select(pl.col("CALLED_INDEX")).to_series().to_list()
    logger.debug(f"found indexes: {len(found_idxs)}")

    # extract only the libraries and indexes that are request or in the file
    if args.library != "all":
        libs_to_use = list(set(found_libs).intersection(set(args.library)))
        libs_missed = list((set(args.library)) - set(found_libs))
        if len(libs_missed) != 0:
            logger.warning(f"failed to find libraries {libs_missed}")
    else:
        libs_to_use = found_libs
    if len(libs_to_use) == 0:
        raise ValueError("no matching libraries found in dataset")

    if args.index != "all":
        idxs_to_use = list(set(found_idxs).intersection(set(args.index + [args.control])))
        idxs_missed = list((set(args.index)) - set(found_idxs))
        if len(idxs_missed) != 0:
            logger.warning(f"failed to find indexes {idxs_missed}")
    else:
        idxs_to_use = found_idxs
    if len(idxs_to_use) == 0:
        raise ValueError("no matching indexes found in dataset")

    report.index_used = idxs_to_use
    report.lib_used = libs_to_use

    if args.normalize and (args.control not in idxs_to_use):
        raise ValueError(f"cannot find control index {args.control}")

    # begin the loop to generate the library-based cube files
    cube_outfiles = []
    for lib in libs_to_use:
        # load in the enumerated library data
        enumerated_lib = pl.read_csv(os.path.join(MODULE_DIR, "data", "enumerated_libraries", f"{lib}_library.csv.gz"))
        logger.debug(f"read in enumerated library for {lib} from "
                     f'{os.path.join(MODULE_DIR, "data", "enumerated_libraries", f"{lib}_library.csv.gz")}')
        outfile = os.path.join(args.outdir, args.prefix + f"{lib}_cube.csv")
        cube_outfiles.append((outfile, lib))

        # do umi clustering if need be #TODO BROKEN???
        if args.umi_cluster:
            logger.debug("conducting UMI clustering")
            df_lib = (
                df.filter(pl.col("CALLED_LIB") == lib)
                .group_by(["CALLED_INDEX", "DEL_ID"])
                .agg(pl.col("UMI").count().cast(pl.Int32).alias("raw"),
                     pl.col("UMI").n_unique().cast(pl.Int32).alias("corrected"),
                     pl.col("UMI").map_elements(cluster_umis, return_dtype=pl.Int64).add(1).alias("clustered"))
                .pivot(index="DEL_ID", columns="CALLED_INDEX", values=["raw", "corrected", "clustered"])
                .fill_null(0)
            )
        else:
            df_lib = (
                df.filter(pl.col("CALLED_LIB") == lib)
                .group_by(["CALLED_INDEX", "DEL_ID"])
                .agg(pl.col("UMI").count().cast(pl.Int32).alias("raw"),
                     pl.col("UMI").n_unique().cast(pl.Int32).alias("corrected"))
                .pivot(index="DEL_ID", columns="CALLED_INDEX", values=["raw", "corrected"])
                .fill_null(0)
            )

        # subtract off the control value if requested
        if args.normalize:
            logger.debug("conducting NTC normalization")
            if args.umi_cluster:
                _col_pat = ("^(corrected).*$", f"corrected_CALLED_INDEX_{args.control}")
            else:
                _col_pat = ("^(clustered).*$", f"clustered_CALLED_INDEX_{args.control}")

            df_lib = (
                df_lib
                .with_columns([(pl.col(_col_pat[0]) - pl.col(_col_pat[1]))
                              .cast(pl.UInt32, strict=False)
                              .name.prefix('normed_')])
                .fill_null(0)
            )

        # sort the columns and add in the enumerated library SMILES
        df_lib = (
            df_lib.select(["DEL_ID"] + sorted(df_lib.columns[1:],
                                              key=lambda x: int(x.split("_")[-1].split("index")[-1])))
            .join(enumerated_lib, on="DEL_ID", how="left")
        )

        # collect monosynthon enrichment if requested
        if args.monosynthon:
            for _mono in ["ID_A", "ID_B", "ID_C"]:
                mono_df = (
                    df_lib.group_by(pl.col(_mono))
                    .agg(pl.col("^(raw).*$").cast(pl.Int32).sum().name.prefix(_mono.split("_")[-1] + "_"),
                         pl.col("^(corrected).*$").cast(pl.Int32).sum().name.prefix(_mono.split("_")[-1] + "_"),
                         pl.col("^(normed).*$").cast(pl.Int32).sum().name.prefix(_mono.split("_")[-1] + "_"),
                         pl.col("^(clustered).*$").cast(pl.Int32).sum().name.prefix(_mono.split("_")[-1] + "_"))
                )
                df_lib = df_lib.join(mono_df, on=_mono, how="outer")
        # collect disynthon enrichment if requested
        if args.disynthon:
            df_lib = df_lib.with_columns(pl.col("DEL_ID").str.replace(r"-C\d\d\d", "").alias("Disynthon_AB"),
                                         pl.col("DEL_ID").str.replace(r"-B\d\d\d", "").alias("Disynthon_AC"),
                                         pl.col("DEL_ID").str.replace(r"A\d\d\d-", "").alias("Disynthon_BC"))
            for _id in ["Disynthon_AB", "Disynthon_AC", "Disynthon_BC"]:
                di_df = (
                    df_lib.group_by(pl.col(_id))
                    .agg(pl.col("^(raw).*$").cast(pl.Int32).sum().name.prefix(_id.split("_")[-1] + "_"),
                         pl.col("^(corrected).*$").cast(pl.Int32).sum().name.prefix(_id.split("_")[-1] + "_"),
                         pl.col("^(normed).*$").cast(pl.Int32).sum().name.prefix(_id.split("_")[-1] + "_"),
                         pl.col("^(clustered).*$").cast(pl.Int32).sum().name.prefix(_id.split("_")[-1] + "_"))
                )

                df_lib = df_lib.join(di_df, on=_id, how="outer")

        # rename cols to clean things up for future analysis
        (
            df_lib.rename({key: key.replace("_CALLED_INDEX", "") for key in df_lib.columns})
            .write_csv(outfile)
        )

        logger.info(f"generated cube for library {lib} with indexes {idxs_to_use} at {outfile}")

    # rename the numerical index to the experiment name using the user provided json configuration
    if args.index_name is not None:
        _convert = json.load(open(args.index_name, "r"))

        def _rename(name):
            # helper func to handle renaming cols
            for key, val in _convert.items():
                if name.endswith(key):
                    return name.replace(key, val)
            return name
    else:
        _convert = None

        def _rename(name):
            # dummy func to do nothing if no renaming is needed
            return name

    # loop through each index and
    for idx in idxs_to_use:
        if idx == args.control:
            continue

        # make the regex statement to match the desired index columns
        _match_regex = idx if args.control is None else f"({idx}|{args.control})"

        n = idx if _convert is None else _convert.get(idx, idx)  # convert index name if convertable

        # loop through all exist library cube generated earlier
        _process_libraries = []
        for _df, _lib in cube_outfiles:
            _df = pl.read_csv(_df)
            # skip the index if it turns out nothing with that index was called
            if f"raw_{idx}" not in _df.columns:
                continue
            _process_libraries.append(
                _df.select(~cs.contains("index") | cs.matches(f"^.*{_match_regex}$"))
                .with_columns(Library=pl.lit(_lib))
                .filter(pl.col(f"raw_{idx}") > 0)
                .select(pl.all().name.map(_rename))
                .select(~cs.matches("^*._right$"))
                .with_columns(DEL_ID=f"{_lib}-" + pl.col("DEL_ID"))
            )

        pl.concat(_process_libraries).write_csv(str(os.path.join(args.outdir, args.prefix + f"{n}_cube.csv")))
