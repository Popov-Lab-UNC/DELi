from Levenshtein import distance

from deli.constants import INDEX, LIBS, LIBS_TO_MAKE, DeliConfigError, LIBRARY_RISK_NORMAL_THRESHOLD, \
    INDEX_RISK_NORMAL_THRESHOLD


def lib_to_make(lib: str) -> str:
    """
    Determines the "make" of the passed DEL library (if it exists, otherwise raises and exception)

    Parameters
    ----------
    lib: str
        the library ID to get the "make" for

    Returns
    -------
    make: str
        the make ID of the library

    Raises
    ------
    DeliConfigError
        "cannot find 'make' for library {lib} in DELi configs"
    """
    _tmp = LIBS_TO_MAKE.get(lib, None)
    if _tmp is None:
        raise DeliConfigError(f"cannot find 'make' for library {lib} in DELi configs")
    return _tmp


def make_barcode_pattern(barcode_length) -> dict:
    """
    This function will take a dictionary defining the lengths of each DEL barcode section and convert it
    into a dictionary containing the start and stop index of the barcode sections (starting from 0)

    Notes
    -----
    Assume that the dictionary is ordered, such that the order of the keys is the order of the barcode sections

    Parameters
    ----------
    barcode_length: dict
        dictionary defining the lengths of each DEL barcode section

    Returns
    -------
    barcode_pattern: dict
        dictionary of the (inclusive) start and (non-inclusive) stop index of each section: `dict[key] = [start, stop)`

    """
    barcode_pattern = {}
    running_total = 0
    for barcode_section_key, val in barcode_length.items():
        barcode_pattern[barcode_section_key] = (running_total, val + running_total)
        running_total += val
    return barcode_pattern


def get_min_index_distance(included_index: list[str] = None) -> int:
    """
    Given a list of index IDs, determine the minimum Levenshtein distance between all of them

    Parameters
    ----------
    included_index: list[str] default = None
        the index_ids to use
        if left as None will return 0

    Returns
    -------
    min_distance: int
        the minimum Levenshtein distance between the passed indexes
    """
    if len(included_index) == 1:
        return INDEX_RISK_NORMAL_THRESHOLD * 2  # if only one index just give it a constant threshold
    if included_index is None:
        return 0
    _index_sequences = [INDEX[i] for i in included_index]
    return min([distance(s1, s2) for s1 in _index_sequences for s2 in _index_sequences if s1 != s2])


def get_dist_from_index(seq: str, included_index: list[str] = None) -> dict[str: int]:
    """
    Given a sequence and list of index IDs, determines the Levenshtein distance between each index and the sequence

    Parameters
    ----------
    seq: str
        the DNA sequence to compare each index to
    included_index: list[str], default = None
        a list of the index IDs to get distances for
        if not passed assumes all index ids in DELi config need to be used

    Returns
    -------
    dict[str: int]:
        a map of index ID to distance from the passed sequence
    """

    if included_index is None:
        return {key: distance(val, seq) for key, val in INDEX.items()}
    picked_index = {idx: val for idx, val in INDEX.items() if idx in included_index}
    return {key: distance(val, seq) for key, val in picked_index.items()}


def get_min_lib_distance(included_lib: list[str] = None) -> int:
    """
    Given a list of library IDs, determine the minimum Levenshtein distance between all of them

    Parameters
    ----------
    included_lib: list[str] default = None
        the index_ids to use
        if left as None will return 0

    Returns
    -------
    min_distance: int
        the minimum Levenshtein distance between the passed libraries
    """
    if len(included_lib) == 1:
        return LIBRARY_RISK_NORMAL_THRESHOLD * 2  # if only one index just give it a constant threshold
    if included_lib is None:
        return 0
    _lib_sequences = [LIBS[i] for i in included_lib]
    return min([distance(s1, s2) for s1 in _lib_sequences for s2 in _lib_sequences if s1 != s2])


def get_dist_from_lib(seq: str, included_libs: list[str] = None) -> dict[str: int]:
    """
    Given a sequence and list of library IDs, determines the Levenshtein distance between each library and the sequence

    Parameters
    ----------
    seq: str
        the DNA sequence to compare each index to
    included_libs: list[str], default = None
        a list of the library IDs to get distances for
        if not passed assumes all library ids in DELi config need to be used

    Returns
    -------
    dict[str: int]:
        a map of library ID to distance from the passed sequence
    """
    if included_libs is None:
        return {key: distance(val, seq) for key, val in LIBS.items()}
    picked_libs = {lib: val for lib, val in LIBS.items() if lib in included_libs}
    return {key: distance(val, seq) for key, val in picked_libs.items()}
