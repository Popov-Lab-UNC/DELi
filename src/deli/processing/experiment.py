"""handles DEL decoding experiments definition"""

from os import PathLike
from typing import Dict, List, Optional, Self, Tuple, Union

from deli.dels import BarcodeSchema, DELibrary, DELibraryGroup, Index, IndexSet


class DELExperimentError(Exception):
    """raised when there is an issue with a DEL experiment"""

    pass


class _BaseExperiment:
    def __init__(self, libraries: DELibraryGroup, indexes: Optional[IndexSet] = None):
        """
        Initialize the DELExperiment object

        Notes
        -----
        DEL experiments hold the libraries (and indexes) that were used in a DEL run
        This info is needed both for decoding DEL sequences and analysis of the DELs.

        When

        Parameters
        ----------
        libraries
        indexes
        """
        self.libraries = libraries
        if indexes is None:
            indexes = IndexSet([])
        self.indexes = indexes
        self._check_validity()

    def _check_validity(self):
        """Check that a DEL experiment is valid"""
        if len(self.libraries) < 1:
            raise DELExperimentError("DELExperiment must have at least one library")

        # check for too many indexes
        if len(self.indexes) > 1 and all(
            [not lib.barcode_schema.has_index() for lib in self.libraries]
        ):
            raise DELExperimentError(
                f"barcode schema(s) lacks index section but multiple indexes "
                f"were passed to the experiment: {self.indexes}"
            )

        # check for too many libraries
        if len(self.libraries) > 1 and all(
            [not lib.barcode_schema.has_index() for lib in self.libraries]
        ):
            raise DELExperimentError(
                f"barcode schema(s) lacks library_tag section but multiple libraries "
                f"were passed to the experiment: {self.libraries}"
            )


class PrimerDELExperiment(_BaseExperiment):
    """A DEL experiment where all libraries are have the same primer"""

    def __init__(self, libraries: DELibraryGroup, indexes: Optional[IndexSet] = None):
        """
        Initialize a PrimerDELExperiment object

        Notes
        -----
        Primer DEL Experiments are what a single matching run needs

        You should avoid initializing this class manually
        It should be created by calling `break_into_matching_experiments`
        on a DELExperiment object (that you create manually)

        Parameters
        ----------
        libraries: DELibraryGroup
            libraries in this experiment
        indexes
        """
        super().__init__(libraries=libraries, indexes=indexes)

        self.library_schema_groups = self.libraries.break_into_schema_groups()
        self.primer = self.libraries[0].barcode_schema["primer"]

    def _check_validity(self):
        """Check that a Primer DEL experiment is valid"""
        super()._check_validity()
        _lib_1 = self.libraries[0]
        for _library in self.libraries[1:]:
            if not _library.barcode_schema.is_experiment_compatible(_lib_1.barcode_schema):
                raise DELExperimentError(
                    "all libraries in primer experiment must have " "the same primer regions"
                )

    def get_min_max_in_front(self) -> Tuple[int, int]:
        """
        Get the minimum and maximum distance from the start to front of primer

        Returns
        -------
        min, max: Tuple[int, int]
        """
        _min_dist = float("inf")
        _max_dist = float("-inf")

        for _lib in self.libraries:
            _dist_from_front = _lib.barcode_schema.barcode_spans["primer"][0]

            if _dist_from_front < _min_dist:
                _min_dist = _dist_from_front
            if _dist_from_front > _max_dist:
                _max_dist = _dist_from_front
        return int(_min_dist), int(_max_dist)

    def get_min_max_behind(self) -> Tuple[int, int]:
        """
        Get the minimum and maximum distance from the end to back of primer

        Returns
        -------
        min, max: Tuple[int, int]
        """
        _min_dist = float("inf")
        _max_dist = float("-inf")

        for _lib in self.libraries:
            _front = _lib.barcode_schema.barcode_spans["primer"][1]

            _back = len(_lib.barcode_schema.full_barcode)

            _dist_from_front = _back - _front

            if _dist_from_front < _min_dist:
                _min_dist = _dist_from_front
            if _dist_from_front > _max_dist:
                _max_dist = _dist_from_front
        return int(_min_dist), int(_max_dist)


class DELExperiment(_BaseExperiment):
    """Class to hold info on a given DEL experiment"""

    @classmethod
    def load_experiment(cls, path: Union[str, PathLike]) -> Self:
        """
        Load a DEL experiment from a file

        Parameters
        ----------
        path: str or PathLike
            the path the experiment file

        Returns
        -------
        Self
        """
        # read the file
        info: Dict[str, List[str]] = {}
        with open(path, "r") as file:
            _current_section: str = ""
            for line in file:
                _line = line.strip()
                if _line.startswith("[") and _line.endswith("]"):
                    _current_section = _line[1:-1]
                    info[_current_section] = list()
                    continue
                if _line != "":
                    info[_current_section].append(_line)

        # load in the library groups first
        library_group = DELibraryGroup(list())
        if "library_group" in info.keys():
            for library_group_name in info["library_group"]:
                library_group += DELibraryGroup.load(library_group_name)

        # load in the libraries
        _libs: list[DELibrary] = list()
        if "libraries" in info.keys():
            for _library in info["libraries"]:
                if library_group.has_library_with_name(_library):
                    _libs.append(library_group.get_library_with_name(_library))
                else:
                    _libs.append(DELibrary.load(_library))
            library = DELibraryGroup(_libs)
        else:
            if len(library_group) == 0:
                raise DELExperimentError(
                    f"cannot find the 'library_group' or 'libraries' "
                    f"section in the experiment file: {path}"
                )
            library = library_group

        # load in the index set
        index_set = IndexSet(list())
        if "index_set" in info.keys():
            for index_set_name in info["index_set"]:
                index_set += IndexSet.load(index_set_name)

        # load in the libraries
        _indexes: list[Index] = list()
        if "indexes" in info.keys():
            for _index_info in info["indexes"]:
                splits = _index_info.split(":")
                _sample_name = splits[0].strip()
                _sample_indexes = [_idx.strip() for _idx in splits[1].split(",")]

                if len(_sample_indexes) > 1:
                    _sample_names = [
                        _sample_name + f"_rep{i+1}" for i in range(len(_sample_indexes))
                    ]
                else:
                    _sample_names = [_sample_name]

                for _samp_name, _index in zip(_sample_names, _sample_indexes):
                    if index_set.has_index_with_name(_index):
                        _tmp = index_set.get_index_with_name(_index)
                    else:
                        _tmp = Index.load(_index)
                    _tmp.sample_name = _samp_name
                    _indexes.append(_tmp)
            index = IndexSet(_indexes)
        else:
            if len(index_set) > 0:
                raise DELExperimentError(
                    "must manually map indexes to sample names; "
                    "see 'Defining DEL experiments for more info'"
                )
            index = index_set

        return cls(libraries=library, indexes=index)

    def break_into_matching_experiments(self) -> List[PrimerDELExperiment]:
        """
        Breaks the DELExperiment into the sub PrimerDELExperiments

        Notes
        -----
        PrimerDELExperiments are the library sets that exist
        within the DELExperiment that have the same primer region
        of their barcode

        For example if I had 5 libraries, 2 of which have primer
        ZZZZ and 3 had ZZZZZZZZZZZZ then we have 2 PrimerDELExperiment
        groups to break into. A PrimerDELExperiment has all the
        indexes that its parent DELExperiment does

        Returns
        -------
        List[PrimerDELExperiment]
        """
        _sub_experiments: List[Tuple[BarcodeSchema, List[DELibrary]]] = list()
        for _lib in self.libraries:
            _found_sub_experiment: bool = False
            for i, (existing_sub_experiment, _) in enumerate(_sub_experiments):
                if _lib.barcode_schema.is_experiment_compatible(existing_sub_experiment):
                    _sub_experiments[i][1].append(_lib)
                    _found_sub_experiment = True
                    break
            if not _found_sub_experiment:
                _sub_experiments.append(
                    (
                        _lib.barcode_schema,
                        [_lib],
                    )
                )

        return [
            PrimerDELExperiment(DELibraryGroup(_lib[1]), self.indexes) for _lib in _sub_experiments
        ]
