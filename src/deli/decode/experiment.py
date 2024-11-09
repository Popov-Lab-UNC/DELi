"""handles DEL decoding experiments definition"""

from os import PathLike
from typing import Dict, List, Optional, Self, Tuple, Union

from deli.dels import BarcodeSchema, DELibrary, DELibraryGroup, Index, IndexSet
from deli.settings import CallerSettings, CubeGenSettings, MatcherSettings


class DELExperimentError(Exception):
    """raised when there is an issue with a DEL experiment"""

    pass


class _BaseExperiment:
    def __init__(
        self,
        libraries: DELibraryGroup,
        indexes: Optional[IndexSet] = None,
        matching_settings: Optional[MatcherSettings] = None,
        caller_settings: Optional[CallerSettings] = None,
        cube_settings: Optional[CubeGenSettings] = None,
    ):
        """
        Initialize the DELExperiment object

        Notes
        -----
        DEL experiments hold the libraries (and indexes) that were used in a DEL run
        This info is needed both for decoding DEL sequences and analysis of the DELs.

        Parameters
        ----------
        libraries: DELibraryGroup
            the libraries to be used for calling
        indexes: Optional[IndexSet]
            the indexes to be used for de-multiplexing (if turned on)
        matching_settings: Optional[MatcherSettings]
            the settings to use for matching
            if None will use default settings
        caller_settings: Optional[CallerSettings]
            the settings to use for calling
            if None will use default settings
        cube_settings: Optional[CubeGenSettings]
            the settings to use for CubeGen
            if None will use default settings
        """
        self.libraries = libraries
        if indexes is None:
            indexes = IndexSet([])
        self.indexes = indexes
        self.matching_settings = (
            matching_settings if matching_settings is not None else MatcherSettings()
        )
        self.caller_settings = caller_settings if caller_settings is not None else CallerSettings()
        self.cubegen_settings = cube_settings if cube_settings is not None else CubeGenSettings()
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

    def to_experiment_file(self, path: Union[str, PathLike]):
        with open(path, "w") as f:
            # libraries
            f.write("[libraries]\n")
            for lib in self.libraries:
                f.write(str(lib.loaded_from))
                f.write("\n")
            f.write("\n")

            # indexes
            f.write("[indexes]\n")
            _index_samples = dict()
            for idx in self.indexes:
                _samp_name = idx.sample_name
                if _samp_name == "" or _samp_name is None:
                    _index_samples[idx.index_id] = [str(idx.loaded_from)]
                else:
                    _index_samples[_samp_name] = [str(idx.loaded_from)]
            for samp_name, index_list in _index_samples.items():
                _str = f"{samp_name}: " + ",".join(index_list)
                f.write(_str)
                f.write("\n")
            f.write("\n")

            # settings
            f.write("[matching_settings]\n")
            f.write(
                "\n".join(
                    [f"{key}: {val}" for key, val in self.matching_settings.__dict__.items()]
                )
            )
            f.write("\n\n")

            f.write("[call_settings]\n")
            f.write(
                "\n".join([f"{key}: {val}" for key, val in self.caller_settings.__dict__.items()])
            )
            f.write("\n\n")

            f.write("[cube_settings]\n")
            f.write(
                "\n".join([f"{key}: {val}" for key, val in self.cubegen_settings.__dict__.items()])
            )
            f.write("\n\n")

    def get_checksum(self, path: Union[str, PathLike]):
        # TODO in issue #36
        raise NotImplementedError


class PrimerDELExperiment(_BaseExperiment):
    """A DEL experiment where all libraries are have the same primer"""

    def __init__(
        self,
        libraries: DELibraryGroup,
        indexes: Optional[IndexSet] = None,
        matching_settings: Optional[MatcherSettings] = None,
        caller_settings: Optional[CallerSettings] = None,
        cube_settings: Optional[CubeGenSettings] = None,
    ):
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
        indexes: Optional[IndexSet]
            the indexes to be used for de-multiplexing (if turned on)
        matching_settings: Optional[MatcherSettings]
            the settings to use for matching
            if None will use default settings
        caller_settings: Optional[CallerSettings]
            the settings to use for calling
            if None will use default settings
        cube_settings: Optional[CubeGenSettings]
            the settings to use for CubeGen
            if None will use default settings
        """
        super().__init__(
            libraries=libraries,
            indexes=indexes,
            caller_settings=caller_settings,
            matching_settings=matching_settings,
            cube_settings=cube_settings,
        )

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
        DELExperiment
        """
        # read the file
        info: Dict[str, List[Tuple[str, int]]] = {}
        with open(path, "r") as file:
            _current_section: str = ""
            for i, line in enumerate(file):
                _line = line.strip()
                if _line.startswith("[") and _line.endswith("]"):
                    _current_section = _line[1:-1]
                    info[_current_section] = list()
                    continue
                if _line != "":
                    info[_current_section].append((_line, i))

        # load in the libraries
        _libs: list[DELibrary] = list()
        if "libraries" in info.keys():
            for _library, _line_number in info["libraries"]:
                _libs.append(DELibrary.load(_library))
            library = DELibraryGroup(_libs)
        else:
            raise DELExperimentError(
                "requires section 'libraries' is missing from experiment file"
            )

        # load in the indexes
        _indexes: list[Index] = list()
        if "indexes" in info.keys():
            for _index_info, _line_number in info["indexes"]:
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
                    _tmp = Index.load(_index)
                    _tmp.sample_name = _samp_name
                    _indexes.append(_tmp)
            index = IndexSet(_indexes)
        else:
            index = IndexSet([])

        # load in optional matching settings
        _matching_settings = dict()
        if "matching_settings" in info.keys():
            for _match_setting, _line_number in info["matching_settings"]:
                _name, _val = _match_setting.replace(" ", "").split(":")
                _matching_settings[_name] = _val
        matching_settings = MatcherSettings.from_dict(_matching_settings)

        # load in optional calling settings
        _calling_settings = dict()
        if "call_settings" in info.keys():
            for _call_setting, _line_number in info["call_settings"]:
                _name, _val = _call_setting.replace(" ", "").split(":")
                _calling_settings[_name] = _val
        calling_settings = CallerSettings.from_dict(_calling_settings)

        # load in optional cubegen settings
        _cube_settings = dict()
        if "cube_settings" in info.keys():
            for _cube_setting, _line_number in info["cube_settings"]:
                _name, _val = _cube_setting.replace(" ", "").split(":")
                _cube_settings[_name] = _val
        cube_settings = CubeGenSettings.from_dict(_cube_settings)

        return cls(
            libraries=library,
            indexes=index,
            caller_settings=calling_settings,
            cube_settings=cube_settings,
            matching_settings=matching_settings,
        )

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
            PrimerDELExperiment(
                DELibraryGroup(_lib[1]),
                self.indexes,
                caller_settings=self.caller_settings,
                cube_settings=self.cubegen_settings,
                matching_settings=self.matching_settings,
            )
            for _lib in _sub_experiments
        ]
