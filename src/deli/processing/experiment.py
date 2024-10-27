"""handles DEL decoding experiments definition"""

from os import PathLike
from typing import Dict, List, Optional, Self, Union

from deli.dels import BarcodeSchema, DELibrary, DELibraryGroup, IndexSet


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

    def _check_validity(self):
        """Check that a Primer DEL experiment is valid"""
        super()._check_validity()
        _lib_1 = self.libraries[0]
        for _library in self.libraries[1:]:
            if not _library.barcode_schema.is_experiment_compatible(_lib_1.barcode_schema):
                raise DELExperimentError(
                    "all libraries in primer experiment must have " "the same primer regions"
                )


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
                info[_current_section].append(line)

        # load in the libraries
        if "library" not in info.keys():
            raise DELExperimentError(
                f"cannot find the 'library' section in the experiment file: {path}"
            )
        _libs: list[DELibrary] = list()
        for _library in info["library"]:
            _libs.append(DELibrary.load(_library))
        library = DELibraryGroup(_libs)

        # load in the index if included
        index = IndexSet(list())
        if "index" in info.keys():
            for index_set in info["index"]:
                index += IndexSet.load(index_set)

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
        _sub_experiments: Dict[BarcodeSchema, List[DELibrary]] = dict()
        for _lib in self.libraries:
            _found_sub_experiment: bool = False
            for existing_sub_experiment in _sub_experiments.keys():
                if _lib.barcode_schema.is_experiment_compatible(existing_sub_experiment):
                    _sub_experiments[existing_sub_experiment].append(_lib)
                    _found_sub_experiment = True
                    break
            if not _found_sub_experiment:
                _sub_experiments[_lib.barcode_schema] = [_lib]

        return [
            PrimerDELExperiment(DELibraryGroup(_lib), self.indexes)
            for _lib in _sub_experiments.values()
        ]
