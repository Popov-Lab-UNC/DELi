"""
The DELi data directory is a filesystem that contains user specified DEL data.

The data directory is expected to have the following structure:
    <data_dir>/
        libraries/
            <library_id>.json
        building_blocks/
            <building_block_set_id>.csv
        tool_compounds/
            <tool_compound_id>.json
        reactions/
            <reaction_id>.json

Connection to the direcotry is file system agnositc via fsspec, so the data directory can
be local or remote (e.g. S3, GCS, etc.) as long as the correct filesystem implementation is used.
The default is to build use a local filesystem, but this can be configured by the user via the DELi config.
"""
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Final, Literal

from fsspec import AbstractFileSystem
from fsspec.implementations.local import LocalFileSystem


DELI_DATA_SUB_DIRS: Final = (
    ("libraries", ".json"),
    ("building_blocks", ".csv"),
    ("reactions", ".txt"),
    ("tool_compounds", ".json")
)


class DELiDataDirError(Exception):
    """Base class for errors related to the DELi data directory"""

    pass


class DELiDataDir:
    """Class representing the DELi data directory, which is a filesystem containing user specified DEL data."""

    def __init__(self, filesystem: AbstractFileSystem, validate: Literal["skip", "weak", "strict"] = "weak"):
        """
        Initialize the DELiDataDir object.

        Parameters
        ----------
        filesystem: AbstractFileSystem
            The fsspec filesystem to use for connecting to the data directory.
        validate: Literal["skip", "weak", "strict"], default "weak"
            Validation level for the data directory.
            "skip": skip all validation checks
            "weak": check and warn about file conflicts, but do not check file contents.
            "strict": check and fail if any file conflicts are detected.
        """
        self._filesystem = filesystem

        # check connection
        try:
            self._filesystem.ls("/")
        except Exception as e:
            raise DELiDataDirError("Unable to connect to DELi data directory with provided filesystem") from e

        # initalize and validate
        self._initalize_subdirs()

        if validate != "skip":
            conflicts = self.detect_file_conflicts()
            if conflicts:
                if validate == "weak":
                    warnings.warn(
                        f"File conflicts detected in DELi data directory sub-dirs {conflicts.keys()}",
                        UserWarning, stacklevel=1
                    )
                else:
                    raise DELiDataDirError(
                        f"File conflicts detected in DELi data directory sub-dirs: {conflicts.keys()}"
                    )

    def _initalize_subdirs(self) -> None:
        """Initialize the expected subdirectories in the data directory if they do not already exist."""
        for sub_dir in DELI_DATA_SUB_DIRS:
            if not self._filesystem.exists(sub_dir):
                self._filesystem.mkdir(sub_dir)
            else:
                if not self._filesystem.isdir(sub_dir):
                    raise DELiDataDirError(f"Expected subdirectory '{sub_dir}' is not a directory")

    def detect_file_conflicts(self) -> dict[str, dict[str, list[str]]]:
        """Detect if there is a conflict with an existing file in the data directory."""
        conflicting_files: dict[str, dict[str, list[str]]] = {}  # subdir -> file_stem -> list of conflicting files
        for sub_dir, file_ext in DELI_DATA_SUB_DIRS:
            all_files = self._filesystem.glob(f"{sub_dir}/**/*{file_ext}")
            file_stems = [Path(f).stem for f in all_files]

            if not len(set(file_stems)) == len(file_stems):
                _counter: defaultdict[str, list[str]] = defaultdict(list)
                for f, f_stem in zip(all_files, file_stems, strict=True):
                    _counter[f_stem].append(f)
                conflicting_stems = {stem: paths for stem, paths in _counter.items() if len(paths) > 1}
                conflicting_files[sub_dir] = conflicting_stems
        return conflicting_files

    def read_subdir_file_by_id(self, sub_dir: str, file_id: str) -> str:
        """Get the content of a file in the data directory."""
        if sub_dir not in [sd[0] for sd in DELI_DATA_SUB_DIRS]:
            raise ValueError(
                f"Invalid subdirectory '{sub_dir}'. Expected one of {[sd[0] for sd in DELI_DATA_SUB_DIRS]}"
            )
        matches = self._filesystem.glob(f"{sub_dir}/**/{file_id}")
        if not matches:
            raise DELiDataDirError(f"File '{file_id}' not found in subdirectory '{sub_dir}'")
        elif len(matches) > 1:
            raise DELiDataDirError(
                f"Multiple files found with ID '{file_id}' in subdirectory '{sub_dir}': {matches}. "
                "This indicates a file conflict in the data directory."
            )
        path = matches[0]
        return self._filesystem.read_text(path)

    def read_library_file(self, library_id: str) -> str:
        """Get the content of a library file in the data directory by library ID."""
        return self.read_subdir_file_by_id("libraries", library_id)

    def read_building_block_set_file(self, building_block_set_id: str) -> str:
        """Get the content of a building block set file in the data directory by building block set ID."""
        return self.read_subdir_file_by_id("building_blocks", building_block_set_id)

    def read_reaction_file(self, reaction_id: str) -> str:
        """Get the content of a reaction file in the data directory by reaction ID."""
        return self.read_subdir_file_by_id("reactions", reaction_id)

    def read_tool_compound_file(self, tool_compound_id: str) -> str:
        """Get the content of a tool compound file in the data directory by tool compound ID."""
        return self.read_subdir_file_by_id("tool_compounds", tool_compound_id)


def _parse_filesystem_config(filesystem_type: str, filesystem_config: dict) -> AbstractFileSystem:
    # New file systems can be added here
    if filesystem_type == "local":
        return _parse_local_filesystem(filesystem_config)
    else:
        raise DELiDataDirError(f"Unsupported filesystem type: {filesystem_type}")


def _parse_local_filesystem(filesystem_config: dict) -> LocalFileSystem:
    """Parse a local filesystem config and return a LocalFileSystem object."""
    root_dir = filesystem_config.get("dir", None)
    if root_dir is None:
        raise DELiDataDirError("Missing 'dir' key in local filesystem config")
    return LocalFileSystem(root_dir=root_dir)
