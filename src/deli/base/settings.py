"""base classes for setting objects"""

from os import PathLike
from typing import Self

import yaml


class BaseSettings:
    """
    Base class for all setting objects

    Setting object are used to organize the possible
    settings that DELi CLI commands might use.

    DELi uses YAML syntax for all setting files and
    only supports safe loading of yaml files.
    This means all variables saved in settings
    must be primitive objects in python
    """

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def to_file(self, path: str | PathLike):
        """
        Save settings to a yaml file

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

        Parameters
        ----------
        path : str | PathLike
            Path to yaml file

        Returns
        -------
        Self
        """
        try:
            return cls(**yaml.safe_load(open(path, "r")))
        except Exception as e:
            raise RuntimeError(f"Failed to load settings from {path}") from e
