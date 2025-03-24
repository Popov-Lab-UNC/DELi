"""base classes for setting objects"""

import abc
from os import PathLike
from typing import Self


class BaseSettings(abc.ABC):
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

    @abc.abstractmethod
    def to_file(self, path: str | PathLike):
        """Save settings to a yaml file"""
        raise NotImplementedError()

    @classmethod
    @abc.abstractmethod
    def from_file(cls, path: str | PathLike) -> Self:
        """Load settings from a yaml file"""
        raise NotImplementedError()
