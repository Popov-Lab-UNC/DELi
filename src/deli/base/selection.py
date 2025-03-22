"""defines base DEL selection objects"""

import abc
import random
import time
from os import PathLike
from typing import Self

from deli.dels import DELibraryPool


class BaseSelection(abc.ABC):
    """
    Base class for all children that deal with selection data/settings

    This can go beyond just a "selection".
    For example, all decoding and analysis experiment definitions are
    children of selection, as they are a selection with extra setting.
    """

    def __init__(
        self,
        library_pool: DELibraryPool,
        target_id: str,
        selection_id: str | None = None,
    ):
        """
        Initialize the BaseSelection object

        Notes
        -----
        DELi assumes all data is already demultiplex at the target level
        so every selection should only have one target ID

        DELi assumes uniqueness of the target ID, but it is never
        strictly enforced. It would not impact anything with decoding
        but it could mess with some downstream analysis

        Parameters
        ----------
        library_pool: DELibraryPool
            the libraries used in the selection
        target_id: str
            the id of the target used in the selection
        selection_id: str, optional
            the id of the selection.
            if `None` will default to a random number based on
            timestamp the object was created
        """
        self.library_pool = library_pool
        self.target_id = target_id
        self.selection_id: str = (
            selection_id
            if selection_id
            else (str(time.time()).replace(".", "") + f"{random.randint(0, 1000000):06d}")
        )

    @abc.abstractmethod
    def to_file(self, out_path: str | PathLike | None) -> PathLike:
        """Covert the object to a human readable file"""
        raise NotImplementedError()

    @classmethod
    @abc.abstractmethod
    def load_file(cls, path: str | PathLike | None) -> Self:
        """Load the object from a human readable file"""
        raise NotImplementedError()
