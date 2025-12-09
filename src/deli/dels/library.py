"""Hold base level classes for DEL compounds"""
import abc
from typing import TYPE_CHECKING, Generic, TypeVar, Sequence, Iterator


if TYPE_CHECKING:
    from deli.dels.compound import Compound


C = TypeVar("C", bound=Compound)
class Library(Generic[C], abc.ABC):
    """Base class for all libraries"""
    library_size: int

    def __init__(self, library_id: str):
        self.library_id: str = library_id

    def len(self) -> int:
        """Return the size of the library"""
        return self.library_size

    def __hash__(self):
        """Return the hash of the object."""
        return hash(self.library_id)

    def __eq__(self, other):
        """Check if two Library objects are equal (share the same library_id)"""
        if not isinstance(other, Library):  # any child of this class can be equal
            return False
        return self.library_id == other.library_id

    def __repr__(self):
        """Return a string representation of the object."""
        return f"{self.__class__.__name__}({self.library_id})"

    @abc.abstractmethod
    def get_compound(self, *args, **kwargs) -> C:
        """
        Given some type of compound id input, return that compound object.

        Notes
        -----
        Not all libraries will use the same type of input to identify compounds.
        For example, some libraries may use a simple string compound ID,
        while others may require a list of building block IDs.
        The specific parameters for this method should be defined in the subclasses.


        Returns
        -------
        Compound
            The corresponding Compound object.

        Raises
        ------
        KeyError
            If the compound ID is not found in the library.
        """

        raise NotImplementedError()


LibType = TypeVar("LibType", bound=Library)
class LibraryCollection(Generic[LibType]):
    """
    base class for any class that holds a group of DEL libraries
    """

    def __init__(self, libraries: Sequence[LibType]):
        """
        Initialize a DELibrarySchemaGroup object

        Parameters
        ----------
        libraries: List[CombinatorialLibrary]
            libraries to include in the library schema group
        """
        self.libraries: Sequence[LibType] = libraries
        self._library_map = {lib.library_id: lib for lib in self.libraries}

        self.collection_size = sum([lib.library_size for lib in self.libraries])

        ### VALIDATE ###
        _ids: list[str] = []
        for _library in self.libraries:
            # check id uniqueness
            if _library.library_id in _ids:
                raise KeyError(f"multiple libraries share identical `library_id` '{_library.library_id}'")
            else:
                _ids.append(_library.library_id)

    def __len__(self) -> int:
        """Return the number of libraries in the library collection"""
        return len(self.libraries)

    def __iter__(self) -> Iterator[LibType]:
        """Iterate through all libraries in the library collection"""
        return iter(self.libraries)

    def __contains__(self, item) -> bool:
        """Check if a library with the given ID is in the collection"""
        return item in self._library_map

    def __getitem__(self, item) -> LibType:
        """Get a library by its ID from the collection"""
        return self.get_library(item)

    def get_library(self, library_id: str) -> LibType:
        """
        Return the library from the collection with the same ID

        Parameters
        ----------
        library_id: str
            id of the library to get

        Returns
        -------
        CombinatorialLibrary

        Raises
        ------
        KeyError
            if `library_id` not in the collection
        """
        try:
            return self._library_map[library_id]
        except KeyError as e:
            raise KeyError(KeyError(f"cannot find library with id '{library_id}' in collection")) from e
