"""Hold base level classes for DEL compounds"""


class Library:
    """Base class for all libraries"""

    def __init__(self, library_id: str):
        self.library_id: str = library_id

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
