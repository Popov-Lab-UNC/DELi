"""Hold base level classes for DEL compounds"""


class _Library:
    """Base class for all libraries"""

    def __init__(self, library_id: str):
        self.library_id: str = library_id
