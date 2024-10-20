"""Hold generic mixin and abstract classes for use by DEL objects"""

import abc

from deli.configure import accept_deli_data


class DeliDataLoadableMixin(abc.ABC):
    """Mixin for objects that can be loaded from DeliDataDir"""

    @classmethod
    @accept_deli_data("barcodes", ".json")
    @abc.abstractmethod
    def load(cls, path: str):
        """
        Load in an object from Deli Data
        """
        raise NotImplementedError()
