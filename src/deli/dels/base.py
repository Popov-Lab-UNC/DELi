"""Hold generic mixin and abstract classes for use by DEL objects"""

import abc

from deli.configure import accept_deli_data_name


class DeliDataLoadableMixin(abc.ABC):
    """Mixin for objects that can be loaded from DeliDataDir"""

    @classmethod
    @abc.abstractmethod
    @accept_deli_data_name(sub_dir="barcodes", extension="json")
    def load(cls, name_or_path: str):
        """Load the file into the object"""
        raise NotImplementedError()
