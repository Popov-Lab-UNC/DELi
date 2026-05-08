"""Linkers for DELs"""

from typing import Generic, Sequence, TypeVar

from deli.utils import SmilesMixin

from ._base import DecodableObject, DecodableObjectSet, TaggedDecodableObject, TaggedDecodableObjectSet


class Linker(SmilesMixin, DecodableObject):
    """
    A linker that connects a small molecule to a DNA tag in a DEL.

    Parameters
    ----------
    obj_id : str
        The ID of the linker.
    smiles : str | None, optional
        The SMILES string representation of the linker.
    truncated_smiles : str | None, optional
        A truncated SMILES string representation of the linker.
        Used to generate more small molecule like compounds during library enumeration.
    """

    def __init__(self, obj_id: str, smiles: str | None, truncated_smiles: str | None = None, **kwargs):
        super().__init__(smiles=smiles, obj_id=obj_id)
        self.truncated_smiles: str | None = truncated_smiles


class TaggedLinker(Linker, TaggedDecodableObject):
    """
    A linker that connects a small molecule to a DNA tag in a DEL, along with its associated tag information.

    Parameters
    ----------
    linker : Linker
        The linker object.
    tag_info : dict
        A dictionary containing the tag information associated with the linker.
    """

    def __init__(self, obj_id: str, smiles: str | None, truncated_smiles: str | None, tags: list[str] | str):
        super().__init__(obj_id=obj_id, smiles=smiles, truncated_smiles=truncated_smiles, tags=tags)


T = TypeVar("T", bound=Linker)


class LinkerSet(DecodableObjectSet[T], Generic[T]):
    """
    A set of linkers that can be used in a DEL.

    Parameters
    ----------
    objects : list[Linker]
        A list of Linker objects to include in the set.
    """

    def __init__(self, linkers: Sequence[T], **kwargs):
        super().__init__(objects=linkers, **kwargs)


class TaggedLinkerSet(LinkerSet[TaggedLinker], TaggedDecodableObjectSet[TaggedLinker]):
    """
    A set of linkers that can be used in a DEL, along with their associated tag information.

    Parameters
    ----------
    objects : list[TaggedLinker]
        A list of TaggedLinker objects to include in the set.
    """

    def __init__(self, linkers: Sequence[TaggedLinker], **kwargs):
        super().__init__(linkers=linkers, **kwargs)
