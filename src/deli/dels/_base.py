"""
Private base classes for DEL objects.

These classes are only used inside the `deli.dels` module and are not intended
to be used directly by users of the library.
"""
from typing import Generic, Literal, Sequence, TypeVar, overload

from .exceptions import DecodableSetConflict


class DecodableObject:
    """Base class for objects that can be decoded from a dictionary."""

    def __init__(self, obj_id: str, **kwargs):
        super().__init__(**kwargs)
        self._obj_id = obj_id

    def __hash__(self):
        """Hash of the object is based on its ID"""
        return hash(self._obj_id)

    @property
    def id(self) -> str:
        """The ID of the object."""
        return self._obj_id


class TaggedDecodableObject(DecodableObject):
    """Base class for objects that can be decoded from a dictionary and have tags."""

    def __init__(self, obj_id: str, tags: list[str] | str, **kwargs):
        super().__init__(obj_id, **kwargs)
        self._tags = tags if isinstance(tags, list) else [tags]

    @property
    def tags(self) -> list[str]:
        """The tags associated with the object."""
        return self._tags


T = TypeVar("T", bound=DecodableObject)
TaggedT = TypeVar("TaggedT", bound=TaggedDecodableObject)

class DecodableObjectSet(Generic[T]):
    """Base class for sets of decodable objects."""

    def __init__(self, objects: Sequence[T], **kwargs):
        super().__init__(**kwargs)
        self._objects: tuple[T, ...] = tuple(objects)

        self._id_lookup: dict[str, T] = {}
        for obj in objects:
            if obj._obj_id in self._id_lookup:
                raise DecodableSetConflict(
                    f"ID '{obj._obj_id}' is already used by another {type(obj).__name__} in the set."
                )
            self._id_lookup[obj._obj_id] = obj

    def __len__(self):
        """Number of objects in the set"""
        return len(self._objects)

    def __iter__(self):
        """Iterate over the objects"""
        return iter(self._objects)

    @overload
    def get_obj_by_id(self, query_id: str, fail_on_missing: Literal[False]) -> T | None: ...

    @overload
    def get_obj_by_id(self, query_id: str, fail_on_missing: Literal[True]) -> T: ...

    @overload
    def get_obj_by_id(self, query_id: str, fail_on_missing: bool = False) -> T | None: ...

    def get_obj_by_id(self, query_id: str, fail_on_missing: bool = False) -> T | None:
        """
        Given an ID, search for corresponding object for that ID

        Notes
        -----
        Will return `None` if no matching object is found *unless* `fail_on_missing` is
        set to `True`, in which case a `KeyError` will be raised.

        Parameters
        ----------
        query_id: str
            ID to query
        fail_on_missing: bool, default False
            if `True` raise a KeyError is no match is found
            else return `None`

        Returns
        -------
        T or None
            will be `None` if no matching object is found

        Raises
        ------
        KeyError
            If `fail_on_missing` is `True` and no matching object is found
        """
        _obj = self._id_lookup.get(query_id, None)
        if _obj is None and fail_on_missing:
            raise KeyError(f"id '{query_id}' not found in {type(self).__name__}")
        return _obj


class TaggedDecodableObjectSet(DecodableObjectSet[TaggedT], Generic[TaggedT]):
    """Base class for sets of decodable objects that have tags."""

    def __init__(self, objects: list[TaggedT], **kwargs):
        super().__init__(objects, **kwargs)

        self.tag_length = len(self._objects[0].tags[0])
        for _obj in self._objects[1:]:
            for _tag in _obj.tags:
                if len(_tag) != self.tag_length:
                    raise DecodableSetConflict(
                        f"All tags in a {type(self).__name__} must have the same length; "
                        f"observed '{self.tag_length}' has length {len(_tag)}"
                    )

        self._tag_lookup_table: dict[str, TaggedT] = {}
        for obj in objects:
            for tag in obj._tags:
                if tag in self._tag_lookup_table:
                    if self._tag_lookup_table[tag].id != obj.id:
                        raise DecodableSetConflict(
                            f"Tag '{tag}' is already used by another {type(obj).__name__} "
                            f"(ID: {self._tag_lookup_table[tag].id}) in the set."
                        )
                self._tag_lookup_table[tag] = obj
        self.num_tags = len(self._tag_lookup_table)

    @property
    def tag_map(self) -> dict[str, TaggedT]:
        """Mapping of tags to objects in the set."""
        return self._tag_lookup_table

    @overload
    def get_obj_by_tag(self, query: str, fail_on_missing: Literal[False]) -> TaggedT | None: ...

    @overload
    def get_obj_by_tag(self, query: str, fail_on_missing: Literal[True]) -> TaggedT: ...

    @overload
    def get_obj_by_tag(self, query: str, fail_on_missing: bool = False) -> TaggedT | None: ...

    def get_obj_by_tag(self, query: str, fail_on_missing: bool = False) -> TaggedT | None:
        """
        Given a tag, search for corresponding object for that tag

        Notes
        -----
        Will return `None` if no matching object is found *unless* `fail_on_missing` is
        set to `True`, in which case a `KeyError` will be raised.

        Parameters
        ----------
        query: str
            tag to query
        fail_on_missing: bool, default False
            if `True` raise a KeyError is no match is found
            else return `None`

        Returns
        -------
        Optional[TaggedT]
            will be `None` if no matching object is found
            else the matching object

        Raises
        ------
        KeyError
            If `fail_on_missing` is `True` and no matching object is found
        """
        _obj = self._tag_lookup_table.get(query, None)
        if _obj is None and fail_on_missing:
            raise KeyError(f"Tag '{query}' not found in {type(self).__name__}")
        return _obj
