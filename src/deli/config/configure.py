"""Code for handling loading/validating DELi configs for decoding"""

from dataclasses import dataclass, field
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from .deli_data_dir import DELiDataDir
    from .settings import DecodingSettings


_ID_RESERVED_TOKENS: frozenset[str] = frozenset({",", "."})


class DELiConfigError(Exception):
    """Exception raised for errors in the DELi config."""

    pass


@dataclass(frozen=True, slots=True)
class _DeliConfig:
    """
    Struct to hold info on DELi settings

    Parameters
    ----------
    deli_data_dir: Path | None
        where the DELi data directory is located
    bb_mask: str
        the 3-char long token to replace any masked building blocks
        masked building blocks result from synthon-based analysis
    nuc_2_int: dict[str, int]
        the nucleotide-to-integer mapping
    """

    # token settings
    comp_id_sep: str

    # building_block settings
    bb_mask: str
    bb_id_column: str
    bb_smiles_column: str
    bb_tag_column: str
    bb_subset_id_column: str
    bb_is_null_column: str
    bb_set_subset_id_sep: str

    # generic settings
    nuc_2_int: dict[str, int]

    # deli data directory
    data_directory_filesystem: DELiDataDir

    # decode settings
    decode_settings: DecodingSettings

    # source
    source: str = field(default="memory")

    # auto loaded settings
    id_reserved_tokens: frozenset[str] = field(default_factory=frozenset, init=False, repr=False)

    def __post_init__(self):
        """Validate the config values after initialization."""
        # validate comp_id_sep
        if self.comp_id_sep in _ID_RESERVED_TOKENS:
            raise DELiConfigError(
                f"'comp_id_sep' cannot be any of the reserved tokens: {_ID_RESERVED_TOKENS}"
            )
        if len(self.comp_id_sep) != 1:
            raise DELiConfigError("'comp_id_sep' must be a single character")

        # set the reserved token set for this instance
        object.__setattr__(self, "id_reserved_tokens", frozenset(_ID_RESERVED_TOKENS.union({self.comp_id_sep})))

        # validate nuc_2_int
        if set(self.nuc_2_int.keys()) != {"A", "T", "G", "C"}:
            raise DELiConfigError(
                f"'nuc_2_int' must contain the the nucleotides 'A', 'T', 'G', 'C'; "
                f"found nucleotides '{set(self.nuc_2_int.keys())}'"
            )
        if set(self.nuc_2_int.values()) != {0, 1, 2, 3}:
            raise DELiConfigError(
                f"'nuc_2_int' must map nucleotides to values 0, 1, 2, 3; found values '{set(self.nuc_2_int.values())}'"
            )

        # validate bb_mask
        if any(token in self.bb_mask for token in self.id_reserved_tokens):
            raise DELiConfigError(
                f"'bb_mask' cannot contain any of the reserved tokens: {self.id_reserved_tokens}"
            )

        # check columns are unique for building block files
        bb_columns = [
            self.bb_id_column,
            self.bb_smiles_column,
            self.bb_tag_column,
            self.bb_subset_id_column,
            self.bb_is_null_column,
        ]
        if len(set(bb_columns)) != len(bb_columns):
            raise DELiConfigError("Building block file column names must be unique")

    def __eq__(self, other):
        """Equality check based on config values (ignoring source)"""
        if not isinstance(other, _DeliConfig):
            return NotImplemented
        return (
            self.comp_id_sep == other.comp_id_sep and
            self.bb_mask == other.bb_mask and
            self.nuc_2_int == other.nuc_2_int and
            self.data_directory_filesystem == other.data_directory_filesystem and
            self.decode_settings == other.decode_settings
        )

    def check_id_for_reserved_tokens(self, query_id: str) -> None:
        """
        Check that an ID contains none of this config instance's reserved tokens.

        Parameters
        ----------
        query_id: str
            the ID to check

        Raises
        ------
        ValueError
            If the ID contains any reserved tokens
        """
        for token in self.id_reserved_tokens:
            if token in query_id:
                raise ValueError(
                    f"ID '{query_id}' contains a reserved token/phrase: {self.id_reserved_tokens}"
                )

    def check_if_bb_id_is_reserved(self, bb_id: str) -> None:
        """
        Check that a building block ID is valid (i.e. does not contain reserved tokens)

        This extends `check_id_for_reserved_tokens` by also checking that the building
        block ID is not the same as the `bb_mask`.

        Parameters
        ----------
        bb_id: str
            the building block ID to check

        Raises
        ------
        ValueError
            If the building block ID is invalid
        """
        self.check_id_for_reserved_tokens(bb_id)
        if bb_id == self.bb_mask:
            raise ValueError(
                f"building block id '{bb_id}' is the same as the null mask '{self.bb_mask}'"
            )

    def check_if_bb_set_id_is_reserved(self, bb_set_id: str) -> None:
        """
        Check that a building block set ID is valid

        A building block set ID is considered invalid if it contains the
        bb_set_subset_id_sep

        Parameters
        ----------
        bb_set_id: str
            the building block set ID to check

        Raise
        -----
        ValueError
            If the building block set ID is invalid
        """
        if self.bb_set_subset_id_sep in bb_set_id:
            raise ValueError(
                f"building block set id '{bb_set_id}' contains the subset id separator '{self.bb_set_subset_id_sep}'"
            )
