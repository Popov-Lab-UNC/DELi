"""utility functions related to molecules"""

from typing import Literal, Union, overload

from rdkit import Chem
from rdkit.Chem import Mol
from rdkit.Chem.MolStandardize.rdMolStandardize import LargestFragmentChooser
from rdkit.rdBase import BlockLogs


_LARGEST_FRAGMENT_CHOOSER = LargestFragmentChooser(preferOrganic=True)

Molable = Union[str, Chem.Mol]


class ChemicalObjectError(Exception):
    """raised if there is an issue with chemical objects"""


class SmilesMixin:
    """
    Mixin for Compound objects that are expected to have a fully enumerated SMILES

    This class is used as the main way to check if a compound can be used
    in downstream applications that require a SMILES strings.
    """

    _smiles: None | str = None
    _mol: None | Mol = None

    def has_smiles(self) -> bool:
        """Returns True if the compound has a non-null SMILES string"""
        return self._smiles is not None

    @property
    def smi(self) -> str:
        """The enumerated SMILES of the compound"""
        if self._smiles is not None:
            return self._smiles
        else:
            raise ChemicalObjectError(f"{self.__class__.__name__} object missing SMILES string")

    @smi.setter
    def smi(self, value):
        """Cannot set the SMILES for a compound, once created it is immutable"""
        raise ChemicalObjectError(
            f"Cannot set SMILES for {self.__class__.__name__} object directly; "
            f"SMILES can only be set at initialization"
        )

    @smi.deleter
    def smi(self):
        """Cannot delete the SMILES for a compound, once created it is immutable"""
        raise ChemicalObjectError(f"Cannot delete SMILES for {self.__class__.__name__} object")

    @property
    def mol(self) -> Mol:
        """The RDKit Mol object for the compound; will cache it after first access"""
        if self._mol is None:
            if self._smiles is not None:
                try:
                    _mol = to_mol(self._smiles, fail_on_error=True)
                    self._mol = _mol
                except ValueError as e:
                    raise ChemicalObjectError(
                        f"Cannot create RDKit Mol from SMILES: {self._smiles}"
                    ) from e
            else:
                raise ChemicalObjectError(
                    f"{self.__class__.__name__} object missing SMILES string, "
                    f"cannot generate mol attribute"
                )
        return self._mol

    @mol.setter
    def mol(self, value):
        """No support for setting mol directly; any modification should be done outside DELi"""
        raise ChemicalObjectError(
            f"Cannot change `mol` for {self.__class__.__name__} object directly; "
            f"Derived from `smi` which can only be set at initialization\n"
        )

    @mol.deleter
    def mol(self):
        """Delete the cache of the Mol object (if cached)"""
        self._mol = None


@overload
def to_mol(smi: Molable, fail_on_error: Literal[False]) -> Chem.Mol: ...


@overload
def to_mol(smi: Molable, fail_on_error: Literal[True]) -> Chem.Mol | None: ...


def to_mol(smi: Molable, fail_on_error: bool = True) -> Chem.Mol | None:
    """
    Given an object, attempt to convert it to a Chem.Mol object

    Notes
    -----
    Can only covert valid SMILES str and rdkit.Mol objects.
    If a rdkit.Mol object is passed, the same object will be returned.

    Parameters
    ----------
    smi: Any
        object to convert to a Chem.Mol
    fail_on_error: bool
        whether to raise an exception when converting fails.
        if True, will return None when a conversion fails

    Returns
    -------
    Chem.Mol

    Raises
    ------
    ValueError
        if the SMILES cannot be parsed by rdkit
    TypeError
        if the passed object is not a type that can be converted to a Chem.Mol
    """
    _lock = BlockLogs()  # this turns off the rdkit logger
    if isinstance(smi, Chem.Mol):
        return smi
    elif isinstance(smi, str):
        _mol = Chem.MolFromSmiles(smi)
        if _mol is None:
            if fail_on_error:
                raise ValueError(f"SMILES {smi} cannot be parsed by RDKit")
        return _mol
    else:
        if fail_on_error:
            raise TypeError(f"cannot convert type {type(smi)} to type rdkit.Mol")
        else:
            return None


def to_smi(mol: Chem.Mol) -> str:
    """
    Given a Chem.Mol object, convert it to a SMILES

    Notes
    -----
    It is my understanding that the Chem.MolToSmiles function cannot fail
    for sanitized Mols. This function expects all Mols passed are sanitized

    Parameters
    ----------
    mol: Any
        object to convert to a Chem.Mol

    Returns
    -------
    str
    """
    _lock = BlockLogs()  # this turns off the rdkit logger
    return Chem.MolToSmiles(mol)


def check_valid_smiles(smi: str) -> bool:
    """
    Checks if a SMILES string is valid (readable by RDKit)

    Parameters
    ----------
    smi: str
        SMILES string to check

    Returns
    -------
    bool
    """
    _block = BlockLogs()
    return Chem.MolFromSmiles(smi) is not None


def get_largest_fragment(mol: Chem.Mol) -> Chem.Mol:
    """
    Given a Chem.Mol, return its largest fragment as a new Mol object

    Parameters
    ----------
    mol: Chem.Mol
        Mol to find the largest fragment in

    Returns
    -------
    Chem.Mol
        the largest fragment as a new Mol object
    """
    return _LARGEST_FRAGMENT_CHOOSER.choose(mol)
