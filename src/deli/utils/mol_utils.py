"""utility functions related to molecules"""

from typing import Literal, Optional, Union, overload

from rdkit import Chem
from rdkit.Chem.MolStandardize.rdMolStandardize import LargestFragmentChooser
from rdkit.rdBase import BlockLogs


_LARGEST_FRAGMENT_CHOOSER = LargestFragmentChooser(preferOrganic=True)

Molable = Union[str, Chem.Mol]


@overload
def to_mol(smi: Molable, fail_on_error: Literal[False]) -> Chem.Mol: ...


@overload
def to_mol(smi: Molable, fail_on_error: Literal[True]) -> Optional[Chem.Mol]: ...


def to_mol(smi: Molable, fail_on_error: bool = True) -> Optional[Chem.Mol]:
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
