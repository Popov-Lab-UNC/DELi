"""utility functions related to molecules"""

from typing import Any, Literal, Optional, overload

from rdkit import Chem
from rdkit.Chem.MolStandardize import LargestFragmentChooser
from rdkit.rdBase import BlockLogs


_LARGEST_FRAGMENT_CHOOSER = LargestFragmentChooser(preferOrganic=True)


@overload
def to_mol(smi: str, fail_on_error: Literal[False]) -> Chem.Mol: ...


@overload
def to_mol(smi: str, fail_on_error: Literal[True]) -> Optional[Chem.Mol]: ...


def to_mol(smi: Any, fail_on_error: bool = False) -> Optional[Chem.Mol]:
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
