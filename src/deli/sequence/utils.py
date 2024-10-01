"""Sequence utilities."""


def reverse_compliment(seq: str) -> str:
    """
    Get the reverse compliment of a DNA sequence

    Parameters
    ----------
    seq: str
        the DNA sequence to generate the reverse compliment for

    Returns
    -------
    rev_comp: str
        the reverse compliment of the passed DNA sequence
    """
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement.get(base, base) for base in reversed(seq))
