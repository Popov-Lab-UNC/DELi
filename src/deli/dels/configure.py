from typing import Optional, Literal
import os

from deli.constants import DELI_DATA_DIR


def check_file_path(
        filepath: str,
        sub_dir: Optional[Literal[
            "building_blocks", "libraries", "indexes", "barcodes"
        ]] = None) -> str:
    """
    Check for file existence

    Notes
    -----
    First asks if file exists, if not will ask if a file
    with that name exists in the DELI_DATA_DIR

    Parameters
    ----------
    filepath: str
        name of file if in DELI_DATA_DIR,
        OR path to file
    sub_dir: "building_blocks", "libraries", "indexes", or "barcodes"
        name of sub-directory in DELI_DATA_DIR to search for file in
        if None, look in DELI_DATA_DIR

    Returns
    -------
    str
        absolute path to file
    """
    if not os.path.exists(filepath):
        if DELI_DATA_DIR is not None:
            if sub_dir is None:
                filepath = os.path.join(DELI_DATA_DIR, filepath)
            else:
                filepath = os.path.join(DELI_DATA_DIR, sub_dir, filepath)

            if not os.path.exists(filepath):
                raise FileNotFoundError(
                    f"cannot locate bb file {os.path.basename(filepath)} "
                    f"in DELI_DATA_DIR"
                )
        else:
            raise FileNotFoundError(
                f"cannot find bb file {filepath}; DELI_DATA_DIR not set"
            )
    return filepath
