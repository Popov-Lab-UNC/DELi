from functools import wraps
from os import PathLike
from pathlib import Path
from typing import Any, Callable, Literal, Optional, ParamSpec, TypeVar, Union
import inspect

from deli.constants import DELI_DATA_DIR

P = ParamSpec("P")
R = TypeVar("R")

DELI_DATA_DIR = Path(DELI_DATA_DIR) if DELI_DATA_DIR else None

def validate_file_path(
        func: Optional[Callable] = None,
        *,
        sub_dir: Optional[Literal[
            "building_blocks", "libraries", "indexes", "barcodes"
        ]] = None) -> Callable:

    def build_path(file_path: Union[str, PathLike]) -> Path:
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
        file_path = Path(file_path)
        if not file_path.exists():
            if DELI_DATA_DIR is not None:
                if sub_dir is None:
                    file_path = DELI_DATA_DIR/file_path
                else:
                    file_path = DELI_DATA_DIR/sub_dir/file_path
                if not file_path.exists():
                    raise FileNotFoundError(
                        f"cannot locate bb file {file_path.name} "
                        f"in DELI_DATA_DIR"
                    )
            else:
                raise FileNotFoundError(
                    f"cannot find bb file {file_path}; DELI_DATA_DIR not set"
                )
        return file_path

    decorator = build_validation_decorator(build_path, 'file_path')
    if func:
        return decorator(func)
    return decorator


def build_validation_decorator(validator_func: Callable[[Any], Any], target_arg_name: str):
    """
    General decorator generator for validating a given argument with a validator function.
    """
    def outer(func: Callable[P, R]) -> Callable[P, R]:
        @wraps(func)
        def inner(*args: P.args, **kwargs: P.kwargs) -> R:
            pos = None
            if target_arg_name in kwargs:
                target_arg_value = kwargs[target_arg_name]
            else:
                signature = inspect.signature(func)
                pos = list(signature.parameters.keys()).index(target_arg_name)
                target_arg_value = args[pos]
            new_arg = validator_func(target_arg_value)
            if pos:
                args = (*args[:pos], new_arg, *args[pos+1:])
            else:
                kwargs[target_arg_name] = new_arg
            return func(*args, **kwargs)
        return inner
    return outer
