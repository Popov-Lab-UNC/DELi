"""deli init file"""

from importlib.metadata import PackageNotFoundError, version


try:
    __version__ = version("deli-chem")
except PackageNotFoundError:
    __version__ = "unknown"
