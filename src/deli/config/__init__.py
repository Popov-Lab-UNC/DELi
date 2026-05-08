"""Configuration handling for DELi and is various modules

This will load the config settings from a TOML file and provide a
global config object via the `get_config` function which lazily loads
the config, caching it for subsequent calls.
Config settings are validated upon loading.

The config file is expected to be in TOML format and contain the following sections:
- [deli_data]: settings related to the DELi data directory
- [tokens]: settings related to tokenization of DEL compounds
- [dna]: settings related to DNA encoding/decoding
- [decode]: settings related to the decoding process

The default location for the config file is ~/.deli/deli_config.toml,
can be overridden using the environment variable DELI_CONFIG_PATH.
It can also be overridden by specifying a path when calling `get_config`,
which will take precedence over the environment variable.

The DELi Data Directory is also defined in this module.
It uses fsspec to provide a filesystem-agnostic interface to a directory
containing user-specified DEL data, which can be local or remote (e.g. S3, GCS, etc.).
The default is to use a local filesystem, rooted at ~/.deli/deli_data.
"""
import os
import warnings
from pathlib import Path

from .configure import DELiConfigError, _DeliConfig
from .deli_data_dir import DELI_DATA_SUB_DIRS, DELiDataDirError
from .parsing import _DEFAULT_CONFIG_PATH, _load_config
from .settings import DecodingSettings


_DELI_CONFIG: _DeliConfig | None = None

def get_deli_config(config_path: os.PathLike[str] | None = None, force_reload: bool = False) -> "_DeliConfig":
    """
    Get the DELi config, loading it lazily if not already loaded

    This method uses an LRU cache to avoid reloading the config multiple times.
    If a new config_path is provided that differs from the previously defined
    config, the new config will be loaded.

    Notes
    -----
    The config will *not* update if the file contents are modified after it has
    been loaded. If you need to reload the config after.
    """
    global _DELI_CONFIG

    config_path = config_path or Path(os.getenv("DELI_CONFIG_PATH", _DEFAULT_CONFIG_PATH))

    # lazy loading of the deli config
    if _DELI_CONFIG is None:
        _DELI_CONFIG = _load_config(
            config_path=config_path,
            create_if_missing=config_path == _DEFAULT_CONFIG_PATH
        )
    else:
        if force_reload:
            _DELI_CONFIG = _load_config(config_path=config_path, create_if_missing=False)
        elif config_path and (config_path != _DELI_CONFIG.source):
            new_config = _load_config(config_path=config_path, create_if_missing=False)
            if new_config != _DELI_CONFIG:
                warnings.warn(
                    "A new DELi config was loaded that differs from the previously loaded config. "
                    "This could lead to behavior changes mid-execution.",
                    stacklevel=1
                )
            _DELI_CONFIG = new_config
    return _DELI_CONFIG


__all__ = [
    "get_config",
    "DecodingSettings",
    "DELiConfigError",
    "DELiDataDirError",
    "DELI_DATA_SUB_DIRS"
]
