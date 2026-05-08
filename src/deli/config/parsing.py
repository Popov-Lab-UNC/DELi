"""Parse DELi config files"""

import os
import tomllib
import warnings
from pathlib import Path

from .configure import DELiConfigError, _DeliConfig
from .deli_data_dir import DELiDataDir, _parse_filesystem_config
from .settings import DecodingSettings


_DEFAULT_CONFIG_PATH = Path.home() / ".deli" / "config.toml"

def _write_default_config(config_path: os.PathLike[str]) -> None:
    """Write the default DELi config to a TOML file at the specified path"""
    default_config_string = """
[deli_data]
filesystem = "local"

[deli_data.local]
dir = "~/.deli/deli-data"

[dna]
A = 0
G = 1
T = 2
C = 3

[tokens]
comp_id_sep = "-"

[building_blocks]
bb_mask = "###"
bb_id_column = "id"
bb_smiles_column = "smiles"
bb_tag_column = "tag"
bb_subset_id_column = "subset_id"
bb_is_null_column = "is_null"
bb_set_subset_id_sep = ":::" # separator for building block set ids and subset ids (used in reactions)

[decode]
ignore_tool_compounds = true
demultiplexer_algorithm = "regex" # "cutadapt", "regex", "full"
demultiplexer_mode = "single" # "library", "single", "flanking"
realign = false
library_error_tolerance = 1
library_error_correction_mode_str = "levenshtein_dist:2,asymmetrical"
min_library_overlap = 8
revcomp = true
library_wiggle = false
wiggle = false
decode_matching_approach = "first_perfect" # "greedy", "first_best", "first_perfect", "search_all"
max_read_length = inf
min_read_length = -inf
default_error_correction_mode_str = "levenshtein_dist:1,asymmetrical"
    """
    config_path = Path(config_path)
    config_path.parent.mkdir(parents=True, exist_ok=True)
    with open(config_path, "w") as f:
        f.write(default_config_string)


def _load_config(config_path: os.PathLike[str], create_if_missing: bool = True) -> _DeliConfig:
    """Load the DELi config from a TOML file"""
    _config_path = Path(config_path)

    if not _config_path.exists():
        if create_if_missing:
            warnings.warn(
                f"No DELi config found at default path '{config_path}'; writing a default config to this path.",
                UserWarning, stacklevel=1
            )
            _write_default_config(_config_path)
        else:
            raise FileNotFoundError(f"Config file not found at '{config_path}'")
    else:
        raise FileNotFoundError(f"DELi Config file not found at '{config_path}'")

    try:
        with open(_config_path, "rb") as f:
            config_data = tomllib.load(f)
    except Exception as e:
        raise DELiConfigError(f"Unable to parse DELi config from path '{config_path}'") from e

    # parse the config data into a _DeliConfig object
    filesystem_type = config_data["deli_data"].get("filesystem", "local")
    filesystem_config = config_data["deli_data"].get(filesystem_type, {})
    deli_data_dir = DELiDataDir(_parse_filesystem_config(filesystem_type, filesystem_config))

    decode_settings = DecodingSettings(**config_data.get("decode", {}))

    return _DeliConfig(
        comp_id_sep=config_data["tokens"].get("comp_id_sep", "-"),
        bb_mask=config_data["building_blocks"].get("bb_mask", "###"),
        bb_id_column = config_data["building_blocks"].get("bb_id_column", "id"),
        bb_smiles_column = config_data["building_blocks"].get("bb_smiles_column", "smiles"),
        bb_tag_column = config_data["building_blocks"].get("bb_tag_column", "tag"),
        bb_subset_id_column = config_data["building_blocks"].get("bb_subset_id_column", "subset_id"),
        bb_is_null_column = config_data["building_blocks"].get("bb_is_null_column", "is_null"),
        bb_set_subset_id_sep = config_data["building_blocks"].get("bb_set_subset_id_sep", ":::"),
        nuc_2_int=config_data.get("dna", {"A": 0, "G": 1, "T": 2, "C": 3}),
        data_directory_filesystem=deli_data_dir,
        decode_settings=decode_settings,
        source=str(config_path),
    )
