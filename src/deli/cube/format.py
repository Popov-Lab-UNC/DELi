import dataclasses
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from deli.decode.decoder import DecodedCompound, DecodedDELCompound, DecodedToolCompound


@dataclasses.dataclass
class CubeFormat:
    """
    Data class representing the format options for cube data export.

    Attributes
    ----------
    include_bb_ids : bool, default=False
        Whether to include building block IDs in the output.
        Will include each cycle as a separate column.
    include_bb_smiles : bool, default=False
        Whether to include building block SMILES in the output.
        Will include each cycle as a separate column.
    include_library_id : bool, default=False
        Whether to include library IDs in the output.
    enumerate_smiles : bool, default=False
        Whether to enumerate SMILES strings in the output.
        Will enumerate the SMILES from building blocks if
        a reaction schema is provided with the library.
        Will override any existing SMILES
    include_raw_counts : bool, default=False
        Whether to include raw counts in the output.
        These are the non-umi corrected counts.
    include_umis : bool, default=False
        Whether to include UMIs in the output.
        Will be ';' separated if multiple UMIs exist for a single
        compound entry.

    """
    include_bb_ids: bool = False
    include_bb_smiles: bool = False
    include_library_id: bool = False
    enumerate_smiles: bool = False
    include_raw_counts: bool = False
    include_umis: bool = False

    @classmethod
    def parse_cube_format_string(cls, format_string: str) -> "CubeFormat":
        """
        Parses a cube format string and returns a dictionary with its components.

        Read the 'Cube File' docs for more information on the format string.

        Parameters
        ----------
        format_string : str
            The cube format string to parse.

        Returns
        -------
        CubeFormat
            An instance of CubeFormat with the parsed options.
        """
        format_string = format_string.strip()

        parsed_options = {
            'include_bb_ids': False,
            'include_bb_smiles': False,
            'include_library_id': False,
            'enumerate_smiles': False,
            'include_raw_counts': False,
            'include_umis': False
        }

        for option in format_string:
            match option:
                case 'b':
                    parsed_options['include_bb_ids'] = True
                case 'B':
                    parsed_options['include_bb_smiles'] = True
                case 'l':
                    parsed_options['include_library_id'] = True
                case 's':
                    parsed_options['enumerate_smiles'] = True
                case 'r':
                    parsed_options['include_raw_counts'] = True
                case 'u':
                    parsed_options['include_umis'] = True
                case _:
                    raise ValueError(f"Unknown cube format option: '{option}'")
        return cls(**parsed_options)


def get_decoded_compound_format_dict(
        format_settings: CubeFormat, compound: "DecodedCompound", skip_compound_id: bool = False
) -> dict[str, str]:
    """
    Convert a DecodedCompound into a dict of fields based on the format options.

    This will only include fields specific to the compound, not its count information

    Parameters
    ----------
    format_settings : CubeFormat
        The format options to use for formatting the compound.
    compound : DecodedCompound
        The compound to format.
    skip_compound_id : bool, default=False
        Whether to skip including the compound ID in the output.

    Returns
    -------
    dict[str, str]
        the formatted compound as a dictionary of fields.
    """
    output_fields = {}

    if format_settings.include_library_id:
        output_fields["lib_id"] = compound.get_library_id()

    if isinstance(compound, DecodedDELCompound):
        if not skip_compound_id:
            output_fields["compound_id"] = compound.compound_id

        if format_settings.include_bb_ids or format_settings.include_bb_smiles:
            for i, bb in enumerate(compound.building_blocks):
                if format_settings.include_bb_ids:
                    output_fields[f"BB{i:02d}_id"] = bb.bb_id
                if format_settings.include_bb_smiles:
                    output_fields[f"BB{i:02d}_smiles"] = bb.smi

    elif isinstance(compound, DecodedToolCompound):
        if not skip_compound_id:
            output_fields["compound_id"] = compound.tool_compound.compound_id

    return output_fields


class CubeFormater:
    """
    Formats DecodedCompounds into cube file format dict based on specified options.

    Parameters
    ----------
    format_options : CubeFormat
        The format options to use for formatting the cube data.
    """

    def __init__(self, format_options: CubeFormat):
        self.format_options = format_options



    def format_compound(self, compound: "DecodedCompound", skip_compound_id: bool = False) -> dict[str, str]:
        pass


