def check_torch_available():
    """
    Check if PyTorch and torch_geometric are available for import.

    Returns:
        bool: True if both torch and torch_geometric can be imported, False otherwise.
    """
    try:
        import torch
        import torch_geometric

        return True
    except ImportError:
        return False


def create_deli_cube(data, id_col, indexes, control_cols=None, lib_size=None, raw_indexes=None):
    """
    Create a DELi_Cube instance with maximum available functionality.

    This factory function automatically detects available dependencies and creates
    the most capable version of DELi_Cube possible. If PyTorch is available,
    creates DELi_Cube_Full with GNN capabilities. Otherwise, creates DELi_Cube_Core
    with basic functionality.

    Parameters:
        data (pd.DataFrame): The DataFrame containing the data.
        id_col (str): The column name for the ID column.
        indexes (dict): A dictionary mapping experimental IDs to index ranges.
        control_cols (dict, optional): A dictionary mapping experimental IDs to control columns.
        lib_size (int, optional): The size of the library.
        raw_indexes (dict, optional): A dictionary mapping experimental IDs to raw index ranges.

    Returns:
        DELi_Cube_Full: If PyTorch dependencies are available.
        DELi_Cube_Core: If PyTorch dependencies are not available.

    Example:
        >>> cube = create_deli_cube(df, 'DEL_ID', indexes_dict)
        Creating DELi_Cube with full functionality (including GNN)
        >>> cube.gnn_classifier()  # Available if torch is installed
    """
    if check_torch_available():
        from .cube_class_full import DELi_Cube_Full

        return DELi_Cube_Full(
            data=data,
            id_col=id_col,
            indexes=indexes,
            control_cols=control_cols,
            lib_size=lib_size,
            raw_indexes=raw_indexes,
        )
    else:
        print("Creating basic DELi_Cube")
        print("   GNN methods unavailable - requires: pip install torch torch_geometric")
        from .cube_class import DELi_Cube_Core

        return DELi_Cube_Core(
            data=data,
            id_col=id_col,
            indexes=indexes,
            control_cols=control_cols,
            lib_size=lib_size,
            raw_indexes=raw_indexes,
        )
