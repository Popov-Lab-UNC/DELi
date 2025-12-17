# def write_cube(
#         self,
#         out_path: str,
#         include_library_id_col: bool = True,
#         include_bb_id_cols: bool = True,
#         include_bb_smi_cols: bool = False,
#         include_raw_count_col: bool = True,
#         enumerate_smiles: bool = False,
#         fail_on_failed_numeration: bool = False,
#         file_format: Literal["csv", "tsv"] = "csv",
# ) -> None:
#     """
#     Write the resulting DEL counts to a human-readable separated value file
#
#     Will write a unique result file for each selection in the experiment.
#
#     Each row will be a unique DEL ID with at least 2 and at most 11 columns.
#     If including BB columns, the max cycle number is the number of columns included,
#     with DELs from libraries with smaller cycle counts set to 'null'.
#     Below are the columns in order of appearance:
#     - DEL_ID
#     - LIBRARY_ID (if `include_library_id_col` is True)
#     - <For [#] of bb cycles in the library>
#         - BB[#]_ID (if `include_bb_id_cols` is True)
#         - BB[#]_SMILES (if `include_bb_smi_cols` is True)
#     - ENUMERATED_SMILES
#     - RAW_COUNT (if `include_raw_count_col` is True)
#     - UMI_CORRECTED_COUNT (same as raw count if not using UMI degen)
#
#     If `enumerate_smiles` is True, the enumerated SMILES will be generated on the
#     fly using the reaction data provided in the library files.
#     If this is missing, the SMILES will be set to 'null'.
#     Also note that on the fly enumeration will have a dramatic impact on runtime.
#     For large DELs, it can be far more efficient to enumerate the SMILES
#     once, save in a database, and then map the enumerated SMILES to the DEL ID
#     in the CSV file with a join.
#
#     Notes
#     -----
#     If running decoding in parallel, a call to `to_csv` should be made after
#     all decoding is complete and the counters have been merged.
#     Merging raw csv files will result is data loss and incorrect counts,
#     especially if UMI degen is used.
#
#     Parameters
#     ----------
#     out_path: str
#         path to save results to
#         will create the directory if it does not exist
#     include_library_id_col: bool, default = False
#         include the library ID in the output
#     include_bb_id_cols: bool, default = False
#         include the building block IDs in the output
#     include_bb_smi_cols: bool, default = False
#         include the building block SMILES in the output
#     include_raw_count_col: bool, default = True
#         include the raw count in the output
#     enumerate_smiles: bool, default = False
#         include the enumerated SMILES in the output
#         Note: this will have a dramatic impact on runtime
#     fail_on_failed_numeration: bool, default = False
#         if `True`, will raise an exception if enumeration fails for any DEL
#         if `False`, will set the SMILES to 'null' if enumeration fails
#         only used if `enumerate_smiles` is `True`
#     file_format: Literal["csv", "tsv"] = "csv"
#         which file format to write to
#         either "csv" or "tsv"
#
#     Warnings
#     --------
#     UserWarning
#         - raised when enumeration is requested but some libraries are
#         missing enumerators
#         - raised when building block smiles are requested but some libraries
#         are missing this info
#     """
#     delimiter: str
#     if file_format == "tsv":
#         delimiter = "\t"
#     elif file_format == "csv":
#         delimiter = ","
#     else:
#         raise ValueError(f"file format '{file_format}' not recognized")
#
#     # check for enumeration
#     if enumerate_smiles:
#         if not self.selection.library_collection.all_libs_can_enumerate():
#             _libraries_missing_enumerators = [
#                 lib.library_id for lib in self.selection.library_collection.libraries if not lib.can_enumerate()
#             ]
#             warnings.warn(
#                 f"Some libraries are missing enumerators: "
#                 f"{_libraries_missing_enumerators}. "
#                 "On the fly enumeration will be skipped for these libraries. "
#                 "Please check the library files to ensure they have the "
#                 "correct reaction data.",
#                 stacklevel=2,
#             )
#
#     # check for building block smiles
#     if include_bb_smi_cols:
#         if not self.selection.library_collection.all_libs_have_building_block_smiles():
#             _libraries_missing_bb_smi = [
#                 lib.library_id
#                 for lib in self.selection.library_collection.libraries
#                 if not lib.building_blocks_have_smi()
#             ]
#             warnings.warn(
#                 f"Some libraries are missing building block smiles: "
#                 f"{_libraries_missing_bb_smi}. "
#                 "Building block smiles will be skipped for these libraries. "
#                 "Please check the building block files to ensure they have "
#                 "the correct smiles.",
#                 stacklevel=2,
#             )
#
#     # get max_cycle_size
#     _max_cycle_size = (
#         self.selection.library_collection.max_cycle_size() if (include_bb_id_cols or include_bb_smi_cols) else 0
#     )
#
#     # set header
#     _header = "DEL_ID"
#     if include_library_id_col:
#         _header += f"{delimiter}LIBRARY_ID"
#     for i in range(_max_cycle_size):
#         if include_bb_id_cols:
#             _header += f"{delimiter}BB{i + 1}_ID"
#         if include_bb_smi_cols:
#             _header += f"{delimiter}BB{i + 1}_SMILES"
#     if enumerate_smiles:
#         _header += f"{delimiter}SMILES"
#     if include_raw_count_col:
#         _header += f"{delimiter}RAW_COUNT"
#     _header += f"{delimiter}UMI_CORRECTED_COUNT\n"
#
#     if self.degen is not None:
#         if os.path.dirname(out_path) != "":
#             os.makedirs(os.path.dirname(out_path), exist_ok=True)
#
#         # write the file
#         with open(out_path, "w") as f:
#             # write header
#             f.write(_header)
#             for library_count_set in self.degen.del_counter.values():
#                 for decoded_compound, counts in library_count_set.items():
#                     _row_dict = decoded_compound.to_cube_row_dict()
#                     _row = _row_dict["compound_id"]
#                     if include_library_id_col:
#                         _row += f"{delimiter}{_row_dict.get('library_id', 'null')}"
#                     for i in range(1, self._max_cycle_size + 1):
#                         if include_bb_id_cols:
#                             _row += f"{delimiter}{_row_dict.get(f'BB{i:02d}_id', 'null')}"
#                         if include_bb_smi_cols:
#                             _row += f"{delimiter}{_row_dict.get(f'BB{i:02d}_smiles', 'null')}"
#
#                     if enumerate_smiles:
#                         try:
#                             _smi = decoded_compound.get_smiles()
#                         except (EnumerationRunError, RuntimeError, ValueError) as e:
#                             if fail_on_failed_numeration:
#                                 raise RuntimeError(
#                                     f"Failed to enumerate SMILES for decoded DEL ID {decoded_compound.compound_id}"
#                                 ) from e
#                             _smi = "null"
#                         _row += f"{delimiter}{_smi}"
#                     if include_raw_count_col:
#                         _row += f"{delimiter}{counts.get_raw_count()}"
#                     _row += f"{delimiter}{counts.get_degen_count()}\n"
#                     f.write(_row)
#     else:
#         raise RuntimeError("Cannot write cube, no degeneration counter found")
