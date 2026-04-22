import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib_venn import venn2, venn3
from PIL import ImageDraw
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from sklearn.dummy import DummyClassifier, DummyRegressor
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import accuracy_score, confusion_matrix, r2_score
from sklearn.model_selection import KFold
from tqdm import tqdm

from deli.analysis.poly_o import PolyO


class DELi_Cube:
    def __init__(
        self,
        data: pd.DataFrame,
        id_col: str,
        indexes: dict,
        control_cols: dict = None,
        lib_size: int = None,
        raw_indexes: dict = None,
    ):
        """
        Initialize the DELi_Cube object with the provided data, indexes, and control columns.

        Parameters
        ----------
            data (pd.DataFrame): The DataFrame containing the data.
            id_col (str): The column name for the ID column.
            indexes (dict): A dictionary mapping experimental IDs to index ranges.
                Example:
                tri_synth_dict = {
                    'protein_replicates': ['corrected_index3', 'corrected_index4', 'corrected_index6'],
                    'protein_withinhibitor_replicates': ['corrected_index7', 'corrected_index8', 'corrected_index9']
                }
            control_cols (dict, optional): A dictionary mapping experimental IDs to control columns.
                Example:
                control_cols = {
                    'protein_replicates': 'control_col1',
                    'protein_withinhibitor_replicates': 'control_col2'
                }
            lib_size (int, optional): The size of the library.
            raw_indexes (dict, optional): Optional; retained for compatibility. NSC/SD use ``indexes`` (UMI-corrected).
        """
        self.data = data
        self.id_col = id_col
        self.indexes = indexes
        self.control_cols = control_cols
        self.lib_size = lib_size
        self.raw_indexes = raw_indexes
        self._exclude_control_cols_from_spotfire = False

        if not all(col in self.data.columns for col in ["ID_A", "ID_B", "ID_C"]):
            first_row_id = self.data[self.id_col].iloc[0]
            parts_count = len(first_row_id.split("-"))

            if parts_count == 2:
                new_cols = self.data[self.id_col].str.split("-", expand=True)
                new_cols.columns = ["ID_A", "ID_B"]
                self.data = pd.concat(
                    [
                        self.data.iloc[:, : self.data.columns.get_loc(self.id_col) + 1],
                        new_cols,
                        self.data.iloc[:, self.data.columns.get_loc(self.id_col) + 1 :],
                    ],
                    axis=1,
                )
            elif parts_count == 3:
                new_cols = self.data[self.id_col].str.split("-", expand=True)
                new_cols.columns = ["ID_A", "ID_B", "ID_C"]
                self.data = pd.concat(
                    [
                        self.data.iloc[:, : self.data.columns.get_loc(self.id_col) + 1],
                        new_cols,
                        self.data.iloc[:, self.data.columns.get_loc(self.id_col) + 1 :],
                    ],
                    axis=1,
                )
            else:
                raise ValueError("DEL_ID must be split into 2 or 3 parts.")

    def SD_min(self) -> tuple:
        """
        This is an empirical sampling depth minimum meant to be an additional QC check for a given DEL run.
        Functionally it uses the standard that an enrichment minimum of 10 will be reproducible across replicates,
        thus if we evaluate our NSC_max we can determine if the minimum sampling depth is met.
        DOI: https://doi.org/10.1177/2472555218757718

        Returns
        -------
            tuple: (NSC_max_dict, SD_min_dict, sampling_depth_dict) where:
            - NSC_max_dict (dict): Dictionary with the highest NSC value for each experiment.
            - SD_min_dict (dict): Dictionary with the minimum sampling depth required for each experiment.
            - sampling_depth_dict (dict): Dictionary with the total sampling depth for each experiment.
        """
        if self.lib_size is None:
            raise ValueError(
                "Library size must be provided during initialization to run SD_min method."
            )
        NSC_max_dict = {}
        SD_min_dict = {}
        sampling_depth_dict = {}

        for exp_name, columns in self.indexes.items():
            row_sum = self.data[columns].sum(axis=1)
            total_sampling_depth = row_sum.sum() / self.lib_size

            exp_row_sum = self.data[columns].sum(axis=1)
            exp_NSC = exp_row_sum / total_sampling_depth

            exp_NSC_max = round(exp_NSC.max(), 2)
            NSC_max_dict[f"{exp_name}_NSC_max"] = exp_NSC_max

            SD_min_dict[f"{exp_name}_SD_min"] = 10 / exp_NSC_max
            sampling_depth_dict[f"{exp_name}_sampling_depth"] = total_sampling_depth

        return NSC_max_dict, SD_min_dict, sampling_depth_dict

    def NSC_values(self) -> pd.DataFrame:
        """
        Enrichment factor using normalized counts for a given library member normalized to mean sequence reads
        for all library members.
        DOI: http://dx.doi.org/10.1002/anie.201410736

        Returns
        -------
            pd.DataFrame: DataFrame with NSC values for each library member.
        """
        if self.lib_size is None:
            raise ValueError(
                "Library size must be provided during initialization to run NSC_values method."
            )
        for exp_name, columns in self.indexes.items():
            row_sum = self.data[columns].sum(axis=1)
            total_sampling_depth = row_sum.sum() / self.lib_size

            exp_row_sum = self.data[columns].sum(axis=1)
            exp_NSC = exp_row_sum / total_sampling_depth

            self.data[f"{exp_name}_sum"] = exp_row_sum
            self.data[f"{exp_name}_NSC"] = exp_NSC

        return self.data

    def NSC_enrichment_intervals(self, nsc_df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate normalized sequence count intervals assuming a Poisson distribution.

        Parameters
        ----------
            nsc_df (pd.DataFrame): DataFrame returned from NSC_values method.

        Returns
        -------
            pd.DataFrame: DataFrame with NSC enrichment intervals for each library member.
        """
        if self.lib_size is None:
            raise ValueError(
                "Library size must be provided during initialization to run NSC_enrichment_intervals method."
            )

        if not isinstance(nsc_df, pd.DataFrame):
            raise ValueError("Input must be a DataFrame returned from NSC_values method.")

        df = nsc_df.copy()
        for exp_name, columns in self.indexes.items():
            sum_col = f"{exp_name}_sum"
            if sum_col not in df.columns:
                df[sum_col] = df[columns].sum(axis=1)
            mean_seq_count = df[sum_col].sum() / self.lib_size
            df[f"{exp_name}_NSC+"] = ((df[sum_col] + 1) ** 0.5 + 1) ** 2 / mean_seq_count
            df[f"{exp_name}_NSC-"] = ((df[sum_col] + 1) ** 0.5 - 1) ** 2 / mean_seq_count
        self.data = df
        return self.data

    def maximum_likelihood_enrichment_ratio(self) -> pd.DataFrame:
        """
        Calculate the maximum likelihood enrichment ratio to avoid division by zero for members without sequence reads in control.

        DOI: https://doi.org/10.1021/acsomega.3c02152

        Returns
        -------
            pd.DataFrame: DataFrame with maximum likelihood enrichment ratio for each library member.
        """
        if self.control_cols is None:
            raise ValueError(
                "Control columns must be provided during initialization to run maximum_likelihood_enrichment_ratio method."
            )

        df = self.data.copy()

        for exp_name, columns in self.indexes.items():
            control_cols = self.control_cols.get(exp_name)
            if control_cols is None:
                raise ValueError(f"Missing control columns for experiment '{exp_name}'.")

            for control_col in control_cols:
                if control_col not in df.columns:
                    raise ValueError(f"Control column '{control_col}' is missing in data.")

                control_total = df[control_col].sum()
                selection_total = df[columns].sum(axis=1).sum()

                df[f"{exp_name}_MLE"] = (control_total / selection_total) * (
                    ((df[columns].sum(axis=1)) + 3 / 8) / (df[control_col] + 3 / 8)
                )

        self.data = df
        return self.data

    def z_score(self):
        """
        Normalized z-score (zn) under a uniform library null: expected fraction p = 1 / lib_size.

        Uses UMI-corrected counts (``indexes``). DOI: 10.1021/acscombsci.8b00116
        """
        if self.lib_size is None or self.lib_size <= 0:
            raise ValueError(
                "Positive library size (lib_size) is required to run z_score method."
            )

        df = self.data.copy()
        p = 1.0 / self.lib_size

        for exp_name, columns in self.indexes.items():
            df[f"{exp_name}_C"] = df[columns].sum(axis=1)
            selection_total = df[columns].sum(axis=1).sum()

            if selection_total == 0:
                raise ValueError(
                    f"Total selection reads for experiment '{exp_name}' is zero, cannot compute z-score."
                )

            df[f"{exp_name}_E"] = selection_total * p
            sigma = np.sqrt(selection_total * p * (1.0 - p))
            df[f"{exp_name}_sigma"] = sigma

            df[f"{exp_name}_z_score"] = (df[f"{exp_name}_C"] - df[f"{exp_name}_E"]) / sigma
            df[f"{exp_name}_norm_z_score"] = df[f"{exp_name}_z_score"] / np.sqrt(selection_total)

        self.data = df
        return self.data

    def z_score_log_data(self):
        """
        Preprocess the data by calculating the log-transformed z-scores for each experimental group.
        This method performs the following steps:
        1. Creates a copy of the original data.
        2. For each experimental group, calculates the sum and average of the specified columns.
        3. Applies a log transformation to the average values of the experimental group and the control group.
        4. Computes the mean and standard deviation of the log-transformed control group.
        5. Calculates the log-transformed z-scores for the experimental group based on the control group's statistics.

        Returns
        -------
            pd.DataFrame: DataFrame with additional columns for the log-transformed sums, averages, and z-scores for each experimental group.
        """
        if not self.control_cols:
            raise ValueError(
                "Control columns must be provided during initialization to run z_score_log_data method."
            )

        df = self.data.copy()

        for exp_name, columns in self.indexes.items():
            control_cols = self.control_cols.get(exp_name)

            if control_cols is None:
                raise ValueError(f"Missing control columns for experiment '{exp_name}'.")

            df[f"{exp_name}_sum"] = df[columns].sum(axis=1)
            df[f"{exp_name}_avg"] = df[f"{exp_name}_sum"] / len(columns)

            # log transformation on experimental group average
            df[f"{exp_name}_log"] = np.log1p(df[f"{exp_name}_avg"])

            missing_controls = [col for col in control_cols if col not in df.columns]
            if missing_controls:
                raise ValueError(f"Missing control columns in data: {missing_controls}")

            for control_col in control_cols:
                control_values = df[control_col]

                if control_values.sum() == 0:
                    raise ValueError(
                        f"The control column '{control_col}' has a sum of zero. Calculation cannot be done."
                    )

                df[f"{exp_name}_control_log"] = np.log1p(control_values)
                df[f"{exp_name}_control_log"] = df[f"{exp_name}_control_log"].replace(
                    0, np.log1p(control_values[control_values > 0].median())
                )

            # mean and standard deviation of the log-transformed control values
            df[f"{exp_name}_control_log_mean"] = df[f"{exp_name}_control_log"].mean()
            df[f"{exp_name}_control_log_std"] = df[f"{exp_name}_control_log"].std()

            # log-transformed z-score while handling zero standard deviation
            df[f"{exp_name}_z_score_log"] = df.apply(
                lambda row: (
                    (row[f"{exp_name}_log"] - row[f"{exp_name}_control_log_mean"])
                    / row[f"{exp_name}_control_log_std"]
                    if row[f"{exp_name}_control_log_std"] > 0
                    else np.nan
                ),
                axis=1,
            )

            # Handle any missing or problematic z-scores (i.e., divide by zero or invalid calculation)
            if df[f"{exp_name}_z_score_log"].isnull().any():
                raise ValueError(
                    f"Some log-transformed z-scores for experiment '{exp_name}' are invalid (NaN values)."
                )

        self.data = df
        return self.data

    def PolyO(self, feature_mode="disynthon") -> pd.DataFrame:
        """
        Calculate PolyO scores and update the DataFrame with the results.

        Returns
        -------
            pd.DataFrame: The updated DataFrame with PolyO scores.
        """
        if self.lib_size is None or self.lib_size <= 0:
            raise ValueError(
                "Positive library size (lib_size) is required to run PolyO method."
            )
        polyo = PolyO(self.data, self.indexes, self.lib_size, feature_mode=feature_mode)
        polyo.calculate_polyOraw()
        polyo.find_Ccpd()
        polyo.calculate_Cread()
        polyo.poly_o_base()
        polyo.calculate_polyO_score()
        self.data = polyo.data.replace([np.inf, -np.inf], np.nan).fillna(0)
        return self.data

    def normalize(self) -> pd.DataFrame:
        """
        Normalize specified columns by subtracting the control column in the DataFrame,
        ensuring no values go below zero.

        Returns
        -------
            pd.DataFrame: The DataFrame with normalized columns.
        """
        if self.control_cols is None:
            raise ValueError(
                "Control columns must be provided during initialization to run normalize method."
            )

        normalized_df = self.data.copy()

        control_cols_to_drop = set()
        for exp_name, columns in self.indexes.items():
            control_col = self.control_cols.get(exp_name)
            if control_col is None:
                raise ValueError(f"Control column for {exp_name} not provided.")

            if isinstance(control_col, list):
                if not control_col:
                    raise ValueError(f"Control column for {exp_name} is an empty list.")
                control_col = control_col[0]

            if control_col not in normalized_df.columns:
                raise ValueError(f"Control column '{control_col}' not found in the DataFrame")

            print(f"Normalizing {exp_name} with columns: {columns}")

            # Use a temporary DataFrame for normalization
            temp_df = normalized_df[columns].sub(normalized_df[control_col], axis=0)
            temp_df = temp_df.clip(lower=0)
            normalized_df[columns] = temp_df
            normalized_df[f"{exp_name}_avg"] = normalized_df[columns].mean(axis=1)
            # Collect control columns to drop after normalization
            control_cols_to_drop.add(control_col)
            raw_col = control_col.replace("corrected", "raw")
            if raw_col in normalized_df.columns:
                control_cols_to_drop.add(raw_col)

        # Drop all control columns after normalization
        # normalized_df.drop(columns=list(control_cols_to_drop), inplace=True, axis=1)

        self.data = normalized_df
        self._exclude_control_cols_from_spotfire = True
        return self.data

    def disynthonize(self, df: pd.DataFrame = None, synthon_ids: list = None) -> tuple:
        """
        Generate disynthon columns and count data for each experimental group.

        Parameters
        ----------
            df (pd.DataFrame, optional): The DataFrame containing the data. Defaults to self.data.
            synthon_ids (list, optional): List of synthon column names. Defaults to ['ID_A', 'ID_B', 'ID_C'].

        Returns
        -------
            tuple: A tuple containing the updated DataFrame and a dictionary mapping experiment names to disynthon count columns.
        """
        if df is None:
            df = self.data

        if synthon_ids is None:
            synthon_ids = ["ID_A", "ID_B", "ID_C"]

        if len(synthon_ids) not in [2, 3]:
            raise ValueError("synthon_ids must contain exactly two or three elements.")

        missing_cols = [col for col in synthon_ids if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing synthon columns: {', '.join(missing_cols)}")

        new_cols = {}
        new_cols["AB"] = df[synthon_ids[0]] + "-" + df[synthon_ids[1]]
        if len(synthon_ids) == 3:
            new_cols["AC"] = df[synthon_ids[0]] + "-" + df[synthon_ids[2]]
            new_cols["BC"] = df[synthon_ids[1]] + "-" + df[synthon_ids[2]]

        col_order = list(df.columns)
        if "DEL_ID" in col_order:
            del_idx = col_order.index("DEL_ID") + 1
            for col_name, col_values in reversed(new_cols.items()):
                df.insert(del_idx, col_name, col_values)
        else:
            for col_name, col_values in new_cols.items():
                df[col_name] = col_values

        exp_dict = {}

        for exp_id, index_range in self.indexes.items():
            ab_count = df.groupby("AB")[index_range].sum().reset_index()
            ab_count.columns = ["AB"] + [f"count_AB_{exp_id}_{idx}" for idx in index_range]

            if len(synthon_ids) == 3:
                ac_count = df.groupby("AC")[index_range].sum().reset_index()
                ac_count.columns = ["AC"] + [f"count_AC_{exp_id}_{idx}" for idx in index_range]

                bc_count = df.groupby("BC")[index_range].sum().reset_index()
                bc_count.columns = ["BC"] + [f"count_BC_{exp_id}_{idx}" for idx in index_range]

            df = df.merge(ab_count, on="AB", how="left")
            if len(synthon_ids) == 3:
                df = df.merge(ac_count, on="AC", how="left")
                df = df.merge(bc_count, on="BC", how="left")

            for idx in index_range:
                exp_dict[f"AB_{exp_id}"] = exp_dict.get(f"AB_{exp_id}", []) + [
                    f"count_AB_{exp_id}_{idx}"
                ]
                if len(synthon_ids) == 3:
                    exp_dict[f"AC_{exp_id}"] = exp_dict.get(f"AC_{exp_id}", []) + [
                        f"count_AC_{exp_id}_{idx}"
                    ]
                    exp_dict[f"BC_{exp_id}"] = exp_dict.get(f"BC_{exp_id}", []) + [
                        f"count_BC_{exp_id}_{idx}"
                    ]

        for count_cols in exp_dict.values():
            for col in count_cols:
                df[col].fillna(0, inplace=True)

        self.data = df
        return df, exp_dict

    def get_top_disynthons(
        self,
        disynthon_data: pd.DataFrame,
        exp_name1: str,
        comparison_type: str = "control",
        exp_name2: str = None,
        control_name: str = None,
        comparison_metric: str = "avg",
        top_count: int = 20,
        output_dir: str = ".",
    ) -> None:
        """
        Plots a bar chart comparing enrichment between two experimental sets or a single experiment vs control based on a specified metric.
        Creates a separate plot for each synthon type (e.g., AB, BC, AC).

        Parameters
        ----------
            disynthon_data (pd.DataFrame): The DataFrame containing disynthon data.
            exp_name1 (str): The name of the first experiment.
            comparison_type (str): The type of comparison ('control', 'exp2', or 'none').
            exp_name2 (str, optional): The name of the second experiment for comparison_type 'exp2'.
            control_name (str, optional): The name of the control column for comparison_type 'control'.
            comparison_metric (str): The metric for comparison ('sum' or 'avg').
            top_count (int): The number of top disynthons to display.
            output_dir (str): The directory to save the plot.
        """
        if exp_name1 not in self.indexes:
            raise ValueError(f"Experiment '{exp_name1}' not found in the indexes.")

        if comparison_type == "exp2":
            if not exp_name2:
                raise ValueError("exp_name2 must be provided for comparison_type 'exp2'.")
            if exp_name2 not in self.indexes:
                raise ValueError(f"Experiment '{exp_name2}' not found in the indexes.")

        elif comparison_type == "control":
            if not control_name:
                raise ValueError("control_name must be provided for comparison_type 'control'.")
            if control_name not in self.control_cols:
                raise ValueError(f"Control column '{control_name}' not found in the DataFrame.")

        elif comparison_type not in ["control", "exp2", "none"]:
            raise ValueError("comparison_type must be either 'control', 'exp2', or 'none'.")

        if comparison_metric not in ["sum", "avg"]:
            raise ValueError("Comparison metric must be either 'sum' or 'avg'.")

        synthon_ids = []
        if any(col.startswith("AB") for col in disynthon_data.columns):
            synthon_ids.append("AB")
        if any(col.startswith("BC") for col in disynthon_data.columns):
            synthon_ids.append("BC")
        if any(col.startswith("AC") for col in disynthon_data.columns):
            synthon_ids.append("AC")

        if not synthon_ids:
            raise ValueError("No valid synthon columns found in the DataFrame.")

        for synthon in synthon_ids:
            enrichment_results = []

            synthon_cols = [col for col in disynthon_data.columns if col.startswith(synthon)]

            index_cols1 = self.indexes[exp_name1]
            grouped_counts1 = disynthon_data.groupby(synthon)[index_cols1].sum().reset_index()
            grouped_counts1["exp"] = exp_name1

            if comparison_metric == "avg":
                grouped_counts1["metric"] = grouped_counts1[index_cols1].mean(axis=1)
            else:
                grouped_counts1["metric"] = grouped_counts1[index_cols1].sum(axis=1)

            grouped_counts2 = pd.DataFrame()

            if comparison_type == "exp2":
                if exp_name2:
                    index_cols2 = self.indexes[exp_name2]
                    grouped_counts2 = (
                        disynthon_data.groupby(synthon)[index_cols2].sum().reset_index()
                    )
                    grouped_counts2["exp"] = exp_name2

                    if comparison_metric == "avg":
                        grouped_counts2["metric"] = grouped_counts2[index_cols2].mean(axis=1)
                    else:
                        grouped_counts2["metric"] = grouped_counts2[index_cols2].sum(axis=1)

            elif comparison_type == "control":
                if control_name:
                    index_cols2 = self.control_cols[control_name]
                    grouped_counts2 = (
                        disynthon_data.groupby(synthon)[index_cols2].sum().reset_index()
                    )
                    grouped_counts2["exp"] = control_name

                    if comparison_metric == "avg":
                        grouped_counts2["metric"] = grouped_counts2[index_cols2].mean(axis=1)
                    else:
                        grouped_counts2["metric"] = grouped_counts2[index_cols2].sum(axis=1)

            elif comparison_type == "none":
                grouped_counts2 = grouped_counts1

            if not grouped_counts2.empty:
                merged_counts = grouped_counts1.merge(
                    grouped_counts2, on=synthon, suffixes=("_exp1", "_exp2")
                )

                if comparison_type != "none":
                    merged_counts["diff"] = (
                        merged_counts["metric_exp1"] - merged_counts["metric_exp2"]
                    )
                    top_disynthons = merged_counts.nlargest(top_count, "diff")
                else:
                    top_disynthons = merged_counts.nlargest(top_count, "metric_exp1")

                enrichment_results.append(top_disynthons)

            if enrichment_results:
                final_results = pd.concat(enrichment_results)

                fig, ax = plt.subplots(figsize=(12, 6))

                x_labels = final_results[synthon].astype(str)
                x = range(len(x_labels))

                width = 0.35

                if comparison_type == "none":
                    ax.bar(
                        x, final_results["metric_exp1"], width, label=exp_name1, color="#56B4E9"
                    )
                    ax.set_xticks(x)
                    ax.set_xticklabels(x_labels, rotation=90, fontsize=14)
                    ax.set_xlabel(f"{synthon} Disynthon", fontsize=16, labelpad=10)
                    ax.set_ylabel("Enrichment", fontsize=16, labelpad=10)
                    ax.tick_params(axis='y', labelsize=14)
                    ax.set_title(
                        f"Top {top_count} {synthon} Disynthons: {exp_name1}",
                        fontsize=18,
                        fontweight="bold",
                        pad=20
                    )
                    ax.legend(fontsize=14)

                else:
                    x1 = [pos - width / 2 for pos in x]
                    x2 = [pos + width / 2 for pos in x]

                    ax.bar(
                        x1, final_results["metric_exp1"], width, label=exp_name1, color="#56B4E9"
                    )
                    ax.bar(
                        x2,
                        final_results["metric_exp2"],
                        width,
                        label=exp_name2
                        if comparison_type == "exp2"
                        else control_name + "_control",
                        color="#E69F00",
                    )

                    ax.set_xticks(x)
                    ax.set_xticklabels(x_labels, rotation=90, fontsize=14)
                    ax.set_xlabel(f"{synthon} Disynthon", fontsize=16, labelpad=10)
                    ax.set_ylabel("Enrichment", fontsize=16, labelpad=10)
                    ax.tick_params(axis='y', labelsize=14)
                    ax.set_title(
                        f"Top {top_count} {synthon} Disynthons: {exp_name1}{' vs ' + (exp_name2 if comparison_type == 'exp2' else control_name + '_control') if comparison_type != 'none' else ''}",
                        fontsize=18,
                        fontweight="bold",
                        pad=20
                    )
                    ax.legend(fontsize=14)

                plt.tight_layout()

                plot_path = f"{output_dir}/top_{synthon}_disynthons_{exp_name1}{' vs ' + (exp_name2 if comparison_type == 'exp2' else control_name) if comparison_type != 'none' else ''}.svg"
                plt.savefig(plot_path, dpi=300)
                print(f"Plot saved to {plot_path}")

            else:
                print(f"No valid results for {synthon} disynthon.")

    def simple_spotfire_version(self) -> pd.DataFrame:
        """
        Create a Spotfire-friendly version of the DataFrame using normalized data if available
        with corrected indexes, renaming based on experiment IDs and adding an average column for each experiment.
        Ensures that '_sum' columns are present for each experiment and calculates them if missing.

        Returns
        -------
            pd.DataFrame: The Spotfire-friendly DataFrame.
        """
        df = self.data.copy()
        smiles_index = df.columns.get_loc("SMILES")
        cols_to_keep = df.columns[: smiles_index + 1].tolist() + [
            col
            for col in df.columns[smiles_index + 1 :]
            if any(
                pattern in col
                for pattern in [
                    "_sum",
                    "_NSC",
                    "_MLE",
                    "_z_score",
                    "_norm_z_score",
                    "_PolyO_score",
                ]
            )
        ]

        if getattr(self, "_exclude_control_cols_from_spotfire", False) and self.control_cols:
            exclude = set()
            for _exp, ctrl in self.control_cols.items():
                if ctrl is None:
                    continue
                for col in ctrl if isinstance(ctrl, (list, tuple)) else [ctrl]:
                    if not col:
                        continue
                    raw_col = str(col).replace("corrected", "raw")
                    for name in df.columns:
                        if name == col or name == raw_col or name.endswith(f"_{col}") or name.endswith(
                            f"_{raw_col}"
                        ):
                            exclude.add(name)
            cols_to_keep = [c for c in cols_to_keep if c not in exclude]

        for exp_name, index_range in self.indexes.items():
            sum_col_name = f"{exp_name}_sum"
            if sum_col_name not in df.columns:
                df[sum_col_name] = df[index_range].sum(axis=1)
                cols_to_keep.append(sum_col_name)

        # Round relevant columns
        for col in cols_to_keep:
            if any(
                pattern in col
                for pattern in [
                    "_sum",
                    "_NSC",
                    "_MLE",
                    "_z_score",
                    "_norm_z_score",
                    "_PolyO_score",
                ]
            ):
                df[col] = df[col].round(2)

        spotfire_df = df[cols_to_keep]

        for exp_name, index_range in self.indexes.items():
            for column in index_range:
                spotfire_df.rename(columns={column: f"{exp_name}_{column}"}, inplace=True)
            avg_col_name = f"{exp_name}_avg"
            if avg_col_name not in spotfire_df.columns:
                spotfire_df[avg_col_name] = df[index_range].mean(axis=1).round(2)

        return spotfire_df

    def trisynthon_overlap(self, output_dir=None, normalized_data=None, threshold=0.0):
        """
        Create overlap diagrams for tri-synth experiments using normalized data and corrected indexes.

        Parameters
        ----------
            normalized_data (pd.DataFrame, optional): The DataFrame containing normalized data.
               If None, it defaults to using self.data.
            threshold (float): The threshold for filtering counts. Default is 0.0.
        """
        if normalized_data is None:
            normalized_data = self.data

        if output_dir is None:
            output_dir = "."

        valid_experiments = 0  # Track how many experiments have at least two valid indices

        for exp_name, indices in self.indexes.items():
            valid_indices = [idx for idx in indices if idx in normalized_data.columns]
            num_valid_indices = len(valid_indices)

            if num_valid_indices < 2:
                print(
                    f"Skipping experiment '{exp_name}' - not enough valid indices ({num_valid_indices} found)."
                )
                continue
            elif num_valid_indices < len(indices):
                print(
                    f"Warning: Experiment '{exp_name}' has missing indices. Proceeding with {num_valid_indices} valid indices."
                )

            valid_experiments += 1

            filtered_data = normalized_data[
                (normalized_data[valid_indices[0]] > threshold)
                | (normalized_data[valid_indices[1]] > threshold)
            ]

            set_a = set(filtered_data[filtered_data[valid_indices[0]] > threshold]["DEL_ID"])
            set_b = set(filtered_data[filtered_data[valid_indices[1]] > threshold]["DEL_ID"])

            # Use consistent figure size/DPI and margins for uniform scaling
            fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
            if num_valid_indices == 3:
                set_c = set(filtered_data[filtered_data[valid_indices[2]] > threshold]["DEL_ID"])
                v = venn3(
                    [set_a, set_b, set_c],
                    set_labels=(valid_indices[0], valid_indices[1], valid_indices[2]),
                )
                plt.title(f"Overlap Diagram for {exp_name} (Three Indices)", fontsize=18, pad=20)
                for text in v.set_labels:
                    if text: text.set_fontsize(14)
                for text in v.subset_labels:
                    if text: text.set_fontsize(12)
            else:
                v = venn2([set_a, set_b], set_labels=(valid_indices[0], valid_indices[1]))
                plt.title(f"Overlap Diagram for {exp_name} (Two Indices)", fontsize=18, pad=20)
                for text in v.set_labels:
                    if text: text.set_fontsize(14)
                for text in v.subset_labels:
                    if text: text.set_fontsize(12)

            # Enforce equal aspect and fixed margins; avoid tight bbox to keep scale constant
            ax.set_aspect("equal", adjustable="box")
            fig.subplots_adjust(left=0.08, right=0.92, top=0.88, bottom=0.12)
            plt.savefig(f"{output_dir}/{exp_name}_overlap_diagram.svg", dpi=300)
            plt.close()

        if valid_experiments == 0:
            raise ValueError(
                "No experiments had at least two valid indices. Cannot generate overlap diagrams."
            )

    def disynthon_overlap(
        self,
        output_dir: str = None,
        disynthon_data: pd.DataFrame = None,
        disynth_exp_dict: dict = None,
        threshold=0.0,
    ) -> None:
        """
        Create overlap diagrams for disynth experiments using normalized data and corrected indexes.

        Parameters
        ----------
            output_dir (str, optional): The directory to save the overlap diagrams. Defaults to the current directory.
            disynthon_data (pd.DataFrame, optional): The DataFrame containing disynthon data. Defaults to None.
            disynth_exp_dict (dict, optional): A dictionary mapping experiment names to the corresponding corrected index columns. Defaults to None.
            threshold (float, optional): The threshold for filtering counts. Default is 0.0.

        Example of disynth_exp_dict:
            disynth_exp_dict = {
                'AB_NSP2': ['count_AB_corrected_index3', 'count_AB_corrected_index4', 'count_AB_corrected_index6'],
                'AB_NSP2_2034': ['count_AB_corrected_index7', 'count_AB_corrected_index8', 'count_AB_corrected_index9'],
            }
        """
        if disynthon_data is None or disynth_exp_dict is None:
            raise ValueError("disynthon_data and disynth_exp_dict must be provided.")

        disynthon_data = disynthon_data.copy()

        if output_dir is None:
            output_dir = "."

        valid_experiments = 0  # Track number of valid experiments
        overlap_stats = []

        for exp_name, indices in disynth_exp_dict.items():
            if isinstance(threshold, dict) and "mode" in threshold:
                mode = str(threshold.get("mode", "manual")).lower()
                if mode == "auto":
                    percentile = float(threshold.get("percentile", 90))
                    vals = disynthon_data[indices].to_numpy().ravel()
                    vals = vals[np.isfinite(vals)]
                    vals = vals[vals > 0]
                    exp_threshold = float(np.percentile(vals, percentile)) if len(vals) else 0.0
                else:
                    exp_threshold = float(threshold.get("value", 0.0))
            elif isinstance(threshold, dict):
                exp_threshold = float(threshold.get(exp_name, 0.0))
            elif isinstance(threshold, str):
                if threshold.lower() == "auto":
                    vals = disynthon_data[indices].to_numpy().ravel()
                    vals = vals[np.isfinite(vals)]
                    vals = vals[vals > 0]
                    exp_threshold = float(np.percentile(vals, 90)) if len(vals) else 0.0
                else:
                    exp_threshold = float(threshold)
            else:
                exp_threshold = float(threshold)
            # Filter valid indices that exist in the dataframe
            valid_indices = [idx for idx in indices if idx in disynthon_data.columns]
            num_valid_indices = len(valid_indices)

            if num_valid_indices < 2:
                print(
                    f"Skipping experiment '{exp_name}' - not enough valid indices ({num_valid_indices} found)."
                )
                continue
            elif num_valid_indices < len(indices):
                print(
                    f"Warning: Experiment '{exp_name}' has missing indices. Proceeding with {num_valid_indices} valid indices."
                )

            valid_experiments += 1

            mask = disynthon_data[valid_indices].gt(exp_threshold).any(axis=1)
            filtered_data = disynthon_data[mask]
            corr_df = filtered_data[valid_indices].corr()
            corr_vals = corr_df.values[np.triu_indices_from(corr_df.values, k=1)]
            mean_autocorr = float(np.nanmean(corr_vals)) if len(corr_vals) else np.nan
            overlap_stats.append(
                {
                    "experiment": exp_name,
                    "threshold": exp_threshold,
                    "mean_pairwise_autocorr": mean_autocorr,
                    "n_rows": int(len(filtered_data)),
                }
            )

            disynthon_type = exp_name.split("_")[0]
            if disynthon_type not in filtered_data.columns:
                print(
                    f"Skipping experiment '{exp_name}' - column '{disynthon_type}' not found in the DataFrame."
                )
                continue

            set_a = set(
                filtered_data[filtered_data[valid_indices[0]] > exp_threshold][disynthon_type]
            )
            set_b = set(
                filtered_data[filtered_data[valid_indices[1]] > exp_threshold][disynthon_type]
            )

            # Use consistent figure size/DPI and margins across disynthon plots
            fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
            if num_valid_indices == 3:
                set_c = set(
                    filtered_data[filtered_data[valid_indices[2]] > exp_threshold][disynthon_type]
                )
                v = venn3(
                    [set_a, set_b, set_c],
                    set_labels=(valid_indices[0], valid_indices[1], valid_indices[2]),
                )
                plt.title(
                    f"Venn Diagram for {exp_name} (Three Indices) \n Threshold = {exp_threshold}, Mean autocorr = {mean_autocorr:.3f}",
                    fontsize=18,
                    pad=20
                )
                for text in v.set_labels:
                    if text: text.set_fontsize(14)
                for text in v.subset_labels:
                    if text: text.set_fontsize(12)
            else:
                v = venn2([set_a, set_b], set_labels=(valid_indices[0], valid_indices[1]))
                plt.title(
                    f"Venn Diagram for {exp_name} (Two Indices) \n Threshold = {exp_threshold}, Mean autocorr = {mean_autocorr:.3f}",
                    fontsize=18,
                    pad=20
                )
                for text in v.set_labels:
                    if text: text.set_fontsize(14)
                for text in v.subset_labels:
                    if text: text.set_fontsize(12)
            
            # Center the Venn diagram
            plt.subplots_adjust(left=0.1, right=0.9, top=0.85, bottom=0.1)
            
            # Enforce equal aspect and fixed margins; avoid tight bbox to keep scale constant
            ax.set_aspect('equal', adjustable='box')
            plt.savefig(f"{output_dir}/{exp_name}_venn_diagram.svg", dpi=300)
            plt.close()

        if valid_experiments == 0:
            raise ValueError(
                "No experiments had at least two valid indices. Cannot generate overlap diagrams."
            )
        pd.DataFrame(overlap_stats).to_csv(
            os.path.join(output_dir, "disynthon_overlap_stats.csv"), index=False
        )

    def ml_fingerprints_to_RF(self, output_dir: str = None) -> None:
        """
        Trains a Random Forest regressor and a Dummy regressor on trisynthon SMILES fingerprints
        to predict the average normalized enrichment using a subset of 100 molecules.

        The method performs the following steps:
        1. Converts SMILES strings to Morgan fingerprints.
        2. Selects the top 50 molecules with the highest average enrichment and 50 random molecules.
        3. Trains the models using 5-fold cross-validation.
        4. Evaluates the models using R² scores.
        5. Plots the true vs. predicted average enrichment for the last fold.

        Parameters
        ----------
            output_dir (str, optional): Directory to save the scatter plot. Defaults to None.

        Returns
        -------
            None
        """

        def smiles_to_fingerprint(smiles, radius=3, n_bits=2048):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                try:
                    Chem.SanitizeMol(mol)
                    # Use modern MorganGenerator instead of deprecated functions
                    from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

                    fpgen = GetMorganGenerator(radius=radius, fpSize=n_bits)
                    return fpgen.GetFingerprint(mol)
                except:
                    return None
            else:
                return None

        for exp_name, indices in self.indexes.items():
            self.data[f"{exp_name}_average_enrichment"] = self.data[indices].mean(axis=1)
            top_50 = self.data.nlargest(50, f"{exp_name}_average_enrichment")
            remaining_data = self.data.drop(top_50.index)
            random_50 = remaining_data.sample(n=50, random_state=42)
            subset_data = pd.concat([top_50, random_50]).copy()

            subset_data["fingerprints"] = [
                smiles_to_fingerprint(smiles) for smiles in tqdm(subset_data["SMILES"])
            ]
            subset_data.dropna(subset=["fingerprints"], inplace=True)
            X = np.array([list(fp) for fp in subset_data["fingerprints"]])
            y = subset_data[f"{exp_name}_average_enrichment"]

            rf_model = RandomForestRegressor(random_state=42)
            dummy_model = DummyRegressor(strategy="mean")

            kf = KFold(n_splits=5, shuffle=True, random_state=42)
            r2_rf_scores, r2_dummy_scores = [], []

            last_y_test, last_y_pred_rf, last_y_pred_dummy = None, None, None

            for train_idx, test_idx in kf.split(X):
                X_train, X_test = X[train_idx], X[test_idx]
                y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

                rf_model.fit(X_train, y_train)
                dummy_model.fit(X_train, y_train)

                y_pred_rf = rf_model.predict(X_test)
                y_pred_dummy = dummy_model.predict(X_test)

                r2_rf_scores.append(r2_score(y_test, y_pred_rf))
                r2_dummy_scores.append(r2_score(y_test, y_pred_dummy))

                last_y_test, last_y_pred_rf, last_y_pred_dummy = y_test, y_pred_rf, y_pred_dummy

            avg_r2_rf = np.mean(r2_rf_scores)
            avg_r2_dummy = np.mean(r2_dummy_scores)

            plt.figure(figsize=(8, 6))
            plt.scatter(
                last_y_test,
                last_y_pred_rf,
                label=f"Random Forest (Last Fold R² = {r2_rf_scores[-1]:.2f})",
                alpha=0.7,
                s=50,
                color="#0072B2"
            )
            plt.scatter(
                last_y_test,
                last_y_pred_dummy,
                label=f"Dummy Regressor (Last Fold R² = {r2_dummy_scores[-1]:.2f})",
                alpha=0.7,
                marker="x",
                s=50,
                color="#D55E00"
            )
            plt.plot(
                [min(last_y_test), max(last_y_test)],
                [min(last_y_test), max(last_y_test)],
                color="black",
                linestyle="--",
                linewidth=2
            )
            plt.title(
                f"{exp_name} RF Regression: K-Fold Avg R² = {avg_r2_rf:.2f}, Dummy = {avg_r2_dummy:.2f}",
                fontsize=16,
                pad=20
            )
            plt.xlabel("True Average Enrichment", fontsize=14, labelpad=10)
            plt.ylabel("Predicted Average Enrichment", fontsize=14, labelpad=10)
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)
            plt.legend(loc="upper left", fontsize=12)

            if output_dir is None:
                output_dir = "."
            plt.savefig(f"{output_dir}/{exp_name}_RF_regression.svg", dpi=300)
            plt.close()

    def ml_fingerprints_to_classifier(self, threshold: int = 2, output_dir: str = None) -> None:
        """
        Trains a Random Forest classifier and a Dummy classifier on trisynthon SMILES fingerprints
        to predict the enrichment status using a subset of 100 molecules.

        The method performs the following steps:
        1. Converts SMILES strings to Morgan fingerprints.
        2. Selects the top 50 molecules with the highest average enrichment and 50 random molecules.
        3. Trains the models using 5-fold cross-validation.
        4. Evaluates the models using accuracy scores and confusion matrices.
        5. Plots the confusion matrices for the last fold.

        Parameters
        ----------
            threshold (int, optional): Threshold value for enrichment status. Defaults to 2.
            output_dir (str, optional): Directory to save the classifier image. Defaults to None.

        Returns
        -------
            None
        """

        def smiles_to_fingerprint(smiles, radius=3, n_bits=2048):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                try:
                    Chem.SanitizeMol(mol)
                    # Use modern MorganGenerator instead of deprecated functions
                    from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

                    fpgen = GetMorganGenerator(radius=radius, fpSize=n_bits)
                    return fpgen.GetFingerprint(mol)
                except:
                    return None
            else:
                return None

        for exp_name, indices in self.indexes.items():
            self.data[f"{exp_name}_average_enrichment"] = self.data[indices].mean(axis=1)
            top_50 = self.data.nlargest(100, f"{exp_name}_average_enrichment")
            remaining_data = self.data.drop(top_50.index)
            random_50 = remaining_data.sample(n=100, random_state=42)
            subset_data = pd.concat([top_50, random_50]).copy()
            subset_data["fingerprints"] = [
                smiles_to_fingerprint(smiles) for smiles in tqdm(subset_data["SMILES"])
            ]
            subset_data.dropna(subset=["fingerprints"], inplace=True)
            subset_data["target"] = (
                subset_data[f"{exp_name}_average_enrichment"] > threshold
            ).astype(int)
            X = np.array([list(fp) for fp in subset_data["fingerprints"]])
            y = subset_data["target"]
            rf_model = RandomForestClassifier(random_state=42)
            dummy_model = DummyClassifier(strategy="most_frequent")
            kf = KFold(n_splits=5, shuffle=True, random_state=42)
            acc_rf_scores, acc_dummy_scores = [], []
            last_y_test, last_y_pred_rf, last_y_pred_dummy = None, None, None

            for train_idx, test_idx in kf.split(X):
                X_train, X_test = X[train_idx], X[test_idx]
                y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
                rf_model.fit(X_train, y_train)
                dummy_model.fit(X_train, y_train)
                y_pred_rf = rf_model.predict(X_test)
                y_pred_dummy = dummy_model.predict(X_test)
                acc_rf_scores.append(accuracy_score(y_test, y_pred_rf))
                acc_dummy_scores.append(accuracy_score(y_test, y_pred_dummy))
                last_y_test, last_y_pred_rf, last_y_pred_dummy = y_test, y_pred_rf, y_pred_dummy

            avg_acc_rf = np.mean(acc_rf_scores)
            avg_acc_dummy = np.mean(acc_dummy_scores)
            cm_rf = confusion_matrix(last_y_test, last_y_pred_rf)
            cm_dummy = confusion_matrix(last_y_test, last_y_pred_dummy)

            plt.figure(figsize=(12, 5))
            plt.subplot(1, 2, 1)
            sns.heatmap(
                cm_rf,
                annot=True,
                fmt="d",
                cmap="Blues",
                annot_kws={"size": 14},
                xticklabels=["Not Enriched", "Enriched"],
                yticklabels=["Not Enriched", "Enriched"],
            )
            plt.title(
                f"{exp_name} Random Forest Classifier (Avg Acc = {avg_acc_rf:.2f}), Threshold = {threshold} \n"
                f"Last Fold Confusion Matrix Accuracy= {acc_rf_scores[-1]:.2f}",
                fontsize=12,
                pad=15
            )
            plt.xlabel("Predicted", fontsize=14, labelpad=10)
            plt.ylabel("True", fontsize=14, labelpad=10)
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)

            plt.subplot(1, 2, 2)
            sns.heatmap(
                cm_dummy,
                annot=True,
                fmt="d",
                cmap="Blues",
                annot_kws={"size": 14},
                xticklabels=["Not Enriched", "Enriched"],
                yticklabels=["Not Enriched", "Enriched"],
            )
            plt.title(
                f"{exp_name} Dummy Classifier (Avg Acc = {avg_acc_dummy:.2f}), Threshold = {threshold} \n"
                f"Last Fold Confusion Matrix Accuracy= {acc_dummy_scores[-1]:.2f}",
                fontsize=12,
                pad=15
            )
            plt.xlabel("Predicted", fontsize=14, labelpad=10)
            plt.ylabel("True", fontsize=14, labelpad=10)
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)
            plt.tight_layout()

            if output_dir is None:
                output_dir = "."
            plt.savefig(f"{output_dir}/{exp_name}_classifier.svg", dpi=300)
            plt.close()

    def gnn_classifier(
        self,
        threshold: int = 2,
        output_dir: str = ".",
        arch: str = "GAT",
        num_layers: int = 3,
        encoding: str = "embedding",
    ) -> None:
        """
        Train a Graph Neural Network (GNN) classifier to predict the enrichment status of compounds.

        Parameters
        ----------
            threshold (int, optional): Threshold value for enrichment status, by default 2.
            output_dir (str, optional): Directory to save the classifier image, by default ".".
            arch (str, optional): GNN architecture (GAT or GCN), by default "GAT".
            num_layers (int, optional): Number of GNN layers, by default 3.
            encoding (str, optional): Node encoding method (embedding or onehot), by default "embedding".
        """
        try:
            from deli.analysis.analyzers.gnn_analyzer import GNNAnalyzer

            analyzer = GNNAnalyzer(self.data, self.indexes)
            analyzer.classify(
                threshold=threshold,
                output_dir=output_dir,
                arch=arch,
                num_layers=num_layers,
                encoding=encoding,
            )
        except ImportError as e:
            raise ImportError("GNN classifier requires PyTorch") from e

    def top_n_compounds(self, n: int = 20, metric: str = "sum", output_dir: str = None) -> None:
        """
        Display the top N compounds based on a specified metric and save the image to the specified directory.

        Parameters
        ----------
            n (int, optional): Number of top compounds to display. Defaults to 20.
            metric (str, optional): Metric to rank the compounds. Defaults to 'sum'.
            output_dir (str, optional): Directory to save the image. Defaults to None (current directory).
        """
        if metric == "avg":
            for exp_name, index_range in self.indexes.items():
                avg_col_name = f"{exp_name}_avg"
                if avg_col_name not in self.data.columns:
                    self.data[avg_col_name] = self.data[index_range].mean(axis=1).round(2)

        for exp_name, indices in self.indexes.items():
            # Ensure the compound is present at least once in each replicate
            filtered_data = self.data[(self.data[indices] > 0).all(axis=1)]
            top_n = filtered_data.nlargest(n, f"{exp_name}_{metric}")
            mols = [Chem.MolFromSmiles(smiles) for smiles in top_n["SMILES"]]
            legends = [
                f"{row['DEL_ID']}\n\n{metric}: {row[f'{exp_name}_{metric}']:.2f}"
                for _, row in top_n.iterrows()
            ]

            dopts = rdMolDraw2D.MolDrawOptions()
            dopts.legendFontSize = 24
            dopts.bondLineWidth = 1.5
            dopts.minFontSize = 12
            #dopts.annotationFontScale = 0.35  # Smaller atom labels for better readability
            dopts.addAtomIndices = False

            try:
                svg_string = Draw.MolsToGridImage(
                    mols, molsPerRow=4, subImgSize=(300, 300), legends=legends, drawOptions=dopts, useSVG=True
                )
            except AttributeError:
                pass

            if output_dir is None:
                output_dir = "."
            file_name = f"{exp_name}_top_{n}_compounds.svg"

            try:
                with open(f"{output_dir}/{file_name}", "w") as f:
                    f.write(svg_string)
            except (AttributeError, IOError):
                pass

            # Optional notebook display - lazy import to avoid hard dependency on IPython/Jupyter
            try:
                from IPython.display import display, SVG

                print(f"Top {n} compounds for {exp_name} based on {metric}")
                display(SVG(svg_string))
            except ImportError:
                pass

    def top_delta_compounds(self, comparisons: list[dict], output_dir: str = None) -> None:
        """
        Plot top compounds by delta metric between experiment and control.
        """
        if output_dir is None:
            output_dir = "."
        os.makedirs(output_dir, exist_ok=True)

        for spec in comparisons:
            exp_name = spec.get("exp_name")
            control_name = spec.get("control_name")
            metric = spec.get("metric", "avg")
            top_count = int(spec.get("top_count", 20))

            if exp_name not in self.indexes:
                raise ValueError(f"Experiment '{exp_name}' not found in indexes.")
            if control_name not in self.control_cols:
                raise ValueError(f"Control mapping '{control_name}' not found in control_cols.")

            if metric in ["sum", "avg"]:
                exp_vals = self.data[self.indexes[exp_name]].sum(axis=1)
                if metric == "avg":
                    exp_vals = exp_vals / max(len(self.indexes[exp_name]), 1)
            else:
                exp_metric_col = f"{exp_name}_{metric}"
                if exp_metric_col not in self.data.columns:
                    raise ValueError(
                        f"Metric column '{exp_metric_col}' not found. Run that analysis step first."
                    )
                exp_vals = self.data[exp_metric_col]

            control_cols = self.control_cols[control_name]
            control_vals = self.data[control_cols].sum(axis=1)
            if metric == "avg":
                control_vals = control_vals / max(len(control_cols), 1)

            delta_col = f"{exp_name}_minus_{control_name}_{metric}_delta"
            plot_df = self.data[[self.id_col, "SMILES"]].copy()
            plot_df[delta_col] = exp_vals - control_vals
            top_n = plot_df.nlargest(top_count, delta_col)

            mols = [Chem.MolFromSmiles(smiles) for smiles in top_n["SMILES"]]
            legends = [
                f"{row[self.id_col]}\n\nDelta {metric}: {row[delta_col]:.2f}"
                for _, row in top_n.iterrows()
            ]

            dopts = rdMolDraw2D.MolDrawOptions()
            dopts.legendFontSize = 24
            dopts.bondLineWidth = 1.5
            dopts.minFontSize = 12
            dopts.addAtomIndices = False

            svg_string = Draw.MolsToGridImage(
                mols,
                molsPerRow=4,
                subImgSize=(300, 300),
                legends=legends,
                drawOptions=dopts,
                useSVG=True,
            )
            file_name = f"{exp_name}_vs_{control_name}_top_{top_count}_{metric}_deltas.svg"
            with open(f"{output_dir}/{file_name}", "w") as f:
                f.write(svg_string)

    def positive_control_finder(self, positive_control_ID=None):
        """_summary_

        Parameters
        ----------
        positive_control_ID : _type_, optional
            _description_, by default None
        """
        # if positive_control_ID is None:
        #     raise ValueError("Positive control ID must be provided.")

        # if positive_control_ID not in self.data['DEL_ID'].values:
        #     raise ValueError("Positive control ID not found in the data.")

        # if "A***-B***-C***" in self.data['DEL_ID'].values:
        #     print(f'Searching for Trisynthon with ID: {positive_control_ID}')
        #     positive_control_trisynthon = self.data[self.data['DEL_ID'] == positive_control_ID]
        #     #print enrichment value across all indexes for all exp groups
        #     print(positive_control_trisynthon[self.indexes.keys()])
        #     #show structure of positive control
        #     mol = Chem.MolFromSmiles(positive_control_trisynthon['SMILES'].values[0])
        #     img = Chem.Draw.MolToImage(mol)
        #     DrawMoltoGrid([mol], molsPerRow=1, subImgSize=(200, 200), legends=[positive_control_ID])
        # elif "A***-B***" in self.data['DEL_ID'].values:
        #     print(f'Searching for AB Disynthon with ID: {positive_control_ID}')
        # elif "B***-C***" in self.data['DEL_ID'].values:
        #             print(f'Searching for BC Disynthon with ID: {positive_control_ID}')
        # elif "A***" in self.data['DEL_ID'].values:
        #     print(f'Searching for Monosynthon with ID: {positive_control_ID}')
        # else:
        #     print(f'ID incorrect format: {positive_control_ID}')

        # A048-B080-C096
        pass

    def top_trisynthons_diff_from_competitor(self):
        pass

    def top_disynthons_diff_from_competitor(self):
        pass

    def monosynthon_chemical_space(self, output_dir: str = None, experiments: list[str] = None) -> None:
        """
        Create t-SNE chemical space visualizations for monosynthons (A, B, C) colored by average enrichment.

        Parameters
        ----------
            output_dir (str, optional): Directory to save the HTML files. Defaults to None.
        """
        try:
            import plotly.express as px
            import plotly.graph_objects as go
            from sklearn.manifold import TSNE
            from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
        except ImportError as e:
            raise ImportError("monosynthon_chemical_space requires plotly and scikit-learn") from e

        if output_dir is None:
            output_dir = "."

        # Create clusters directory
        clusters_dir = os.path.join(output_dir, "clusters")
        if not os.path.exists(clusters_dir):
            os.makedirs(clusters_dir)

        def smiles_to_fingerprint(smiles, radius=3, n_bits=2048):
            """Convert SMILES to Morgan fingerprint."""
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                try:
                    Chem.SanitizeMol(mol)
                    fpgen = GetMorganGenerator(radius=radius, fpSize=n_bits)
                    return fpgen.GetFingerprint(mol)
                except:
                    return None
            else:
                return None

        experiment_list = experiments if experiments else list(self.indexes.keys())
        # Process each monosynthon type (A, B, C) per experiment
        for synthon_type in ["A", "B", "C"]:
            synthon_col = f"ID_{synthon_type}"

            if synthon_col not in self.data.columns:
                print(f"Skipping {synthon_type} - column {synthon_col} not found")
                continue

            # Get unique monosynthons for this type
            unique_synthons = self.data[synthon_col].unique()

            if len(unique_synthons) < 2:
                print(
                    f"Skipping {synthon_type} - not enough unique monosynthons ({len(unique_synthons)} found)"
                )
                continue

            for exp_name in experiment_list:
                if exp_name not in self.indexes:
                    continue
                synthon_data = []
                for synthon in unique_synthons:
                    synthon_compounds = self.data[self.data[synthon_col] == synthon]
                    if len(synthon_compounds) == 0:
                        continue
                    exp_sum = synthon_compounds[self.indexes[exp_name]].sum(axis=1)
                    avg_enrichment = exp_sum.mean()
                    representative_smiles = synthon_compounds["SMILES"].iloc[0]
                    synthon_data.append(
                        {
                            "synthon": synthon,
                            "smiles": representative_smiles,
                            "avg_enrichment": avg_enrichment,
                            "count": len(synthon_compounds),
                            "exp_name": exp_name,
                        }
                    )
                if len(synthon_data) < 2:
                    continue
                df_synthons = pd.DataFrame(synthon_data)
                print(
                    f"Generating fingerprints for {len(df_synthons)} {synthon_type} monosynthons in {exp_name}..."
                )
                df_synthons["fingerprint"] = [
                    smiles_to_fingerprint(smiles) for smiles in tqdm(df_synthons["smiles"])
                ]
                df_synthons = df_synthons.dropna(subset=["fingerprint"])
                if len(df_synthons) < 2:
                    continue
                X = np.array([list(fp) for fp in df_synthons["fingerprint"]])
                print(f"Running t-SNE for {synthon_type} monosynthons in {exp_name}...")
                tsne = TSNE(
                    n_components=2, random_state=42, perplexity=min(30, len(df_synthons) - 1)
                )
                tsne_results = tsne.fit_transform(X)

                fig = px.scatter(
                    x=tsne_results[:, 0],
                    y=tsne_results[:, 1],
                    color=df_synthons["avg_enrichment"],
                    hover_data={
                        "synthon": df_synthons["synthon"],
                        "avg_enrichment": df_synthons["avg_enrichment"],
                        "count": df_synthons["count"],
                        "exp_name": df_synthons["exp_name"],
                    },
                    color_continuous_scale="viridis",
                    title=f"{synthon_type} Monosynthon Chemical Space ({exp_name})",
                    labels={"x": "t-SNE 1", "y": "t-SNE 2", "color": "Avg Enrichment"},
                )

                # Update hover template
                fig.update_traces(
                    hovertemplate="<b>%{customdata[0]}</b><br>"
                    + "Avg Enrichment: %{customdata[1]:.2f}<br>"
                    + "Count: %{customdata[2]}<br>"
                    + "<extra></extra>"
                )

                # Update layout
                fig.update_layout(
                    width=800,
                    height=600,
                    showlegend=False,
                    plot_bgcolor="white",
                    paper_bgcolor="white",
                    xaxis=dict(showticklabels=False, showgrid=False),
                    yaxis=dict(showticklabels=False, showgrid=False),
                )

                # Save HTML file
                output_file = os.path.join(
                    clusters_dir, f"monosynthon_{synthon_type}_{exp_name}_chemical_space.html"
                )
                fig.write_html(output_file)
                print(f"Saved {synthon_type} monosynthon chemical space to {output_file}")

                # Create Matplotlib figure for SVG
                plt.figure(figsize=(10, 8))
                sc = plt.scatter(
                    tsne_results[:, 0],
                    tsne_results[:, 1],
                    c=df_synthons["avg_enrichment"],
                    cmap="viridis",
                    alpha=0.7,
                )
                plt.colorbar(sc, label="Avg Enrichment")
                plt.title(f"{synthon_type} Monosynthon Chemical Space ({exp_name})")
                plt.xlabel("t-SNE 1")
                plt.ylabel("t-SNE 2")
                plt.xticks([])
                plt.yticks([])
                plt.tight_layout()

                output_file_svg = os.path.join(
                    clusters_dir, f"monosynthon_{synthon_type}_{exp_name}_chemical_space.svg"
                )
                plt.savefig(output_file_svg, dpi=300)
                plt.close()
                print(f"Saved {synthon_type} monosynthon chemical space SVG to {output_file_svg}")

        print(f"Monosynthon chemical space visualizations saved to {clusters_dir}")

    def trisynthon_comparison_space(
        self, comparisons: list[dict], output_dir: str = None, max_points: int = 3000
    ) -> None:
        try:
            import plotly.express as px
            from sklearn.manifold import TSNE
            from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
        except ImportError as e:
            raise ImportError("trisynthon_comparison_space requires plotly and scikit-learn") from e

        if output_dir is None:
            output_dir = "."
        clusters_dir = os.path.join(output_dir, "clusters")
        os.makedirs(clusters_dir, exist_ok=True)

        def smiles_to_fingerprint(smiles, radius=3, n_bits=2048):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                try:
                    Chem.SanitizeMol(mol)
                    fpgen = GetMorganGenerator(radius=radius, fpSize=n_bits)
                    return fpgen.GetFingerprint(mol)
                except Exception:
                    return None
            return None

        for comp in comparisons:
            exp1 = comp.get("exp1")
            exp2 = comp.get("exp2")
            metric = comp.get("metric", "avg")
            n_points = int(comp.get("max_points", max_points))
            if exp1 not in self.indexes or exp2 not in self.indexes:
                continue

            work = self.data[[self.id_col, "SMILES"] + self.indexes[exp1] + self.indexes[exp2]].copy()
            if metric == "sum":
                work["exp1_enrichment"] = work[self.indexes[exp1]].sum(axis=1)
                work["exp2_enrichment"] = work[self.indexes[exp2]].sum(axis=1)
            else:
                work["exp1_enrichment"] = work[self.indexes[exp1]].mean(axis=1)
                work["exp2_enrichment"] = work[self.indexes[exp2]].mean(axis=1)

            work = work[[self.id_col, "SMILES", "exp1_enrichment", "exp2_enrichment"]]
            if len(work) > n_points:
                work = work.nlargest(n_points, "exp1_enrichment")

            work["fingerprint"] = [smiles_to_fingerprint(sm) for sm in work["SMILES"]]
            work = work.dropna(subset=["fingerprint"])
            if len(work) < 2:
                continue

            X = np.array([list(fp) for fp in work["fingerprint"]])
            tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, len(work) - 1))
            coords = tsne.fit_transform(X)
            work["x"] = coords[:, 0]
            work["y"] = coords[:, 1]

            comp_df = pd.concat(
                [
                    work[[self.id_col, "x", "y", "exp1_enrichment"]]
                    .rename(columns={"exp1_enrichment": "enrichment"})
                    .assign(source=exp1),
                    work[[self.id_col, "x", "y", "exp2_enrichment"]]
                    .rename(columns={"exp2_enrichment": "enrichment"})
                    .assign(source=exp2),
                ],
                ignore_index=True,
            )
            comp_df["size"] = np.log1p(comp_df["enrichment"].clip(lower=0)) + 1

            fig = px.scatter(
                comp_df,
                x="x",
                y="y",
                color="source",
                size="size",
                hover_data={self.id_col: True, "enrichment": True, "source": True},
                title=f"Trisynthon Comparison Space ({exp1} vs {exp2})",
                labels={"x": "t-SNE 1", "y": "t-SNE 2"},
            )
            fig.update_layout(plot_bgcolor="white", paper_bgcolor="white")
            fig.write_html(os.path.join(clusters_dir, f"trisynthon_comparison_{exp1}_vs_{exp2}.html"))

            plt.figure(figsize=(10, 8))
            for src in [exp1, exp2]:
                src_df = comp_df[comp_df["source"] == src]
                plt.scatter(src_df["x"], src_df["y"], s=src_df["size"] * 8, alpha=0.6, label=src)
            plt.title(f"Trisynthon Comparison Space ({exp1} vs {exp2})")
            plt.xlabel("t-SNE 1")
            plt.ylabel("t-SNE 2")
            plt.legend()
            plt.tight_layout()
            plt.savefig(
                os.path.join(clusters_dir, f"trisynthon_comparison_{exp1}_vs_{exp2}.svg"),
                dpi=300,
            )
            plt.close()
