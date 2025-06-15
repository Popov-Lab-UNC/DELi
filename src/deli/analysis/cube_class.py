import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.dummy import DummyRegressor, DummyClassifier
from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import accuracy_score, r2_score, confusion_matrix
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdFingerprintGenerator
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np
from tqdm import tqdm
import seaborn as sns
import torch
from torch_geometric.loader import DataLoader
import torch.optim as optim
import torch.nn as nn
from deli.analysis.gnn import Final_Network, smi_to_pyg
from PIL import ImageDraw, ImageFont
from IPython.display import display
from deli.analysis.poly_o import PolyO

class DELi_Cube:
    def __init__(self, data: pd.DataFrame, id_col: str, indexes: dict, control_cols: dict = None, lib_size: int = None, raw_indexes: dict = None):
        """
        Initialize the DELi_Cube object with the provided data, indexes, and control columns.

        Parameters:
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
            raw_indexes (dict, optional): A dictionary mapping experimental IDs to raw index ranges, used for SD calculation.
        """
        self.data = data
        self.id_col = id_col
        self.indexes = indexes
        self.control_cols = control_cols
        self.lib_size = lib_size
        self.raw_indexes = raw_indexes

        if not all(col in self.data.columns for col in ['ID_A', 'ID_B', 'ID_C']):
            first_row_id = self.data[self.id_col].iloc[0]
            parts_count = len(first_row_id.split('-'))

            if parts_count == 2:
                new_cols = self.data[self.id_col].str.split('-', expand=True)
                new_cols.columns = ['ID_A', 'ID_B']
                self.data = pd.concat([self.data.iloc[:, :self.data.columns.get_loc(self.id_col) + 1], new_cols, self.data.iloc[:, self.data.columns.get_loc(self.id_col) + 1:]], axis=1)
            elif parts_count == 3:
                new_cols = self.data[self.id_col].str.split('-', expand=True)
                new_cols.columns = ['ID_A', 'ID_B', 'ID_C']
                self.data = pd.concat([self.data.iloc[:, :self.data.columns.get_loc(self.id_col) + 1], new_cols, self.data.iloc[:, self.data.columns.get_loc(self.id_col) + 1:]], axis=1)
            else:
                raise ValueError("DEL_ID must be split into 2 or 3 parts.")



    def SD_min(self) -> tuple:
        """
        This is an empirical sampling depth minimum meant to be an additional QC check for a given DEL run.
        Functionally it uses the standard that an enrichment minimum of 10 will be reproducible across replicates,
        thus if we evaluate our NSC_max we can determine if the minimum sampling depth is met.
        DOI: https://doi.org/10.1177/2472555218757718

        Returns:
            tuple: (NSC_max_dict, SD_min_dict, sampling_depth_dict) where:
            - NSC_max_dict (dict): Dictionary with the highest NSC value for each experiment.
            - SD_min_dict (dict): Dictionary with the minimum sampling depth required for each experiment.
            - sampling_depth_dict (dict): Dictionary with the total sampling depth for each experiment.
        """
        if self.lib_size is None:
            raise ValueError("Library size must be provided during initialization to run SD_min method.")
        if self.raw_indexes is None:
            raise ValueError("Raw indexes must be provided during initialization to run SD_min method.")

        NSC_max_dict = {}
        SD_min_dict = {}
        sampling_depth_dict = {}

        for exp_name, columns in self.raw_indexes.items():
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

        Returns:
            pd.DataFrame: DataFrame with NSC values for each library member.
        """
        if self.lib_size is None:
            raise ValueError("Library size must be provided during initialization to run NSC_values method.")
        if self.raw_indexes is None:
            raise ValueError("Raw indexes must be provided during initialization to run NSC_values method.")

        for exp_name, columns in self.raw_indexes.items():
            row_sum = self.data[columns].sum(axis=1)
            total_sampling_depth = row_sum.sum() / self.lib_size

            exp_row_sum = self.data[columns].sum(axis=1)
            exp_NSC = exp_row_sum / total_sampling_depth

            # self.data[f"{exp_name}_sum"] = exp_row_sum <- this is raw sum
            self.data[f"{exp_name}_NSC"] = exp_NSC

        return self.data
    
    def NSC_enrichment_intervals(self, nsc_df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate normalized sequence count intervals assuming a Poisson distribution.

        Parameters:
            nsc_df (pd.DataFrame): DataFrame returned from NSC_values method.

        Returns:
            pd.DataFrame: DataFrame with NSC enrichment intervals for each library member.
        """
        if self.lib_size is None:
            raise ValueError("Library size must be provided during initialization to run NSC_enrichment_intervals method.")
        
        if not isinstance(nsc_df, pd.DataFrame):
            raise ValueError("Input must be a DataFrame returned from NSC_values method.")
        
        df = nsc_df.copy()
        for exp_name in self.indexes.keys():
            mean_seq_count = df[f"{exp_name}_sum"].sum(axis=0) / self.lib_size
            df[f"{exp_name}_NSC+"] = ((df[f"{exp_name}_sum"] + 1) ** 0.5 + 1) ** 2 / mean_seq_count
            df[f"{exp_name}_NSC-"] = ((df[f"{exp_name}_sum"] + 1) ** 0.5 - 1) ** 2 / mean_seq_count
        self.data = df
        return self.data

    def maximum_likelihood_enrichment_ratio(self) -> pd.DataFrame:
        """
        Calculate the maximum likelihood enrichment ratio to avoid division by zero for members without sequence reads in control.

        DOI: https://doi.org/10.1021/acsomega.3c02152

        Returns:
            pd.DataFrame: DataFrame with maximum likelihood enrichment ratio for each library member.
        """
        if self.control_cols is None:
            raise ValueError("Control columns must be provided during initialization to run maximum_likelihood_enrichment_ratio method.")
        
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

                df[f"{exp_name}_MLE"] = (
                    (control_total / selection_total) *
                    (((df[columns].sum(axis=1)) + 3/8) / (df[control_col] + 3/8))
                )

        self.data = df
        return self.data


 
   
    def z_score(self):
        """
        Calculate the normalized z-score (zn) for each experimental group compared to the control column.

        The z-score is calculated using the observed count (C) and expected count (E), with 
        standard deviation from a binomial distribution. The final zn is further normalized by sqrt(n).

        DOI: 10.1021/acscombsci.8b00116

        Raises:
            ValueError: If control columns are missing or contain only zeros.

        Returns:
            pd.DataFrame: DataFrame with sum, average, and normalized z-score columns.
        """
        if not self.control_cols:
            raise ValueError("Control columns must be provided during initialization to run z_score method.")

        df = self.data.copy()

        for exp_name, columns in self.indexes.items():
            control_cols = self.control_cols.get(exp_name)

            if control_cols is None:
                raise ValueError(f"Missing control columns for experiment '{exp_name}'.")

            # Compute observed count (C) as the sum across experimental columns
            df[f"{exp_name}_C"] = df[columns].sum(axis=1)

            for control_col in control_cols:
                if control_col not in df.columns:
                    raise ValueError(f"Control column '{control_col}' is missing in data.")

                control_values = df[control_col]
                if control_values.sum() == 0:
                    raise ValueError(f"The control column '{control_col}' has a sum of zero. Calculation cannot be done.")
                else:
                    # If there are non-zero control values, compute the expected count (E)
                    df[f"{exp_name}_E"] = control_values

                    # If there are zero values in the control column, replace them with the median of non-zero control values
                    df[f"{exp_name}_E"] = df[f"{exp_name}_E"].replace(0, control_values[control_values > 0].median())

                n = df[control_col].sum(axis=0)

                if df[f"{exp_name}_C"].sum(axis=0) == 0 or n == 0:
                    raise ValueError(f"Total sum for experiment '{exp_name}' or its control is zero, cannot compute z-score.")

                # binomial standard deviation: σ = sqrt(E * (1 - p_i))
                df[f"{exp_name}_sigma"] = np.sqrt(df[f"{exp_name}_E"] * (1 - (df[f"{exp_name}_E"] / n)))

                # raw z-score: z = (C - E) / σ
                df[f"{exp_name}_z_score"] = (df[f"{exp_name}_C"] - df[f"{exp_name}_E"]) / df[f"{exp_name}_sigma"]

                # normalized z-score: zn = z / sqrt(n)
                df[f"{exp_name}_norm_z_score"] = df[f"{exp_name}_z_score"] / np.sqrt(n)

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

        Returns:
            pd.DataFrame: DataFrame with additional columns for the log-transformed sums, averages, and z-scores for each experimental group.
        """
        if not self.control_cols:
            raise ValueError("Control columns must be provided during initialization to run z_score_log_data method.")

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
                    raise ValueError(f"The control column '{control_col}' has a sum of zero. Calculation cannot be done.")

                df[f"{exp_name}_control_log"] = np.log1p(control_values)
                df[f"{exp_name}_control_log"] = df[f"{exp_name}_control_log"].replace(0, np.log1p(control_values[control_values > 0].median()))

            # mean and standard deviation of the log-transformed control values
            df[f"{exp_name}_control_log_mean"] = df[f"{exp_name}_control_log"].mean()
            df[f"{exp_name}_control_log_std"] = df[f"{exp_name}_control_log"].std()

            # log-transformed z-score while handling zero standard deviation
            df[f"{exp_name}_z_score_log"] = df.apply(
                lambda row: (
                    (row[f"{exp_name}_log"] - row[f"{exp_name}_control_log_mean"]) / row[f"{exp_name}_control_log_std"]
                    if row[f"{exp_name}_control_log_std"] > 0 else np.nan
                ),
                axis=1
            )

            # Handle any missing or problematic z-scores (i.e., divide by zero or invalid calculation)
            if df[f"{exp_name}_z_score_log"].isnull().any():
                raise ValueError(f"Some log-transformed z-scores for experiment '{exp_name}' are invalid (NaN values).")

        self.data = df
        return self.data


    
    def PolyO(self, feature_mode="disynthon") -> pd.DataFrame:
        """
        Calculate PolyO scores and update the DataFrame with the results.

        Returns:
            pd.DataFrame: The updated DataFrame with PolyO scores.
        """
        polyo = PolyO(self.data, self.indexes, self.raw_indexes, self.lib_size, feature_mode=feature_mode)
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

        Returns:
            pd.DataFrame: The DataFrame with normalized columns.
        """
        if self.control_cols is None:
            raise ValueError("Control columns must be provided during initialization to run normalize method.")

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
            raw_col = control_col.replace('corrected', 'raw')
            if raw_col in normalized_df.columns:
                control_cols_to_drop.add(raw_col)

        # Drop all control columns after normalization
        #normalized_df.drop(columns=list(control_cols_to_drop), inplace=True, axis=1)

        self.data = normalized_df
        return self.data


    def disynthonize(self, df: pd.DataFrame = None, synthon_ids: list = None) -> tuple:
        """
        Generate disynthon columns and count data for each experimental group.

        Parameters:
            df (pd.DataFrame, optional): The DataFrame containing the data. Defaults to self.data.
            synthon_ids (list, optional): List of synthon column names. Defaults to ['ID_A', 'ID_B', 'ID_C'].

        Returns:
            tuple: A tuple containing the updated DataFrame and a dictionary mapping experiment names to disynthon count columns.
        """
        if df is None:
            df = self.data

        if synthon_ids is None:
            synthon_ids = ['ID_A', 'ID_B', 'ID_C']

        if len(synthon_ids) not in [2, 3]:
            raise ValueError("synthon_ids must contain exactly two or three elements.")
        
        missing_cols = [col for col in synthon_ids if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing synthon columns: {', '.join(missing_cols)}")
        
        new_cols = {}
        new_cols['AB'] = df[synthon_ids[0]] + '-' + df[synthon_ids[1]]
        if len(synthon_ids) == 3:
            new_cols['AC'] = df[synthon_ids[0]] + '-' + df[synthon_ids[2]]
            new_cols['BC'] = df[synthon_ids[1]] + '-' + df[synthon_ids[2]]

        col_order = list(df.columns)
        if 'DEL_ID' in col_order:
            del_idx = col_order.index('DEL_ID') + 1  
            for col_name, col_values in reversed(new_cols.items()):
                df.insert(del_idx, col_name, col_values)
        else:
            for col_name, col_values in new_cols.items():
                df[col_name] = col_values  

        exp_dict = {}

        for exp_id, index_range in self.indexes.items():
            ab_count = df.groupby('AB')[index_range].sum().reset_index()
            ab_count.columns = ['AB'] + [f'count_AB_{exp_id}_{idx}' for idx in index_range]

            if len(synthon_ids) == 3:
                ac_count = df.groupby('AC')[index_range].sum().reset_index()
                ac_count.columns = ['AC'] + [f'count_AC_{exp_id}_{idx}' for idx in index_range]

                bc_count = df.groupby('BC')[index_range].sum().reset_index()
                bc_count.columns = ['BC'] + [f'count_BC_{exp_id}_{idx}' for idx in index_range]

            df = df.merge(ab_count, on='AB', how='left')
            if len(synthon_ids) == 3:
                df = df.merge(ac_count, on='AC', how='left')
                df = df.merge(bc_count, on='BC', how='left')

            for idx in index_range:
                exp_dict[f'AB_{exp_id}'] = exp_dict.get(f'AB_{exp_id}', []) + [f'count_AB_{exp_id}_{idx}']
                if len(synthon_ids) == 3:
                    exp_dict[f'AC_{exp_id}'] = exp_dict.get(f'AC_{exp_id}', []) + [f'count_AC_{exp_id}_{idx}']
                    exp_dict[f'BC_{exp_id}'] = exp_dict.get(f'BC_{exp_id}', []) + [f'count_BC_{exp_id}_{idx}']

        for count_cols in exp_dict.values():
            for col in count_cols:
                df[col].fillna(0, inplace=True)

        self.data = df
        return df, exp_dict


    def get_top_disynthons(self, disynthon_data: pd.DataFrame, exp_name1: str, comparison_type: str = 'control', exp_name2: str = None, control_name: str = None, comparison_metric: str = 'avg', top_count: int = 20, output_dir: str = '.') -> None:
        """
        Plots a bar chart comparing enrichment between two experimental sets or a single experiment vs control based on a specified metric.
        Creates a separate plot for each synthon type (e.g., AB, BC, AC).

        Parameters:
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

        if comparison_type == 'exp2':
            if not exp_name2:
                raise ValueError("exp_name2 must be provided for comparison_type 'exp2'.")
            if exp_name2 not in self.indexes:
                raise ValueError(f"Experiment '{exp_name2}' not found in the indexes.")

        elif comparison_type == 'control':
            if not control_name:
                raise ValueError("control_name must be provided for comparison_type 'control'.")
            if control_name not in self.control_cols:
                raise ValueError(f"Control column '{control_name}' not found in the DataFrame.")

        elif comparison_type not in ['control', 'exp2', 'none']:
            raise ValueError("comparison_type must be either 'control', 'exp2', or 'none'.")

        if comparison_metric not in ['sum', 'avg']:
            raise ValueError("Comparison metric must be either 'sum' or 'avg'.")

        synthon_ids = []
        if any(col.startswith('AB') for col in disynthon_data.columns):
            synthon_ids.append('AB')
        if any(col.startswith('BC') for col in disynthon_data.columns):
            synthon_ids.append('BC')
        if any(col.startswith('AC') for col in disynthon_data.columns):
            synthon_ids.append('AC')

        if not synthon_ids:
            raise ValueError("No valid synthon columns found in the DataFrame.")

        for synthon in synthon_ids:
            enrichment_results = []

            synthon_cols = [col for col in disynthon_data.columns if col.startswith(synthon)]

            index_cols1 = self.indexes[exp_name1]
            grouped_counts1 = disynthon_data.groupby(synthon)[index_cols1].sum().reset_index()
            grouped_counts1['exp'] = exp_name1

            if comparison_metric == 'avg':
                grouped_counts1['metric'] = grouped_counts1[index_cols1].mean(axis=1)
            else:
                grouped_counts1['metric'] = grouped_counts1[index_cols1].sum(axis=1)

            grouped_counts2 = pd.DataFrame()

            if comparison_type == 'exp2':
                if exp_name2:
                    index_cols2 = self.indexes[exp_name2]
                    grouped_counts2 = disynthon_data.groupby(synthon)[index_cols2].sum().reset_index()
                    grouped_counts2['exp'] = exp_name2

                    if comparison_metric == 'avg':
                        grouped_counts2['metric'] = grouped_counts2[index_cols2].mean(axis=1)
                    else:
                        grouped_counts2['metric'] = grouped_counts2[index_cols2].sum(axis=1)

            elif comparison_type == 'control':
                if control_name:
                    index_cols2 = self.control_cols[control_name]
                    grouped_counts2 = disynthon_data.groupby(synthon)[index_cols2].sum().reset_index()
                    grouped_counts2['exp'] = control_name

                    if comparison_metric == 'avg':
                        grouped_counts2['metric'] = grouped_counts2[index_cols2].mean(axis=1)
                    else:
                        grouped_counts2['metric'] = grouped_counts2[index_cols2].sum(axis=1)

            elif comparison_type == 'none':
                grouped_counts2 = grouped_counts1

            if not grouped_counts2.empty:
                merged_counts = grouped_counts1.merge(grouped_counts2, on=synthon, suffixes=('_exp1', '_exp2'))

                if comparison_type != 'none':
                    merged_counts['diff'] = merged_counts['metric_exp1'] - merged_counts['metric_exp2']
                    top_disynthons = merged_counts.nlargest(top_count, 'diff')
                else:
                    top_disynthons = merged_counts.nlargest(top_count, 'metric_exp1')

                enrichment_results.append(top_disynthons)

            if enrichment_results:
                final_results = pd.concat(enrichment_results)

                fig, ax = plt.subplots(figsize=(12, 6))

                x_labels = final_results[synthon].astype(str)
                x = range(len(x_labels))

                width = 0.35

                if comparison_type == 'none':
                    ax.bar(x, final_results['metric_exp1'], width, label=exp_name1, color='skyblue')
                    ax.set_xticks(x)
                    ax.set_xticklabels(x_labels, rotation=90, fontsize=12)
                    ax.set_xlabel(f"{synthon} Disynthon")
                    ax.set_ylabel("Enrichment")
                    ax.set_title(f"Top {top_count} {synthon} Disynthons: {exp_name1}", fontsize=14, fontweight='bold')
                    ax.legend()

                else:
                    x1 = [pos - width / 2 for pos in x]
                    x2 = [pos + width / 2 for pos in x]

                    ax.bar(x1, final_results['metric_exp1'], width, label=exp_name1, color='skyblue')
                    ax.bar(x2, final_results['metric_exp2'], width, label=exp_name2 if comparison_type == 'exp2' else control_name + '_control', color='orange')

                    ax.set_xticks(x)
                    ax.set_xticklabels(x_labels, rotation=90, fontsize=12)
                    ax.set_xlabel(f"{synthon} Disynthon")
                    ax.set_ylabel("Enrichment")
                    ax.set_title(f"Top {top_count} {synthon} Disynthons: {exp_name1}{' vs ' + (exp_name2 if comparison_type == 'exp2' else control_name + '_control') if comparison_type != 'none' else ''}", fontsize=14, fontweight='bold')
                    ax.legend()

                plt.tight_layout()

                plot_path = f"{output_dir}/top_{synthon}_disynthons_{exp_name1}{' vs ' + (exp_name2 if comparison_type == 'exp2' else control_name) if comparison_type != 'none' else ''}.png"
                plt.savefig(plot_path)
                print(f"Plot saved to {plot_path}")

            else:
                print(f"No valid results for {synthon} disynthon.")

    def simple_spotfire_version(self) -> pd.DataFrame:
        """
        Create a Spotfire-friendly version of the DataFrame using normalized data if available
        with corrected indexes, renaming based on experiment IDs and adding an average column for each experiment.
        Ensures that '_sum' columns are present for each experiment and calculates them if missing.

        Returns:
            pd.DataFrame: The Spotfire-friendly DataFrame.
        """
        df = self.data.copy()
        smiles_index = df.columns.get_loc('SMILES')
        cols_to_keep = df.columns[:smiles_index + 1].tolist() + [col for col in df.columns[smiles_index + 1:] if any(pattern in col for pattern in ['_sum', '_NSC', '_MLE', '_z_score', '_norm_z_score', "_PolyO_score"])]

        for exp_name, index_range in self.indexes.items():
            sum_col_name = f'{exp_name}_sum'
            if sum_col_name not in df.columns:
                df[sum_col_name] = df[index_range].sum(axis=1)
                cols_to_keep.append(sum_col_name)

        # Round relevant columns
        for col in cols_to_keep:
            if any(pattern in col for pattern in ['_sum', '_NSC', '_MLE', '_z_score', '_norm_z_score', "_PolyO_score"]):
                df[col] = df[col].round(2)

        spotfire_df = df[cols_to_keep]

        for exp_name, index_range in self.indexes.items():
            for column in index_range:
                spotfire_df.rename(columns={column: f'{exp_name}_{column}'}, inplace=True)
            avg_col_name = f'{exp_name}_avg'
            if avg_col_name not in spotfire_df.columns:
                spotfire_df[avg_col_name] = df[index_range].mean(axis=1).round(2)

        return spotfire_df

    
    def trisynthon_overlap(self, output_dir=None, normalized_data=None, threshold=0.0):
        """
        Create overlap diagrams for tri-synth experiments using normalized data and corrected indexes.

        Parameters:
            normalized_data (pd.DataFrame, optional): The DataFrame containing normalized data.
               If None, it defaults to using self.data.
            threshold (float): The threshold for filtering counts. Default is 0.0.
        """
       
        if normalized_data is None:
            normalized_data = self.data

        if output_dir is None:
            output_dir = '.'

        valid_experiments = 0  # Track how many experiments have at least two valid indices

        for exp_name, indices in self.indexes.items():
            valid_indices = [idx for idx in indices if idx in normalized_data.columns]
            num_valid_indices = len(valid_indices)

            if num_valid_indices < 2:
                print(f"Skipping experiment '{exp_name}' - not enough valid indices ({num_valid_indices} found).")
                continue
            elif num_valid_indices < len(indices):
                print(f"Warning: Experiment '{exp_name}' has missing indices. Proceeding with {num_valid_indices} valid indices.")

            valid_experiments += 1

            filtered_data = normalized_data[
                (normalized_data[valid_indices[0]] > threshold) |
                (normalized_data[valid_indices[1]] > threshold)
            ]

            set_a = set(filtered_data[filtered_data[valid_indices[0]] > threshold]['DEL_ID'])
            set_b = set(filtered_data[filtered_data[valid_indices[1]] > threshold]['DEL_ID'])

            plt.figure(figsize=(8, 6))
            if num_valid_indices == 3:
                set_c = set(filtered_data[filtered_data[valid_indices[2]] > threshold]['DEL_ID'])
                venn3([set_a, set_b, set_c], set_labels=(valid_indices[0], valid_indices[1], valid_indices[2]))
                plt.title(f'Overlap Diagram for {exp_name} (Three Indices)')
            else:
                venn2([set_a, set_b], set_labels=(valid_indices[0], valid_indices[1]))
                plt.title(f'Overlap Diagram for {exp_name} (Two Indices)')

            plt.savefig(f'{output_dir}/{exp_name}_overlap_diagram.png')
            plt.close()

        if valid_experiments == 0:
            raise ValueError("No experiments had at least two valid indices. Cannot generate overlap diagrams.")



    def disynthon_overlap(self, output_dir: str = None, disynthon_data: pd.DataFrame = None, disynth_exp_dict: dict = None, threshold: float = 0.0) -> None:
        """
        Create overlap diagrams for disynth experiments using normalized data and corrected indexes.

        Parameters:
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
            output_dir = '.'

        valid_experiments = 0  # Track number of valid experiments

        for exp_name, indices in disynth_exp_dict.items():
            # Filter valid indices that exist in the dataframe
            valid_indices = [idx for idx in indices if idx in disynthon_data.columns]
            num_valid_indices = len(valid_indices)

            if num_valid_indices < 2:
                print(f"Skipping experiment '{exp_name}' - not enough valid indices ({num_valid_indices} found).")
                continue
            elif num_valid_indices < len(indices):
                print(f"Warning: Experiment '{exp_name}' has missing indices. Proceeding with {num_valid_indices} valid indices.")

            valid_experiments += 1

            mask = disynthon_data[valid_indices].gt(threshold).any(axis=1)
            filtered_data = disynthon_data[mask]

            disynthon_type = exp_name.split('_')[0]
            if disynthon_type not in filtered_data.columns:
                print(f"Skipping experiment '{exp_name}' - column '{disynthon_type}' not found in the DataFrame.")
                continue

            set_a = set(filtered_data[filtered_data[valid_indices[0]] > threshold][disynthon_type])
            set_b = set(filtered_data[filtered_data[valid_indices[1]] > threshold][disynthon_type])

            plt.figure(figsize=(8, 6))
            if num_valid_indices == 3:
                set_c = set(filtered_data[filtered_data[valid_indices[2]] > threshold][disynthon_type])
                venn3([set_a, set_b, set_c], set_labels=(valid_indices[0], valid_indices[1], valid_indices[2]))
                plt.title(f'Venn Diagram for {exp_name} (Three Indices) \n Threshold = {threshold}', fontsize=16)
            else:
                venn2([set_a, set_b], set_labels=(valid_indices[0], valid_indices[1]))
                plt.title(f'Venn Diagram for {exp_name} (Two Indices) \n Threshold = {threshold}', fontsize=16)

            plt.savefig(f'{output_dir}/{exp_name}_venn_diagram.png', bbox_inches='tight')
            plt.close()

        if valid_experiments == 0:
            raise ValueError("No experiments had at least two valid indices. Cannot generate overlap diagrams.")


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

        Parameters:
            output_dir (str, optional): Directory to save the scatter plot. Defaults to None.

        Returns:
            None
        """

        def bulk_smiles_to_fingerprint(smiles_list, radius=3, n_bits=2048):
            mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius, n_bits)
            return [mfpgen.GetFingerprint(Chem.MolFromSmiles(smiles)) for smiles in smiles_list]

        for exp_name, indices in self.indexes.items():
            self.data[f'{exp_name}_average_enrichment'] = self.data[indices].mean(axis=1)
            top_50 = self.data.nlargest(50, f'{exp_name}_average_enrichment')
            remaining_data = self.data.drop(top_50.index)
            random_50 = remaining_data.sample(n=50, random_state=42)
            subset_data = pd.concat([top_50, random_50]).copy()

            subset_data['fingerprints'] = bulk_smiles_to_fingerprint(subset_data['SMILES'])
            subset_data.dropna(subset=['fingerprints'], inplace=True)
            X = np.array([list(fp) for fp in subset_data['fingerprints']])
            y = subset_data[f'{exp_name}_average_enrichment']

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
            plt.scatter(last_y_test, last_y_pred_rf, label=f"Random Forest (Last Fold R² = {r2_rf_scores[-1]:.2f})", alpha=0.7)
            plt.scatter(last_y_test, last_y_pred_dummy, label=f"Dummy Regressor (Last Fold R² = {r2_dummy_scores[-1]:.2f})", alpha=0.7, marker='x')
            plt.plot([min(last_y_test), max(last_y_test)], [min(last_y_test), max(last_y_test)], color='black', linestyle='--')
            plt.title(f'{exp_name} RF Regression: K-Fold Avg R² = {avg_r2_rf:.2f}, Dummy = {avg_r2_dummy:.2f}')
            plt.xlabel('True Average Enrichment')
            plt.ylabel('Predicted Average Enrichment')
            plt.legend(loc="upper left")

            if output_dir is None:
                output_dir = '.'
            plt.savefig(f'{output_dir}/{exp_name}_RF_regression.png')
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

        Parameters:
            threshold (int, optional): Threshold value for enrichment status. Defaults to 2.
            output_dir (str, optional): Directory to save the classifier image. Defaults to None.

        Returns:
            None
        """
        def bulk_smiles_to_fingerprint(smiles_list, radius=3, n_bits=2048):
            if isinstance(smiles_list, str):
                smiles_list = [smiles_list]
            mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius, n_bits)
            return [mfpgen.GetFingerprint(Chem.MolFromSmiles(smiles)) for smiles in tqdm(smiles_list)]

        for exp_name, indices in self.indexes.items():
            self.data[f'{exp_name}_average_enrichment'] = self.data[indices].mean(axis=1)
            top_50 = self.data.nlargest(100, f'{exp_name}_average_enrichment')
            remaining_data = self.data.drop(top_50.index)
            random_50 = remaining_data.sample(n=100, random_state=42)
            subset_data = pd.concat([top_50, random_50]).copy()
            subset_data['fingerprints'] = bulk_smiles_to_fingerprint(subset_data['SMILES'])
            subset_data.dropna(subset=['fingerprints'], inplace=True)
            subset_data['target'] = (subset_data[f'{exp_name}_average_enrichment'] > threshold).astype(int)
            X = np.array([list(fp) for fp in subset_data['fingerprints']])
            y = subset_data['target']
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
            sns.heatmap(cm_rf, annot=True, fmt='d', cmap='Blues', xticklabels=["Not Enriched", "Enriched"], yticklabels=["Not Enriched", "Enriched"])
            plt.title(f"{exp_name} Random Forest Classifier (Avg Acc = {avg_acc_rf:.2f}), Threshold = {threshold} \n"
                  f"Last Fold Confusion Matrix Accuracy= {acc_rf_scores[-1]:.2f}")
            plt.xlabel('Predicted')
            plt.ylabel('True')
            plt.subplot(1, 2, 2)
            sns.heatmap(cm_dummy, annot=True, fmt='d', cmap='Blues', xticklabels=["Not Enriched", "Enriched"], yticklabels=["Not Enriched", "Enriched"])
            plt.title(f"{exp_name} Dummy Classifier (Avg Acc = {avg_acc_dummy:.2f}), Threshold = {threshold} \n"
                  f"Last Fold Confusion Matrix Accuracy= {acc_dummy_scores[-1]:.2f}")
            plt.xlabel('Predicted')
            plt.ylabel('True')
            plt.tight_layout()

            if output_dir is None:
                output_dir = '.'
            plt.savefig(f'{output_dir}/{exp_name}_classifier.png')
            plt.close()

    def gnn_classifier(self, threshold: int = 2, output_dir: str = ".", arch: str = "GAT", num_layers: int = 3, encoding: str = "embedding") -> None:
        """
        Train a Graph Neural Network (GNN) classifier to predict the enrichment status of compounds.

        Parameters:
            threshold (int, optional): Threshold value for enrichment status, by default 2.
            output_dir (str, optional): Directory to save the classifier image, by default ".".
            arch (str, optional): GNN architecture (GAT or GCN), by default "GAT".
            num_layers (int, optional): Number of GNN layers, by default 3.
            encoding (str, optional): Node encoding method (embedding or onehot), by default "embedding".
        """
        
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        for exp_name, indices in self.indexes.items():
            self.data[f'{exp_name}_average_enrichment'] = self.data[indices].mean(axis=1)

            top_100 = self.data.nlargest(100, f'{exp_name}_average_enrichment')
            remaining_data = self.data.drop(top_100.index)
            random_200 = remaining_data.sample(n=200, random_state=42)

            subset_data = pd.concat([top_100, random_200]).copy()
            subset_data["graphs"] = [smi_to_pyg(smiles, 0, encoding="embedding") for smiles in tqdm(subset_data['SMILES'])]
            subset_data.dropna(subset=["graphs"], inplace=True)
            subset_data["target"] = (subset_data[f'{exp_name}_average_enrichment'] > threshold).astype(int)

            graphs = list(subset_data["graphs"])
            labels = torch.tensor(subset_data["target"].values, dtype=torch.long)

            train_graphs, test_graphs, train_labels, test_labels = train_test_split(graphs, labels, test_size=0.2, random_state=42)
            train_graphs, val_graphs, train_labels, val_labels = train_test_split(train_graphs, train_labels, test_size=0.1, random_state=42)


            for i, g in enumerate(train_graphs):
                g.y = train_labels[i]
            for i, g in enumerate(test_graphs):
                g.y = test_labels[i]

            train_loader = DataLoader(train_graphs, batch_size=10, shuffle=True)
            val_loader = DataLoader(val_graphs, batch_size=10, shuffle=False)
            test_loader = DataLoader(test_graphs, batch_size=10, shuffle=False)

            node_dim, edge_dim = train_graphs[0].x.shape[1], train_graphs[0].edge_attr.shape[1]
            model = Final_Network(node_dim, edge_dim, arch=arch, num_layers=num_layers, encoding=encoding).to(device)
            optimizer = optim.AdamW(model.parameters(), lr=1e-5)
            scheduler = optim.lr_scheduler.OneCycleLR(optimizer, max_lr=1e-3, steps_per_epoch=len(train_loader), epochs=200)
            criterion = nn.CrossEntropyLoss()

            train_losses = []
            val_losses = []

            best_val_loss = float('inf')

            for epoch in range(200):
                model.train()
                epoch_loss = 0
                for batch in train_loader:
                    batch = batch.to(device)
                    optimizer.zero_grad()
                    out = model(batch.x, batch.edge_index, batch.edge_attr, batch.batch)
                    loss = criterion(out, batch.y)
                    loss.backward()
                    optimizer.step()
                    scheduler.step()
                    epoch_loss += loss.item()
                train_losses.append(epoch_loss / len(train_loader))

                model.eval()
                val_loss = 0
                with torch.no_grad():
                    for batch in val_loader:
                        batch = batch.to(device)
                        out = model(batch.x, batch.edge_index, batch.edge_attr, batch.batch)
                        loss = criterion(out, batch.y)
                        val_loss += loss.item()
                val_loss /= len(val_loader)
                val_losses.append(val_loss)

                # Early stopping (optional)
                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    torch.save(model.state_dict(), f"{output_dir}/{exp_name}_{arch}_best_model.pth")


            # # Plot Training Loss Curve
            # plt.figure(figsize=(6, 5))
            # plt.plot(range(1, 201), train_losses, marker='o', linestyle='-')
            # plt.xlabel("Epoch")
            # plt.ylabel("Training Loss")
            # plt.title(f"{exp_name} {arch} Training Loss (Fold {fold+1})")
            # plt.grid()
            # plt.savefig(f'{output_dir}/{exp_name}_{arch}_loss_fold{fold+1}.png')
            # plt.close()

            model.load_state_dict(torch.load(f"{output_dir}/{exp_name}_{arch}_best_model.pth"))
            model.eval()
            y_true, y_pred = [], []

            with torch.no_grad():
                for batch in test_loader:
                    batch = batch.to(device)
                    out = model(batch.x, batch.edge_index, batch.edge_attr, batch.batch)
                    preds = torch.argmax(out, dim=1)
                    y_true.extend(batch.y.cpu().numpy())
                    y_pred.extend(preds.cpu().numpy())


            acc = accuracy_score(y_true, y_pred)

            cm = confusion_matrix(y_true, y_pred)

            plt.figure(figsize=(8, 6))
            sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=["Not Enriched", "Enriched"], yticklabels=["Not Enriched", "Enriched"])
            plt.title(f"{exp_name} {arch} Classifier (Acc = {acc:.2f}), Threshold = {threshold}")
            plt.xlabel('Predicted')
            plt.ylabel('True')
            plt.tight_layout()
            plt.savefig(f'{output_dir}/{exp_name}_{arch}_classifier.png')
            plt.close()



    def top_n_compounds(self, n: int = 20, metric: str = 'sum', output_dir: str = None) -> None:
        """
        Display the top N compounds based on a specified metric and save the image to the specified directory.

        Parameters:
            n (int, optional): Number of top compounds to display. Defaults to 20.
            metric (str, optional): Metric to rank the compounds. Defaults to 'sum'.
            output_dir (str, optional): Directory to save the image. Defaults to None (current directory).
        """
        if metric == 'avg':
            for exp_name, index_range in self.indexes.items():
                avg_col_name = f'{exp_name}_avg'
                if avg_col_name not in self.data.columns:
                    self.data[avg_col_name] = self.data[index_range].mean(axis=1).round(2)

        for exp_name, indices in self.indexes.items():
            # Ensure the compound is present at least once in each replicate
            filtered_data = self.data[(self.data[indices] > 0).all(axis=1)]
            top_n = filtered_data.nlargest(n, f'{exp_name}_{metric}')
            mols = [Chem.MolFromSmiles(smiles) for smiles in top_n['SMILES']]
            legends = [f"{row['DEL_ID']}\n\n{metric}: {row[f'{exp_name}_{metric}']:.2f}" for _, row in top_n.iterrows()]
            
            dopts = rdMolDraw2D.MolDrawOptions()
            dopts.legendFontSize = 20

            try:
                img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(200, 200), legends=legends, drawOptions=dopts)
                # Add filename as text on top
                draw = ImageDraw.Draw(img)
                file_name = f'{exp_name}_top_{n}_compounds.png'
                draw.text((10, 10), file_name[:-4], fill="black")
            except AttributeError:
                pass

            if output_dir is None:
                output_dir = '.'
            file_name = f'{exp_name}_top_{n}_compounds.png'
            
            try:
                img.save(f'{output_dir}/{file_name}')
            except AttributeError:
                pass

            try:
                print(f"Top {n} compounds for {exp_name} based on {metric}")
                display(img)
            except ImportError:
                pass


    def positive_control_finder(self, positive_control_ID = None):
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
        

        
        #A048-B080-C096
        pass


    def top_trisynthons_diff_from_competitor(self):
        pass


    def top_disynthons_diff_from_competitor(self):
        pass

    


