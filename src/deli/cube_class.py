import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.dummy import DummyRegressor, DummyClassifier
from sklearn.metrics import accuracy_score, r2_score, confusion_matrix
from sklearn.model_selection import KFold
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np
from tqdm import tqdm
import seaborn as sns

class DELi_Cube:
    def __init__(self, data, indexes, control_cols=None, lib_size=None):
        """
        Initialize the DELi_Cube object with the provided data, indexes, and control columns.

        Parameters:
            data (pd.DataFrame): The DataFrame containing the data.
            indexes (dict): A dictionary mapping experimental IDs to index ranges.
                ex. tri_synth_dict = {
                'protein_replicates': ['corrected_index3', 'corrected_index4', 'corrected_index6'],
                'protein_withinhibitor_replicates': ['corrected_index7', 'corrected_index8', 'corrected_index9']
                }
            control_cols (dict): A dictionary mapping experimental IDs to control columns.
                ex. control_cols = {
                'protein_replicates': 'control_col1',
                'protein_withinhibitor_replicates': 'control_col2'
                }
        """

        self.data = data
        self.indexes = indexes 
        self.control_cols = control_cols 
        self.lib_size = lib_size

    def SD_min(self):
        """This is an empirical sampling depth minimum meant to be an additional QC check for a given DEL run.
            Functionally it uses the standard that an enrichment minimum of 10 will be reproducible across replicates,
            thus if we evaluate our NSC_max we can determine if the minimum sampling depth is met.
            DOI: https://doi.org/10.1177/2472555218757718

            Returns:
                tuple: (NSC_max, SD_min) where:
                    - NSC_max is the highest NSC value in the dataset.
                    - SD_min is the minimum sampling depth required to discover enriched compounds using the 
                    empirical threshold of 10.
        """
        if self.lib_size is None:
            raise ValueError("Library size must be provided during initialization to run SD_min method.")

        NSC_max_dict = {}
        SD_min_dict = {}
        sampling_depth_dict = {}

        for exp_name, columns in self.indexes.items():
            row_avg = self.data[columns].mean(axis=1)
            total_sampling_depth = row_avg.sum() / self.lib_size
            exp_NSC = row_avg / total_sampling_depth
            exp_NSC_max = exp_NSC.max()
            NSC_max_dict[f"{exp_name}_NSC_max"] = exp_NSC_max
            SD_min_dict[f"{exp_name}_SD_min"] = 10 / exp_NSC_max
            sampling_depth_dict[f"{exp_name}_sampling_depth"] = total_sampling_depth

        return NSC_max_dict, SD_min_dict, sampling_depth_dict

    def NSC_values(self):
        """Enrichment factor using normalized counts for given library member normalized to mean sequence reads
           for all library members (DOI: http://dx.doi.org/10.1002/anie.201410736)

        Returns:
            pd.DataFrame: DataFrame with NSC values for each library member
        """
        if self.lib_size is None:
            raise ValueError("Library size must be provided during initialization to run NSC_values method.")
        
        df = self.data.copy()
        for exp_name, columns in self.indexes.items():
            df[f"{exp_name}_sum"] = df[columns].sum(axis=1)
            df[f"{exp_name}_NSC"] = df[f"{exp_name}_sum"] / (df[f"{exp_name}_sum"].sum(axis=0)/ self.lib_size)
        self.data = df
        return self.data
    
    def NSC_enrichment_intervals(self, nsc_df):
        """Normalized sequence count interval assuming a poisson distribution. Calculated using
           the equation: NSC +/- = (sqrt(observed count +1) +/- 1)**2)/(mean sequence count) for all library members

        Parameters:
            nsc_df (pd.DataFrame): DataFrame returned from NSC_values method.

        Returns:
            pd.DataFrame: DataFrame with NSC enrichment intervals for each library member
        """

        if self.lib_size is None:
            raise ValueError("Library size must be provided during initialization to run NSC_enrichment_intervals method.")
        
        if not isinstance(nsc_df, pd.DataFrame):
            raise ValueError("Input must be a DataFrame returned from NSC_values method.")
        
        df = nsc_df.copy()
        for exp_name in self.indexes.keys():
            df[f"{exp_name}_NSC+"] = ((df[f"{exp_name}_sum"] + 1) ** 0.5 + 1) ** 2 / (df[f"{exp_name}_sum"].sum(axis=0)/ self.lib_size)
            df[f"{exp_name}_NSC-"] = ((df[f"{exp_name}_sum"] + 1) ** 0.5 - 1) ** 2 / (df[f"{exp_name}_sum"].sum(axis=0)/ self.lib_size)
        self.data = df
        return self.data

    def maximum_likelihood_enrichment_ratio(self):
        """Formulation of enrichment ratio that avoids division by zero for members without sequence reads
           in control. doi: https://doi.org/10.1021/acsomega.3c02152?

           Calculated using: (n_controlsequences/n_selection_sequences) * (c_selection_reads + 3/8)/(c_control_reads + 3/8)
           where c is the number of sequence reads for a given library member and n is total number of sequence reads for all library members.
        
        Returns:
            pd.DataFrame: DataFrame with maximum likelihood enrichment ratio for each library member
        """
        if self.control_cols is None:
            raise ValueError("Control columns must be provided during initialization to run maximum_likelihood_enrichment_ratio method.")
        
        df = self.data.copy()
        for exp_name, columns in self.indexes.items():
            control_col = self.control_cols.get(exp_name)
            if control_col is None:
                raise ValueError(f"Control column for {exp_name} not provided.")
            df[f"{exp_name}_sum"] = df[columns].sum(axis=1)
            df[f"{exp_name}_MLE"] = (df[control_col].sum(axis=0)) / (df[f"{exp_name}_sum"].sum(axis=0)) * ((df[columns].sum(axis=1) + 3/8) / (df[control_col] + 3/8))
        self.data = df
        return self.data

    def z_score(self):
        """Calculate the z-score for each experimental group compared to the control column.
        The z-score is calculated for each experimental group by comparing the proportion of 
        observations in the experimental group to the proportion of observations in the control group,
        using MAD instead of standard deviation for scaling.
        DOI: 10.1021/acscombsci.8b00116

        Raises:
            ValueError: If the control column is not provided during initialization.
            pandas.DataFrame: A DataFrame with additional columns for the sum and z-score of each experimental group.

        Returns:
            pd.DataFrame: DataFrame with z-score values for each library member
        """

        if self.control_cols is None:
            raise ValueError("Control columns must be provided during initialization to run z_score method.")
        
        df = self.data.copy()

        for exp_name, columns in self.indexes.items():
            control_col = self.control_cols.get(exp_name)
            if control_col is None:
                raise ValueError(f"Control column for {exp_name} not provided.")
            df[f"{exp_name}_sum"] = df[columns].sum(axis=1)
            df[f"{exp_name}_avg"] = df[f"{exp_name}_sum"] / len(columns)
            
            p_o = df[f"{exp_name}_avg"] / df[f"{exp_name}_avg"].sum(axis=0)
            p_i = df[control_col] / df[control_col].sum(axis=0)
            
            # Filter out zero control values due to NTC having lots of 0's
            non_zero_control = df[control_col].loc[df[control_col] > 0]
            median_control = np.median(non_zero_control)

            # (absolute deviation from the median with scaling factor)
            mad = (1.4286 * np.median(np.abs(non_zero_control - median_control)))
            
            # calculate the z-score using this MAD
            df[f"{exp_name}_z_score"] = (p_i - p_o) / mad
            
        self.data = df
        return self.data

    def z_score_log_data(self):
        """
        Preprocesses the data by calculating the log-transformed z-scores for each experimental group.
        This method performs the following steps:
        1. Creates a copy of the original data.
        2. For each experimental group, calculates the sum and average of the specified columns.
        3. Applies a log transformation to the average values of the experimental group and the control group.
        4. Computes the mean and standard deviation of the log-transformed control group.
        5. Calculates the log-transformed z-scores for the experimental group based on the control group's statistics.
        
        Returns
        -------
        pd.DataFrame
            A DataFrame with additional columns for the log-transformed sums, averages, and z-scores for each experimental group.
        """
        df = self.data.copy()

        if self.control_cols is None:
            raise ValueError("Control columns must be provided during initialization to run z_score method.")
        
        for exp_name, columns in self.indexes.items():
            control_col = self.control_cols.get(exp_name)
            if control_col is None:
                raise ValueError(f"Control column for {exp_name} not provided.")
            df[f"{exp_name}_sum"] = df[columns].sum(axis=1)
            df[f"{exp_name}_avg"] = df[f"{exp_name}_sum"] / len(columns)

            # Log-transform the experimental group and control group
            df[f"{exp_name}_log"] = np.log1p(df[f"{exp_name}_avg"])  # log(x + 1) to avoid log(0)
            df[f"control_log"] = np.log1p(df[control_col])  # Same for control group
            
            # log-transformed mean (lambda) and standard deviation for control group
            control_mean_log = np.mean(df[f"control_log"])
            control_sd_log = np.std(df[f"control_log"])
            
            # Calculate the log-transformed z-score 
            df[f"{exp_name}_z_score_log"] = (df[f"{exp_name}_log"] - control_mean_log) / control_sd_log

        self.data = df
        return self.data


    def normalize(self):
        """
        Normalize specified columns by subtracting the control column in the DataFrame,
        ensuring no values go below zero.

        Returns:
            pd.DataFrame: The DataFrame with normalized columns.
        """

        if self.control_cols is None:
            raise ValueError("Control columns must be provided during initialization to run normalize method.")

        normalized_df = self.data.copy()

        # Iterating over each experimental group in the dictionary
        for exp_name, columns in self.indexes.items():
            control_col = self.control_cols.get(exp_name)
            if control_col is None:
                raise ValueError(f"Control column for {exp_name} not provided.")
            print(f"Normalizing {exp_name} with columns: {columns}")
            
            # Loop through each column in the experimental group
            for column in columns:
                if column in normalized_df.columns:
                    normalized_df[column] = (normalized_df[column] - normalized_df[control_col]).clip(lower=0)
                else:
                    raise ValueError(f"Column '{column}' not found in the DataFrame")
       
            # normalized_df.drop(columns=list(self.control_cols.values()), inplace=True)
        self.data = normalized_df
        return self.data
        

    def disynthonize(self, df=None, synthon_ids=None):
        """
        Create disynthon pairs, count occurrences, and construct a dictionary for disynthon overlap analysis.

        Parameters:
            df (pd.DataFrame): The input DataFrame containing synthon columns and indices. 
            - Defaults to self.data if not provided.
            synthon_ids (list): A list of synthon ID column names (e.g., ['ID_A', 'ID_B'] or ['ID_A', 'ID_B', 'ID_C']).
            - Defaults to 3-cycle ['ID_A', 'ID_B', 'ID_C'].

        Returns:
            pd.DataFrame: A DataFrame with added count columns for each disynthon pair.
            dict: A dictionary mapping disynthon types and experimental IDs to disynthon count columns.
        """

        if df is None:
            df = self.data

        if synthon_ids is None:
            synthon_ids = ['ID_A', 'ID_B', 'ID_C']

        if len(synthon_ids) not in [2, 3]:
            raise ValueError("synthon_ids must contain exactly two or three elements.")

        
        if len(synthon_ids) == 2:
            df['AB'] = df[synthon_ids[0]] + '-' + df[synthon_ids[1]]
        elif len(synthon_ids) == 3:
            df['AB'] = df[synthon_ids[0]] + '-' + df[synthon_ids[1]]
            df['AC'] = df[synthon_ids[0]] + '-' + df[synthon_ids[2]]
            df['BC'] = df[synthon_ids[1]] + '-' + df[synthon_ids[2]]

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

            # dynamically build exp_dict based on exp_id and index ranges for future fxns involving disynth
            for idx in index_range:
                exp_dict[f'AB_{exp_id}'] = exp_dict.get(f'AB_{exp_id}', []) + [f'count_AB_{exp_id}_{idx}']
                if len(synthon_ids) == 3:
                    exp_dict[f'AC_{exp_id}'] = exp_dict.get(f'AC_{exp_id}', []) + [f'count_AC_{exp_id}_{idx}']
                    exp_dict[f'BC_{exp_id}'] = exp_dict.get(f'BC_{exp_id}', []) + [f'count_BC_{exp_id}_{idx}']

        # Used to fill NaN values with 0 in the count columns
        for count_cols in exp_dict.values():
            for col in count_cols:
                df[col].fillna(0, inplace=True)

        return df, exp_dict



    def simple_spotfire_version(self):
        """
        Create a Spotfire-friendly version of the DataFrame using normalized data if available
        with corrected indexes, renaming based on experiment IDs.

        Returns:
            pd.DataFrame: The Spotfire-friendly DataFrame.
        """
        if self.control_col is not None:
            normalized_data = self.normalize()
        else:
            normalized_data = self.data

           
        corrected_indexes = self.indexes
    
        spotfire_df = normalized_data.copy()

        for exp_id, index_range in self.indexes.items():
            for column in index_range:  #
                spotfire_df.rename(columns={column: f'{exp_id}_{column}'}, inplace=True)
        
        if 'SMILES' in self.data.columns:
            spotfire_df['SMILES'] = self.data['SMILES']
        
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

        for exp_name, indices in self.indexes.items():

            num_indices = len(indices)

            
            if num_indices not in [2, 3]:
                raise ValueError(f"Expected exactly two or three indices for '{exp_name}', got {num_indices}")

            
            filtered_data = normalized_data[
                (normalized_data[indices[0]] > threshold) |
                (normalized_data[indices[1]] > threshold) |
                (normalized_data[indices[2]] > threshold) if num_indices == 3 else
                (normalized_data[indices[0]] > threshold) |
                (normalized_data[indices[1]] > threshold)
            ]

            
            set_a = set(filtered_data[filtered_data[indices[0]] > threshold]['DEL_ID'])
            set_b = set(filtered_data[filtered_data[indices[1]] > threshold]['DEL_ID'])
            if num_indices == 3:
                set_c = set(filtered_data[filtered_data[indices[2]] > threshold]['DEL_ID'])

           
            plt.figure(figsize=(8, 6))
            if num_indices == 3:
                venn3([set_a, set_b, set_c], set_labels=(indices[0], indices[1], indices[2]))
                plt.title(f'Overlap Diagram for {exp_name} (Three Indices)')
            else:
                venn2([set_a, set_b], set_labels=(indices[0], indices[1]))
                plt.title(f'Overlap Diagram for {exp_name} (Two Indices)')
            if output_dir is None:
                output_dir = '.'
            plt.savefig(f'{output_dir}/{exp_name}_overlap_diagram.png')


    def disynthon_overlap(self, output_dir=None, disynthon_data=None, disynth_exp_dict=None, threshold=0.0):
        """
        Create overlap diagrams for disynth experiments using normalized data and corrected indexes.

        Parameters:
            disynthon_data (pd.DataFrame): The DataFrame containing disynthon data.
            disynth_exp_dict (dict): A dictionary mapping experiment names to the corresponding corrected index columns.
            threshold (float): The threshold for filtering counts. Default is 0.0.

        Example of disynth_exp_dict:
            disynth_exp_dict = {
                'AB_NSP2': ['count_AB_corrected_index3', 'count_AB_corrected_index4', 'count_AB_corrected_index6'],
                'AB_NSP2_2034': ['count_AB_corrected_index7', 'count_AB_corrected_index8', 'count_AB_corrected_index9'],
            }
        """
        
        if disynthon_data is None or disynth_exp_dict is None:
            raise ValueError("disynthon_data and disynth_exp_dict must be provided.")

        disynthon_data = disynthon_data.copy()

        for exp_name, indices in disynth_exp_dict.items():
           
            if not all(col in disynthon_data.columns for col in indices):
                raise ValueError(f"One or more columns in {indices} are missing from the DataFrame.")

            
            mask = disynthon_data[indices].gt(threshold).any(axis=1)
            filtered_data = disynthon_data[mask]

            #
            disynthon_type = exp_name.split('_')[0]
            if disynthon_type not in filtered_data.columns:
                raise ValueError(f"Column '{disynthon_type}' not found in the DataFrame.")

            
            set_a = set(filtered_data[filtered_data[indices[0]] > threshold][disynthon_type])
            set_b = set(filtered_data[filtered_data[indices[1]] > threshold][disynthon_type])
            
            if len(indices) == 2:
                plt.figure(figsize=(8, 6))
                venn2([set_a, set_b], set_labels=(indices[0], indices[1]))
                plt.title(f'Venn Diagram for {exp_name}', fontsize=16)
                plt.xlabel('Sets', fontsize=14)
                plt.ylabel('Counts', fontsize=14)
                plt.xticks(fontsize=12)
                plt.yticks(fontsize=12)
                venn3_circles = plt.gca().findobj(match=plt.Text)
                for text in venn3_circles:
                    text.set_fontsize(18)  # Increase font size for numbers in circles and overlaps
                if output_dir is None:
                    output_dir = '.'
                plt.savefig(f'{output_dir}/{exp_name}_venn_diagram.png')

            elif len(indices) == 3:
                set_c = set(filtered_data[filtered_data[indices[2]] > threshold][disynthon_type])
                plt.figure(figsize=(8, 6))
                venn3([set_a, set_b, set_c], set_labels=(indices[0], indices[1], indices[2]))
                plt.title(f'Venn Diagram for {exp_name}', fontsize=16)
                plt.xlabel('Sets', fontsize=14)
                plt.ylabel('Counts', fontsize=14)
                plt.xticks(fontsize=12)
                plt.yticks(fontsize=12)
                venn3_circles = plt.gca().findobj(match=plt.Text)
                for text in venn3_circles:
                    text.set_fontsize(18)  # Increase font size for numbers in circles and overlaps
                if output_dir is None:
                    output_dir = '.'
                plt.savefig(f'{output_dir}/{exp_name}_venn_diagram.png')
                
            else:
                raise ValueError(f"Only 2 or 3 indices are supported for Venn diagrams, got {len(indices)}.")

    def ml_fingerprints_to_RF(self, output_dir=None):
        """Trains a Random Forest regressor and a Dummy regressor on trisynthon SMILES fingerprints 
        to predict the average normalized enrichment using a subset of 100 molecules.
        
        The method performs the following steps:
        1. Converts SMILES strings to Morgan fingerprints.
        2. Selects the top 50 molecules with the highest average enrichment and 50 random molecules.
        3. Trains the models using 5-fold cross-validation.
        4. Evaluates the models using R² scores.
        5. Plots the true vs. predicted average enrichment for the last fold.

        Returns:
            None. Saves a scatter plot of the true vs. predicted average enrichment for the last fold.
        """

        def smiles_to_fingerprint(smiles, radius=3, n_bits=2048):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                try:
                    Chem.SanitizeMol(mol)
                    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
                except:
                    return None
            else:
                return None

        for exp_name, indices in self.indexes.items():
            self.data[f'{exp_name}_average_enrichment'] = self.data[indices].mean(axis=1)
            top_50 = self.data.nlargest(50, f'{exp_name}_average_enrichment')
            remaining_data = self.data.drop(top_50.index)
            random_50 = remaining_data.sample(n=50, random_state=42)
            subset_data = pd.concat([top_50, random_50]).copy()

            subset_data['fingerprints'] = [smiles_to_fingerprint(smiles) for smiles in tqdm(subset_data['SMILES'])]
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
        

    def ml_fingerprints_to_classifier(self, threshold=2, output_dir=None):
        def smiles_to_fingerprint(smiles, radius=3, n_bits=2048):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                try:
                    Chem.SanitizeMol(mol)
                    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
                except:
                    return None
            else:
                return None

        for exp_name, indices in self.indexes.items():
            self.data[f'{exp_name}_average_enrichment'] = self.data[indices].mean(axis=1)
            top_50 = self.data.nlargest(100, f'{exp_name}_average_enrichment')
            remaining_data = self.data.drop(top_50.index)
            random_50 = remaining_data.sample(n=100, random_state=42)
            subset_data = pd.concat([top_50, random_50]).copy()
            subset_data['fingerprints'] = [smiles_to_fingerprint(smiles) for smiles in tqdm(subset_data['SMILES'])]
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



    def top_n_compounds(self, n=20, metric='sum', output_dir=None):
        """Display the top N compounds based on a specified metric and save the image to the specified directory.

        Parameters
        ----------
        n : int, optional
            Number of top compounds to display, by default 20
        metric : str, optional
            Metric to rank the compounds, by default 'average_enrichment'
        output_dir : str, optional
            Directory to save the image, by default None (current directory)
        """
        for exp_name, indices in self.indexes.items():
            top_n = self.data.nlargest(n, f'{exp_name}_{metric}')
            mols = [Chem.MolFromSmiles(smiles) for smiles in top_n['SMILES']]
            legends = [f"{row['DEL_ID']}\n\n{metric}: {row[f'{exp_name}_{metric}']:.2f}" for _, row in top_n.iterrows()]
            
            dopts = rdMolDraw2D.MolDrawOptions()
            dopts.legendFontSize = 20

            img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(200, 200), legends=legends, drawOptions=dopts)
            
            if output_dir is None:
                output_dir = '.'
            img.save(f'{output_dir}/{exp_name}_top_{n}_compounds.png')



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

    


