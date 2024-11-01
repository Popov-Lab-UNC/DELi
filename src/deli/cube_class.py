import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from sklearn.ensemble import RandomForestRegressor
from sklearn.dummy import DummyRegressor
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, ConfusionMatrixDisplay
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

class DELi_Cube:
    def __init__(self, data, indexes, control_col=None):
        """
        Initialize the DELi_Cube object with the provided data, indexes, and control column.

        Parameters:
            data (pd.DataFrame): The DataFrame containing the data.
            indexes (dict): A dictionary mapping experimental IDs to index ranges.
                ex. tri_synth_dict = {
                'protein_replicates': ['corrected_index3', 'corrected_index4', 'corrected_index6'],
                'protein_withinhibitor_replicates': ['corrected_index7', 'corrected_index8', 'corrected_index9']
                }
            control_col (str): The control column to use for normalization.
        """

        self.data = data
        self.indexes = indexes 
        self.control_col = control_col 

    def normalize(self):
        """
        Normalize specified columns by subtracting the control column in the DataFrame,
        ensuring no values go below zero.

        Returns:
            pd.DataFrame: The DataFrame with normalized columns.
        """

        if self.control_col is None:
            raise ValueError("Control column must be provided during initialization to run normalize method.")

        normalized_df = self.data.copy()

        # Iterating over each experimental group in the dictionary
        for exp_name, columns in self.indexes.items():
            print(f"Normalizing {exp_name} with columns: {columns}")
            
            # Loop through each column in the experimental group
            for column in columns:
                if column in normalized_df.columns:
                    normalized_df[column] = (normalized_df[column] - normalized_df[self.control_col]).clip(lower=0)
                else:
                    raise ValueError(f"Column '{column}' not found in the DataFrame")
       
        return normalized_df.drop(columns=[self.control_col])
        

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
    


    def trisynthon_overlap(self, normalized_data=None, threshold=0.0):
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
            plt.show()

    def disynthon_overlap(self, disynthon_data, disynth_exp_dict, threshold=0.0):
        """
        Create overlap diagrams for disynth experiments using normalized data and corrected indexes.

        Parameters:
            disynthon_data (pd.DataFrame): The DataFrame containing disynthon data.
            exp_dict (dict): A dictionary mapping experiment names to the corresponding corrected index columns.
            threshold (float): The threshold for filtering counts. Default is 0.0.

        Example of exp_dict:
            exp_dict = {
                'AB_NSP2': ['count_AB_corrected_index3', 'count_AB_corrected_index4', 'count_AB_corrected_index6'],
                'AB_NSP2_2034': ['count_AB_corrected_index7', 'count_AB_corrected_index8', 'count_AB_corrected_index9'],
            }
        """

        disynthon_data = disynthon_data.copy()

        for exp_name, indices in disynth_exp_dict.items():
            filtered_data = disynthon_data[(disynthon_data[indices[0]] > threshold) | 
                                        (disynthon_data[indices[1]] > threshold) | 
                                        (disynthon_data[indices[2]] > threshold)]


            set_a = set(filtered_data[filtered_data[indices[0]] > threshold]['DEL_ID'])
            set_b = set(filtered_data[filtered_data[indices[1]] > threshold]['DEL_ID'])
            set_c = set(filtered_data[filtered_data[indices[2]] > threshold]['DEL_ID'])


            plt.figure(figsize=(8, 6))
            venn3([set_a, set_b, set_c], set_labels=(indices[0], indices[1], indices[2]))
            plt.title(f'Venn Diagram for {exp_name}')
            plt.show()


    def ml_fingerprints_to_RF(self, experiment_key, test_size=None):
        """
        Train a Random Forest regressor and Dummy regressor on trisynthon SMILES fingerprints 
        to predict average normalized enrichment across the provided indexes.

        Parameters:
            experiment_key (str): Key from indexes dict representing the experiment (e.g., 'protein_replicates').
            test_size (int or None): If None, use the entire dataset. Otherwise, use only the first 'test_size' rows.

        Returns:
            None. Plots the R² values for Random Forest and Dummy models and compares predictions.
        """
        
        def smiles_to_fingerprint(smiles, radius=2, n_bits=2048):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
            else:
                return np.nan

        if test_size is not None:
            data_subset = self.data.head(test_size).copy()
        else:
            data_subset = self.data.copy()

        
        print("Generating fingerprints...")
        data_subset['fingerprints'] = [smiles_to_fingerprint(smiles) for smiles in tqdm(data_subset['SMILES'])]

        # Drop rows with NaN fingerprints
        data_subset.dropna(subset=['fingerprints'], inplace=True)

        # Get experiment indexes and calculate average enrichment
        corrected_indexes = self.indexes[experiment_key]
        data_subset['average_enrichment'] = data_subset[corrected_indexes].mean(axis=1)

        
        X = np.array([list(fp) for fp in data_subset['fingerprints']])
        y = data_subset['average_enrichment']

       
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

        # Initialize Random Forest and Dummy Regressor models
        rf_model = RandomForestRegressor(random_state=42)
        dummy_model = DummyRegressor(strategy="mean")

        
        print("Training Random Forest model...")
        with tqdm(total=100, desc="Training RF Model") as pbar:
            rf_model.fit(X_train, y_train)
            pbar.update(100)
        
        print("Training Dummy model...")
        dummy_model.fit(X_train, y_train)

        
        y_pred_rf = rf_model.predict(X_test)
        y_pred_dummy = dummy_model.predict(X_test)

        
        r2_rf = r2_score(y_test, y_pred_rf)
        r2_dummy = r2_score(y_test, y_pred_dummy)

        
        plt.figure(figsize=(8, 6))
        plt.scatter(y_test, y_pred_rf, label=f"Random Forest (R² = {r2_rf:.2f})", alpha=0.7)
        plt.scatter(y_test, y_pred_dummy, label=f"Dummy Regressor (R² = {r2_dummy:.2f})", alpha=0.7, marker='x')
        plt.plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], color='black', linestyle='--')
        plt.title('Regression: Predictions vs True Values (Average Enrichment)')
        plt.xlabel('True Average Enrichment')
        plt.ylabel('Predicted Average Enrichment')
        plt.legend(loc="upper left")
        plt.show()

        
        print(f"Random Forest R²: {r2_rf:.2f}")
        print(f"Dummy Regressor R²: {r2_dummy:.2f}")



    def ml_fingerprints_to_classifier(self, experiment_key, threshold=2, test_size=None):
        """
        Train a Random Forest classifier and Dummy classifier on trisynthon SMILES fingerprints 
        to classify enriched vs. not enriched based on a given threshold.

        Parameters:
            experiment_key (str): Key from indexes dict representing the experiment (e.g., 'protein_replicates').
            threshold (float): The enrichment threshold for classification. Defaults to 2.
            test_size (int or None): Number of samples to use for testing. If None, uses the entire dataset.

        Returns:
            None. Plots the confusion matrix for both classifiers.
        """

        def smiles_to_fingerprint(smiles, radius=2, n_bits=2048):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
            else:
                return np.nan
        
       
        if test_size is not None:
            data_subset = self.data.head(test_size).copy()
        else:  # Use the entire dataset
            data_subset = self.data.copy()
        
        
        tqdm.pandas(desc="Generating fingerprints...")
        data_subset['fingerprints'] = data_subset['SMILES'].progress_apply(lambda x: smiles_to_fingerprint(x))
        
        
        data_subset.dropna(subset=['fingerprints'], inplace=True)

        
        corrected_indexes = self.indexes[experiment_key]
        data_subset['average_enrichment'] = data_subset[corrected_indexes].mean(axis=1)
        data_subset['target'] = (data_subset['average_enrichment'] > threshold).astype(int)

        
        X = np.array([list(fp) for fp in data_subset['fingerprints']])
        y = data_subset['target']

        
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

        
        rf_model = RandomForestClassifier(random_state=42)
        dummy_model = DummyClassifier(strategy="most_frequent")

        
        rf_model.fit(X_train, y_train)
        dummy_model.fit(X_train, y_train)

        
        y_pred_rf = rf_model.predict(X_test)
        y_pred_dummy = dummy_model.predict(X_test)

        
        cm_rf = confusion_matrix(y_test, y_pred_rf)
        cm_dummy = confusion_matrix(y_test, y_pred_dummy)

        
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))

        ConfusionMatrixDisplay(cm_rf, display_labels=["Not Enriched", "Enriched"]).plot(ax=ax[0])
        ax[0].set_title("Random Forest Classifier")

        ConfusionMatrixDisplay(cm_dummy, display_labels=["Not Enriched", "Enriched"]).plot(ax=ax[1])
        ax[1].set_title("Dummy Classifier")

        plt.tight_layout()
        plt.show()


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

    


