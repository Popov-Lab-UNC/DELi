import math

import numpy as np
import pandas as pd
import scipy.stats
from tqdm import tqdm


class PolyO:
    def __init__(self, data, indexes, raw_indexes, lib_size, feature_mode="disynthon"):
        self.data = data
        self.indexes = indexes
        self.raw_indexes = raw_indexes
        self.lib_size = lib_size
        # self.collection_size = collection_size if collection_size is not None else lib_size
        self.feature_mode = feature_mode

        # if pool_reads is not None:
        #     self.pool_reads = pool_reads
        # else:
        #     for exp_name, columns in self.raw_indexes.items():
        #          raw_row_sums= self.data[columns].sum(axis=1)
        #          self.pool_reads = raw_row_sums.sum(axis=0)

    def calculate_polyOraw(self):
        """Compute PolyOraw scores for each disynthon, per experiment."""
        synthons = (
            ["AB", "BC", "AC"] if self.feature_mode == "disynthon" else ["ID_A", "ID_B", "ID_C"]
        )

        print("Calculating total reads and sampling depth per experiment...")
        total_reads_per_exp = {}
        total_sampling_depth_per_exp = {}

        # total reads and sampling depth for each experiment
        for exp_name, columns in tqdm(self.raw_indexes.items(), desc="Experiments"):
            raw_row_sums = self.data[columns].sum(axis=1)
            total_reads = raw_row_sums.sum()
            total_reads_per_exp[exp_name] = total_reads
            total_sampling_depth_per_exp[exp_name] = total_reads / self.lib_size
            print(f"Total reads for {exp_name}: {total_reads}")
            print(f"Sampling depth for {exp_name}: {total_sampling_depth_per_exp[exp_name]}")

        print("Computing Poisson probabilities per row, per experiment...")
        p_compounds = {}  # Poisson probabilities for each experiment

        for exp_name, columns in tqdm(self.indexes.items(), desc="Experiments (Poisson)"):
            sum_columns = self.data[columns].sum(axis=1)
            poisson_pmf = sum_columns.apply(
                lambda x: scipy.stats.poisson.pmf(x, total_sampling_depth_per_exp[exp_name])
            )
            p_compounds[exp_name] = poisson_pmf
            # Store p_compounds in the DataFrame as well
            self.data[f"p_compound_{exp_name}"] = poisson_pmf
        print("Calculating product for all matching features...")

        # Calculate prod of p_compound for matching features in each exp
        for exp_name in tqdm(self.indexes.keys(), desc="Calculating Products per Disynthon"):
            total_reads = sum(
                self.data[columns].sum().sum() for columns in self.raw_indexes.values()
            )
            S_bar = total_reads / self.lib_size
            self.data["S_bar"] = S_bar

            for synthon in synthons:
                # Group by synthon and calculate the product of p_compound for each experiment
                prod_df = (
                    self.data.groupby(synthon)[f"p_compound_{exp_name}"]
                    .apply(np.prod)
                    .reset_index(name=f"{synthon}_Prod_{exp_name}")
                )
                self.data = pd.merge(self.data, prod_df, on=synthon, how="left")
                self.data[f"PolyO_raw_{synthon}_{exp_name}"] = -np.log10(
                    self.data[f"{synthon}_Prod_{exp_name}"]
                )

                # Debugging: Check the merged products and PolyO_raw values
                # print(f"Calculated product and PolyO_raw for {synthon} and {exp_name}:\n{self.data[[synthon, f'{synthon}_Prod_{exp_name}', f'PolyO_raw_{synthon}_{exp_name}']].head()}")

        # Debugging: Check the final DataFrame with products and PolyO_raw values
        # print(f"Updated DataFrame with products and PolyO_raw values:\n{self.data.head()}")

        return self.data

    def find_Ccpd(self, threshold=0.01):
        """_summary_

        Parameters
        ----------
        threshold : float, optional
            _description_, by default 0.01
        """
        features = (
            ["AB", "BC", "AC"] if self.feature_mode == "disynthon" else ["ID_A", "ID_B", "ID_C"]
        )
        Ccpd_values = {}

        for feature in features:
            same_feat_comps = self.data.groupby(feature)
            parallel_features = len(same_feat_comps)
            cmps_in_DEL = self.lib_size
            c_bar = cmps_in_DEL / parallel_features

            k = math.ceil(c_bar) + 1
            while scipy.stats.poisson.pmf(k, c_bar) > threshold:
                k += 1

            Ccpd_values[feature] = k

        for feature, cutoff in Ccpd_values.items():
            self.data[f"{feature}_cpd_cutoff"] = cutoff

        return Ccpd_values

    def calculate_Cread(self):
        """Calculate the Cread cut-off values for each feature."""
        features = (
            ["AB", "BC", "AC"] if self.feature_mode == "disynthon" else ["ID_A", "ID_B", "ID_C"]
        )
        Cread_values = {}

        for feature in features:
            total_DEL_members = self.lib_size
            total_reads = 0
            for columns in self.raw_indexes.values():
                total_reads += self.data[columns].sum().sum()

            S_bar = total_reads / total_DEL_members
            k = math.ceil(S_bar) + 1
            while scipy.stats.poisson.pmf(k, S_bar) > 0.01:
                k += 1
            Cread_values[feature] = k

        for feature, cutoff in Cread_values.items():
            self.data[f"{feature}_read_cutoff"] = cutoff

        return Cread_values

    def poly_o_base(self):
        """Compute PolyObase scores using Equation 9."""
        synthons = (
            ["AB", "BC", "AC"] if self.feature_mode == "disynthon" else ["ID_A", "ID_B", "ID_C"]
        )

        Ccpd_values = self.find_Ccpd()
        Cread_values = self.calculate_Cread()

        total_reads = sum(self.data[columns].sum().sum() for columns in self.raw_indexes.values())
        S_bar = total_reads / self.lib_size

        print("Calculating PolyObase scores using Eq. 9...")

        for synthon in synthons:
            Ccpd = Ccpd_values[synthon]
            Cread = Cread_values[synthon]

            Cread_factorial = math.factorial(Cread)
            self.data[f"PolyObase_{synthon}"] = -Ccpd * np.log10(
                (S_bar**Cread * np.exp(-S_bar)) / Cread_factorial
            )

            # Debugging: Check the computed PolyObase scores
            # print(f"Computed PolyObase for {synthon}:\n{self.data[[synthon, f'PolyObase_{synthon}']].head()}")

        return self.data

    def calculate_polyO_score(self):
        """Calculate PolyO score for each disynthon, per experiment."""
        synthons = (
            ["AB", "BC", "AC"] if self.feature_mode == "disynthon" else ["ID_A", "ID_B", "ID_C"]
        )

        print("Calculating PolyO score for each disynthon...")

        for synthon in synthons:
            for exp_name in self.indexes.keys():
                polyOraw_col = f"PolyO_raw_{synthon}_{exp_name}"
                polyObase_col = f"PolyObase_{synthon}"

                # Check if PolyOraw and PolyObase columns exist before calculating the PolyO score
                if polyOraw_col not in self.data.columns:
                    print(
                        f"Warning: {polyOraw_col} not found in DataFrame. Skipping calculation for this column."
                    )
                    continue
                if polyObase_col not in self.data.columns:
                    print(
                        f"Warning: {polyObase_col} not found in DataFrame. Skipping calculation for this column."
                    )
                    continue
                self.data[f"{synthon}_{exp_name}_PolyO_score"] = (
                    self.data[polyOraw_col] / self.data[polyObase_col]
                )

        #         # Debugging: Print the first few rows of the calculated PolyO scores
        #         print(f"Calculated PolyO score for {synthon} and {exp_name}:\n{self.data[[synthon, f'PolyO_score_{synthon}_{exp_name}']].head()}")

        # # Debugging: Check the final DataFrame with PolyO scores
        # print(f"Updated DataFrame with PolyO scores:\n{self.data.head()}")

        return self.data


# #test the class
# df = pd.read_csv("/Users/brandonnovy/Desktop/deli_repo/DELi/src/deli/DEL3/cube_data_20250318.csv")
# indexes = {
#     "Cycle1": ["corrected_PTMODD3index4", "corrected_PTMODD3index5", "corrected_PTMODD3index6"],
#     "Cycle2": ["corrected_PTMODD3index1", "corrected_PTMODD3index2", "corrected_PTMODD3index3"]
# }
# raw_indexes = {
#     "Cycle1": ["raw_PTMODD3index4", "raw_PTMODD3index5", "raw_PTMODD3index6"],
#     "Cycle2": ["raw_PTMODD3index1", "raw_PTMODD3index2", "raw_PTMODD3index3"]
# }
# lib_size = 58000
# collection_size = 58000

# polyo = PolyO(df, indexes, raw_indexes, lib_size, collection_size, feature_mode="disynthon")
# polyo.calculate_polyOraw()
# #save the data as test
# polyo.find_Ccpd()
# # polyo.find_Ccpd()
# polyo.calculate_Cread()
# polyo.poly_o_base()
# polyo.calculate_polyO_score()
# polyo.data.to_csv("test.csv")
