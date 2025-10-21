import argparse
import yaml
import pandas as pd
import os
from datetime import datetime
from cube_class import DELi_Cube
import analysis_report_gen as report
import ast

def create_output_dir(output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return output_dir

def create_dated_output_dir(base_output_dir, name_suffix=""):
    output_dir = create_output_dir(base_output_dir)
    date_str = datetime.now().strftime("%Y%m%d")
    dated_output_dir = os.path.join(output_dir, f"{date_str}_{name_suffix}")
    if not os.path.exists(dated_output_dir):
        os.makedirs(dated_output_dir)
    return dated_output_dir

def parse_args():
    parser = argparse.ArgumentParser(description="Parse the necessary arguments for DEL analysis of one library")
    parser.add_argument('--config', type=str, help='Path to the YAML config file')
    return parser.parse_args()

def load_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

def read_dict_from_file(file_path):
    with open(file_path, 'r') as f:
        content = f.read()
        return ast.literal_eval(content)

def main():
    args = parse_args()
    if not args.config:
        raise ValueError("Config file path must be provided with --config")
    
    config = load_config(args.config)
    output_dir_base = config['general'].get('output_dir', 'output')
    output_dir = create_dated_output_dir(output_dir_base, "analysis")

    data_file = config['general'].get('data', '')

    if not data_file:
        raise ValueError("Data file must be specified in the YAML config.")
    df = pd.read_csv(data_file)

    if 'ID_col' in config['general']:
        id_col = config['general']['ID_col']
    else:
        id_col = df.columns[0]  # Default to the first column if ID_col is not specified
    indexes = config.get('indexes', {})
    control_cols = config.get('control_cols', {})
    raw_indexes = config.get('raw_indexes', {})

    indexes = ast.literal_eval(str(indexes))
    control_cols = ast.literal_eval(str(control_cols))
    raw_indexes = ast.literal_eval(str(raw_indexes))

    indexes_path = os.path.join(output_dir, "indexes.txt")
    control_cols_path = os.path.join(output_dir, "control_cols.txt")
    raw_indexes_path = os.path.join(output_dir, "raw_indexes.txt")
    
    with open(indexes_path, "w") as file:
        file.write(str(indexes))
    with open(control_cols_path, "w") as file:
        file.write(str(control_cols))
    with open(raw_indexes_path, "w") as file:
        file.write(str(raw_indexes))

    indexes = read_dict_from_file(indexes_path)
    control_cols = read_dict_from_file(control_cols_path)
    raw_indexes = read_dict_from_file(raw_indexes_path)
    if not indexes:
        indexes = raw_indexes

    cube = DELi_Cube(df, id_col, indexes, control_cols, int(config['general'].get('lib_size', 0)), raw_indexes)

    print(cube.data.head())
    if 'flags' in config:
        flags = config['flags']
        if 'SD_min' in flags and flags.get('SD_min', False):
            nsc_max_dict, sd_min_dict, sampling_depth_dict = cube.SD_min()
        if 'NSC_values' in flags and flags.get('NSC_values', False):
            cube.NSC_values()
        if 'MLE' in flags and flags.get('MLE', False):
            cube.maximum_likelihood_enrichment_ratio()
        if 'Z_score' in flags and flags.get('Z_score', False):
            cube.z_score()
        if 'z_score_log_data' in flags and flags.get('z_score_log_data', False):
            cube.z_score_log_data()
        if 'disynthon_data' in flags and flags.get('disynthon_data', False):
            disynthon_data, disynth_exp_dict = cube.disynthonize()
            cube.data = disynthon_data
            print(disynth_exp_dict)
        if 'polyO' in flags and flags.get('polyO', False):
            cube.PolyO()
        if 'top_disynthons' in flags and flags.get('top_disynthons', False):
            comparison_type = flags.get('top_disynthons', {}).get('comparison', 'control')
            exp_name = flags.get('top_disynthons', {}).get('exp_name', 'None')
            exp2_name = flags.get('top_disynthons', {}).get('exp2_name', 'None')
            control_name = flags.get('top_disynthons', {}).get('control_name', 'None')
            top_count = int(flags.get('top_disynthons', {}).get('top_count', 10))
            comparison_metric = flags.get('top_disynthons', {}).get('comparison_metric', 'avg')
            top_disynthons_dir = create_output_dir(os.path.join(output_dir, "top_disynthons"))
            cube.get_top_disynthons(
                disynthon_data=disynthon_data,
                exp_name1=exp_name,
                comparison_type=comparison_type,
                exp_name2=exp2_name,
                control_name=control_name,
                comparison_metric=comparison_metric,
                top_count=top_count,
                output_dir=top_disynthons_dir
            )
        if 'trisynthon_overlap' in flags and flags.get('trisynthon_overlap', False):
            trisynthon_dir = create_output_dir(os.path.join(output_dir, "trisynthon"))
            cube.trisynthon_overlap(output_dir=trisynthon_dir)
        if 'disynthon_overlap' in flags and flags.get('disynthon_overlap', False):
            disynthon_dir = create_output_dir(os.path.join(output_dir, "disynthon"))
            cube.disynthon_overlap(output_dir=disynthon_dir, disynthon_data=disynthon_data, disynth_exp_dict=disynth_exp_dict, threshold=int(flags.get('disynthon_threshold', 20)))
        if 'normalized_data' in flags and flags.get('normalized_data', False):
            cube.normalize()
        if 'simple_spotfire_version' in flags and flags.get('simple_spotfire_version', False):
            spotfire = cube.simple_spotfire_version()
            today_date = datetime.now().strftime("%Y%m%d")
            spotfire.to_csv(os.path.join(output_dir, f"spotfire_{today_date}.csv"), index=False)
        if 'ml_fingerprints_to_RF_reg' in flags and flags.get('ml_fingerprints_to_RF_reg', False):
            ml_fingerprints_to_RF_dir = create_output_dir(os.path.join(output_dir, "ml_fingerprints_to_RF"))
            cube.ml_fingerprints_to_RF(output_dir=ml_fingerprints_to_RF_dir)
        if 'ml_fingerprints_to_RF_clf' in flags and flags.get('ml_fingerprints_to_RF_clf', False):
            ml_fingerprints_to_RF_clf_dir = create_output_dir(os.path.join(output_dir, "ml_fingerprints_to_clf"))
            cube.ml_fingerprints_to_classifier(output_dir=ml_fingerprints_to_RF_clf_dir, threshold=int(flags.get('clf_thresh', 10)))
        if 'gnn_classifier' in flags and flags.get('gnn_classifier', False):
            gnn_dir = create_output_dir(os.path.join(output_dir, "gnn"))
            cube.gnn_classifier(output_dir=gnn_dir, threshold=int(flags.get('gnn_threshold', 10)), arch=flags.get('gnn_arch', 'GAT'))
        if 'top_hits' in flags:
            top_hits_dir = create_output_dir(os.path.join(output_dir, "top_hits"))
            cube.top_n_compounds(int(flags['top_hits']), flags.get('top_hits_metric', 'sum'), output_dir=top_hits_dir)
        if 'monosynthon_chemical_space' in flags and flags.get('monosynthon_chemical_space', False):
            cube.monosynthon_chemical_space(output_dir=output_dir)
        if 'report' in flags and flags.get('report', False):
            nsc_max_dict = nsc_max_dict if 'SD_min' in flags and flags.get('SD_min', False) else None
            sd_min_dict = sd_min_dict if 'SD_min' in flags and flags.get('SD_min', False) else None
            sampling_depth_dict = sampling_depth_dict if 'SD_min' in flags and flags.get('SD_min', False) else None
            report.generate_report(output_dir_base, indexes, control_cols, nsc_max_dict, sd_min_dict, sampling_depth_dict)
            print("Report generation completed!")
            today_date = datetime.now().strftime("%Y%m%d")
            cube.data.to_csv(os.path.join(output_dir_base, f"cube_data_{today_date}.csv"), index=False)
            print(f"Cube data saved to {os.path.join(output_dir, 'cube_data.csv')}")

if __name__ == '__main__':
    main()
