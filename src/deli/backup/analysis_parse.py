import argparse
import configparser
import pandas as pd
import ast
import os
from datetime import datetime
from cube_class import DELi_Cube
import analysis_report_gen as report

def create_output_dir(output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return output_dir

def save_plot(plot, output_dir, filename):
    plot_path = os.path.join(output_dir, filename)
    plot.savefig(plot_path)
    print(f"Plot saved to {plot_path}")

def generate_report(report_content, output_dir, filename):
    report_path = os.path.join(output_dir, filename)
    with open(report_path, 'w') as report_file:
        report_file.write(report_content)
    print(f"Report saved to {report_path}")

def read_dict_from_file(file_path):
    with open(file_path, 'r') as f:
        content = f.read()
    return ast.literal_eval(content)

def create_dated_output_dir(base_output_dir, name_suffix=""):
    output_dir = create_output_dir(base_output_dir)
    date_str = datetime.now().strftime("%Y%m%d")
    dated_output_dir = os.path.join(output_dir, f"{date_str}_{name_suffix}")
    if not os.path.exists(dated_output_dir):
        os.makedirs(dated_output_dir)
    return dated_output_dir

def parse_args():
    parser = argparse.ArgumentParser(description="Parse the necessary arguments for DEL analysis of one library")
    parser.add_argument('--config', type=str, help='Path to the config.ini file')
    return parser.parse_args()

def load_config(config_path):
    config = configparser.ConfigParser()
    config.read(config_path)
    return config

def main():
    args = parse_args()

    # Load config from config.ini if provided
    config = None
    if args.config:
        config = load_config(args.config)

    # If no config is provided, fall back to argparse input
    output_dir_base = config['general'].get('output_dir', 'output') if config else 'output'
    output_dir = create_dated_output_dir(output_dir_base, "analysis")

    # Read data
    data_file = config['general'].get('data', '')
    if not data_file:
        raise ValueError("Data file must be specified.")
    df = pd.read_csv(data_file)

    indexes = read_dict_from_file(config['general'].get('indexes', ''))
    control_cols = read_dict_from_file(config['general'].get('control_cols', ''))
    raw_indexes = read_dict_from_file(config['general'].get('indexes', ''))

    cube = DELi_Cube(df, indexes, control_cols, int(config['general'].get('lib_size', 0)), raw_indexes)
    print(cube.data.head())

    # Process flags from the config file
    if config and config.has_section('flags'):
        if config['flags'].getboolean('SD_min'):
            nsc_max_dict, sd_min_dict, sampling_depth_dict = cube.SD_min()
        if config['flags'].getboolean('NSC_values'):
            cube.NSC_values()
        if config['flags'].getboolean('MLE'):
            cube.maximum_likelihood_enrichment_ratio()
        if config['flags'].getboolean('Z_score'):
            cube.z_score()
        if config['flags'].getboolean('z_score_log_data'):
            cube.z_score_log_data()
        if config['flags'].getboolean('normalized_data'):
            cube.normalize()
        if config['flags'].getboolean('disynthon_data'):
            disynthon_data, disynth_exp_dict = cube.disynthonize()
            print(disynth_exp_dict)
        if config['flags'].getboolean('trisynthon_overlap'):
            trisynthon_dir = create_output_dir(os.path.join(output_dir, "trisynthon"))
            cube.trisynthon_overlap(output_dir=trisynthon_dir)
        if config['flags'].getboolean('disynthon_overlap'):
            disynthon_dir = create_output_dir(os.path.join(output_dir, "disynthon"))
            cube.disynthon_overlap(output_dir=disynthon_dir, disynthon_data=disynthon_data, disynth_exp_dict=disynth_exp_dict, threshold=int(config['flags'].get('disynthon_threshold', 20)))
        if config['flags'].getboolean('simple_spotfire_version'):
            spotfire = cube.simple_spotfire_version()
            today_date = datetime.now().strftime("%Y%m%d")
            spotfire.to_csv(os.path.join(output_dir, f"spotfire_{today_date}.csv"), index=False)
        if config['flags'].getboolean('ml_fingerprints_to_RF_reg'):
            ml_fingerprints_to_RF_dir = create_output_dir(os.path.join(output_dir, "ml_fingerprints_to_RF"))
            cube.ml_fingerprints_to_RF(output_dir=ml_fingerprints_to_RF_dir)
        if config['flags'].getboolean('ml_fingerprints_to_RF_clf'):
            ml_fingerprints_to_RF_clf_dir = create_output_dir(os.path.join(output_dir, "ml_fingerprints_to_clf"))
            cube.ml_fingerprints_to_classifier(output_dir=ml_fingerprints_to_RF_clf_dir, threshold=int(config['flags'].get('clf_thresh', 10)))
        if config['flags'].get('top_hits', None):
            top_hits_dir = create_output_dir(os.path.join(output_dir, "top_hits"))
            cube.top_n_compounds(int(config['flags']['top_hits']), config['flags'].get('top_hits_metric', 'sum'), output_dir=top_hits_dir)
        if config['flags'].getboolean('report'):
            report.generate_report(output_dir_base, indexes, control_cols, nsc_max_dict, sd_min_dict, sampling_depth_dict)
            print("Report generation completed!")
            today_date = datetime.now().strftime("%Y%m%d")
            cube.data.to_csv(os.path.join(output_dir_base, f"cube_data_{today_date}.csv"), index=False)
            print(f"Cube data saved to {os.path.join(output_dir, 'cube_data.csv')}")

if __name__ == '__main__':
    main()
