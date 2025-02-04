import argparse
from cube_class import DELi_Cube
import pandas as pd
import ast
import os
from datetime import datetime
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
    parser.add_argument('--data', type=str, help='Path to the data file')
    parser.add_argument('--indexes', type=str, help='Path to the file containing dictionary mapping experimental IDs to index ranges')
    parser.add_argument('--control_cols', type=str, help='Path to the file containing dictionary mapping experimental IDs to control columns')
    parser.add_argument('--lib_size', type=int, help='Library size')
    parser.add_argument('--SD_min', action='store_true', help='Print the minimum SD value')
    parser.add_argument('--NSC_values', action='store_true', help='Print the NSC values')
    parser.add_argument('--MLE', action='store_true', help='Print the maximum likelihood enrichment ratio')
    parser.add_argument('--Z_score', action='store_true', help='Print the Z-score values')
    parser.add_argument('--z_score_log_data', action='store_true', help='Print the Z-score values with log data')
    parser.add_argument('--normalized_data', action='store_true', help='Print the normalized data')
    parser.add_argument('--disynthon_data', action='store_true', help='Print the disynthon data')
    parser.add_argument('--trisynthon_overlap', action='store_true', help='Print the trisynthon overlap')
    parser.add_argument('--disynthon_overlap', action='store_true', help='Print the disynthon overlap')
    parser.add_argument('--disynthon_threshold', default=20, type=int, help='Threshold for disynthon overlap')
    parser.add_argument('--ml_fingerprints_to_RF_reg', action='store_true', help='Print the ml fingerprints to RF regression')
    parser.add_argument('--ml_fingerprints_to_RF_clf', action='store_true', help='Print the ml fingerprints to RF classification')
    parser.add_argument('--clf_thresh', type=int, help='Threshold for ml fingerprints to RF classification')
    parser.add_argument('--top_hits', type=int, help='Number of top compounds to select')
    parser.add_argument('--top_hits_metric', type=str, default='sum', help='Metric to use for selecting top compounds')
    parser.add_argument('--report', action='store_true', help='Generate a report')
    parser.add_argument('--output_dir', type=str, default='output', help='Output directory')
    return parser.parse_args()

   
def main():
    args = parse_args()

    output_dir = create_dated_output_dir(args.output_dir, "analysis")

    #need to read in data first
    df = pd.read_csv(args.data)

    indexes = read_dict_from_file(args.indexes)
    print(indexes)
    control_cols = read_dict_from_file(args.control_cols)
    print(control_cols)

    cube = DELi_Cube(df, indexes, control_cols, args.lib_size)
    #print self data
    print(cube.data.head())


    if args.SD_min:
        nsc_max_dict, sd_min_dict, sampling_depth_dict = cube.SD_min()
    if args.NSC_values:
        cube.NSC_values()
    if args.MLE:
        cube.maximum_likelihood_enrichment_ratio()
    if args.Z_score:
        cube.z_score()
    if args.z_score_log_data:
        cube.z_score_log_data()
    if args.normalized_data:
        cube.normalize()
    if args.disynthon_data:
        disynthon_data, disynth_exp_dict = cube.disynthonize()
        print(disynth_exp_dict)
    if args.trisynthon_overlap:
        trisynthon_dir = create_output_dir(os.path.join(output_dir, "trisynthon"))
        cube.trisynthon_overlap(output_dir=trisynthon_dir)
    if args.disynthon_overlap:
        disynthon_dir = create_output_dir(os.path.join(output_dir, "disynthon"))
        cube.disynthon_overlap(output_dir=disynthon_dir, disynthon_data=disynthon_data, disynth_exp_dict=disynth_exp_dict, threshold=args.disynthon_threshold)
    if args.ml_fingerprints_to_RF_reg:
        ml_fingerprints_to_RF_dir = create_output_dir(os.path.join(output_dir, "ml_fingerprints_to_RF"))
        cube.ml_fingerprints_to_RF(output_dir=ml_fingerprints_to_RF_dir)
    #clf
    if args.ml_fingerprints_to_RF_clf:
        ml_fingerprints_to_RF_clf_dir = create_output_dir(os.path.join(output_dir, "ml_fingerprints_to_clf"))
        cube.ml_fingerprints_to_classifier(output_dir=ml_fingerprints_to_RF_clf_dir, threshold=args.clf_thresh)
    if args.top_hits:
        top_hits_dir = create_output_dir(os.path.join(output_dir, "top_hits"))
        cube.top_n_compounds(args.top_hits, args.top_hits_metric, output_dir=top_hits_dir)
    if args.report:
        report.generate_report(args.output_dir, indexes, control_cols, nsc_max_dict, sd_min_dict, sampling_depth_dict)
        print("Report generation completed!")
        today_date = datetime.now().strftime("%Y%m%d")
        cube.data.to_csv(os.path.join(output_dir, f"cube_data_{today_date}.csv"), index=False)
        print(f"Cube data saved to {os.path.join(output_dir, 'cube_data.csv')}")


if __name__ == '__main__':
    main()

# Run the script with the following command:
#python analysis_parse.py --data example_data/1_19_25_DEL3_53bp1_cycle_comp_DEL003_cube.csv --indexes example_data/index.txt --control_cols example_data/control.txt --lib_size 800000 --SD_min --NSC_values --MLE --Z_score --z_score_log_data --normalized_data --disynthon_data --trisynthon_overlap --disynthon_overlap --disynthon_threshold 100 --ml_fingerprints_to_RF_reg --ml_fingerprints_to_RF_clf --clf_thresh 10 --top_hits 10 --report --output_dir test