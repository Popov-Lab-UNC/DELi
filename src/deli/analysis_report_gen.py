import os
import sys
from jinja2 import Environment, FileSystemLoader
from datetime import datetime

def get_plot_files(directory):
    if os.path.exists(directory) and os.path.isdir(directory):
        return [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.png')]
    return []

def get_experiment_info(indexes, control_cols):
    experiment_info = []
    for exp_name, columns in indexes.items():
        experiment_info.append({
            'name': exp_name,
            'index': columns,
            'control_columns': control_cols.get(exp_name, [])
        })
    return experiment_info

def generate_report(base_dir, indexes, control_cols, nsc_max_dict, sd_min_dict, sampling_depth_dict):
    today_date = datetime.today().strftime('%Y%m%d')
    analysis_dir = os.path.join(base_dir, f"{today_date}_analysis")
    disynthon_dir = os.path.join(analysis_dir, "disynthon")
    trisynthon_dir = os.path.join(analysis_dir, "trisynthon")
    top_hits_dir = os.path.join(analysis_dir, "top_hits")
    ml_fingerprints_to_RF_dir = os.path.join(analysis_dir, "ml_fingerprints_to_RF")
    ml_fingerprints_to_clf_dir = os.path.join(analysis_dir, "ml_fingerprints_to_clf")

    disynthon_plots = [os.path.relpath(plot, base_dir) for plot in get_plot_files(disynthon_dir)]
    trisynthon_plots = [os.path.relpath(plot, base_dir) for plot in get_plot_files(trisynthon_dir)]
    top_hits_plots = [os.path.relpath(plot, base_dir) for plot in get_plot_files(top_hits_dir)]
    ml_fingerprints_to_RF_plots = [os.path.relpath(plot, base_dir) for plot in get_plot_files(ml_fingerprints_to_RF_dir)]
    ml_fingerprints_to_clf_plots = [os.path.relpath(plot, base_dir) for plot in get_plot_files(ml_fingerprints_to_clf_dir)]

    print(f"Disynthon plots: {disynthon_plots}")
    print(f"Trisynthon plots: {trisynthon_plots}")
    print(f"Top hits plots: {top_hits_plots}")
    print(f"ML Fingerprints to RF plots: {ml_fingerprints_to_RF_plots}")

    experiment_info = get_experiment_info(indexes, control_cols)
    
    sampling_depth_values = [
        f"{exp_name}: {round(nsc_max, 2)}, SD_min = {round(sd_min_dict.get(exp_name.replace('_NSC_max', '_SD_min'), 'N/A'), 2)}"
        for exp_name, nsc_max in nsc_max_dict.items()
    ]

    sampling_depth_values_only = [
        f"{exp_name.replace('_sampling_depth', '')}: {round(value, 2)}"
        for exp_name, value in sampling_depth_dict.items()
    ]

    template_dir = os.path.dirname(__file__)
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template("report_test_template.html")

    rendered_html = template.render(
        experiment_info=experiment_info,
        sampling_depth_values=sampling_depth_values,
        sampling_depth_only=sampling_depth_values_only,
        disynthon_plots=disynthon_plots,
        trisynthon_plots=trisynthon_plots,
        top_hits_plots=top_hits_plots,
        ml_fingerprints_to_RF_plots=ml_fingerprints_to_RF_plots,
        ml_fingerprints_to_clf_plots=ml_fingerprints_to_clf_plots,
    )

    print("Rendered HTML preview:")
    print(rendered_html[:500])

    output_file = os.path.join(base_dir, f"{today_date}_report_test_report.html")
    with open(output_file, "w") as f:
        f.write(rendered_html)

    print("Report generated successfully!")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python analysis_report_gen.py <base_directory>")
        sys.exit(1)

    base_directory = sys.argv[1]
    indexes = {}
    control_cols = {}
    nsc_max_dict = {}
    sd_min_dict = {}
    sampling_depth_dict = {}

    generate_report(base_directory, indexes, control_cols, nsc_max_dict, sd_min_dict, sampling_depth_dict)
