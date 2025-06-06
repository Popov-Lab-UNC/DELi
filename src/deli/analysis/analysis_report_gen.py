import base64
import importlib.resources as resources
import os
import sys
from datetime import datetime

import jinja2


def encode_image_to_base64(image_path):
    """Convert image to base64-encoded string."""
    with open(image_path, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode("utf-8")


def get_plot_files(directory):
    """Retrieve base64-encoded plots from a given directory."""
    if os.path.exists(directory) and os.path.isdir(directory):
        return [
            encode_image_to_base64(os.path.join(directory, f))
            for f in os.listdir(directory)
            if f.endswith(".png")
        ]
    return []


def get_experiment_info(indexes, control_cols):
    experiment_info = []
    for exp_name, columns in indexes.items():
        experiment_info.append(
            {"name": exp_name, "index": columns, "control_columns": control_cols.get(exp_name, [])}
        )
    return experiment_info


def generate_report(
    base_dir, indexes, control_cols, nsc_max_dict=None, sd_min_dict=None, sampling_depth_dict=None
):
    today_date = datetime.today().strftime("%Y%m%d")
    analysis_dir = os.path.join(base_dir, f"{today_date}_analysis")

    # Define directories for different plot categories
    plot_dirs = {
        "disynthon_plots": os.path.join(analysis_dir, "disynthon"),
        "trisynthon_plots": os.path.join(analysis_dir, "trisynthon"),
        "top_disynthons_plots": os.path.join(analysis_dir, "top_disynthons"),
        "top_hits_plots": os.path.join(analysis_dir, "top_hits"),
        "ml_fingerprints_to_RF_plots": os.path.join(analysis_dir, "ml_fingerprints_to_RF"),
        "ml_fingerprints_to_clf_plots": os.path.join(analysis_dir, "ml_fingerprints_to_clf"),
        "gnn_classifier_plots": os.path.join(analysis_dir, "gnn"),
    }

    # Convert plots to base64
    plot_data = {key: get_plot_files(path) for key, path in plot_dirs.items()}

    experiment_info = get_experiment_info(indexes, control_cols)

    nsc_max_dict = nsc_max_dict or {}
    sd_min_dict = sd_min_dict or {}
    sampling_depth_dict = sampling_depth_dict or {}

    sampling_depth_values = [
        f"{exp_name}: {round(nsc_max, 2)}, SD_min = {round(sd_min_dict.get(exp_name.replace('_NSC_max', '_SD_min'), 'N/A'), 2)}"
        for exp_name, nsc_max in nsc_max_dict.items()
    ]

    sampling_depth_values_only = [
        f"{exp_name.replace('_sampling_depth', '')}: {round(value, 2)}"
        for exp_name, value in sampling_depth_dict.items()
    ]

    with resources.path("deli.templates", "analysis_report.html") as template_path:
        template = jinja2.Template(open(template_path).read())

        rendered_html = template.render(
            experiment_info=experiment_info,
            sampling_depth_values=sampling_depth_values,
            sampling_depth_only=sampling_depth_values_only,
            **plot_data,  # Inject base64 plot data
        )

        # Save the report as an HTML file
        output_file = os.path.join(base_dir, f"{today_date}_analysis_report.html")
        with open(output_file, "w") as f:
            f.write(rendered_html)

    print(f"Report generated successfully: {output_file}")


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

    generate_report(
        base_directory, indexes, control_cols, nsc_max_dict, sd_min_dict, sampling_depth_dict
    )
