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


def get_plot_files(directory, file_extensions=[".svg"]):
    """Retrieve plots from a given directory (SVG only for report embedding)."""
    if os.path.exists(directory) and os.path.isdir(directory):
        plots = []
        files = sorted([f for f in os.listdir(directory) if any(f.endswith(ext) for ext in file_extensions)])
        for f in files:
            file_path = os.path.join(directory, f)
            with open(file_path, "r", encoding="utf-8") as svg_file:
                plots.append({
                    "name": f,
                    "data": svg_file.read(),
                    "type": "svg"
                })
        return plots
    return []


def _plot_order_key(filename: str) -> tuple:
    synthon_order = {"AB": 0, "AC": 1, "BC": 2}
    parts = filename.split("_")
    synthon = next((p for p in parts if p in synthon_order), "ZZ")
    # keep comparison cohorts together while preserving AB/AC/BC ordering
    if filename.startswith("top_"):
        cohort = filename.replace("top_AB_", "").replace("top_AC_", "").replace("top_BC_", "")
    else:
        cohort = filename.replace("AB_", "").replace("AC_", "").replace("BC_", "")
    return (cohort, synthon_order.get(synthon, 99), filename)


def get_html_files(directory, name_prefixes=None, include_fn=None):
    """Retrieve HTML files as text from a given directory."""
    if os.path.exists(directory) and os.path.isdir(directory):
        html_files = []
        files = sorted([f for f in os.listdir(directory) if f.endswith(".html")])
        if name_prefixes:
            files = [f for f in files if any(f.startswith(prefix) for prefix in name_prefixes)]
        if include_fn:
            files = [f for f in files if include_fn(f)]
        for f in files:
            with open(os.path.join(directory, f), 'r', encoding='utf-8') as file:
                html_files.append({
                    "name": f,
                    "content": file.read()
                })
        return html_files
    return []


def get_experiment_info(indexes, control_cols):
    experiment_info = []
    for exp_name, columns in indexes.items():
        # Handle case where control_cols is None or empty
        control_columns = control_cols.get(exp_name, []) if control_cols else []
        experiment_info.append(
            {"name": exp_name, "index": columns, "control_columns": control_columns}
        )
    return experiment_info


def generate_report(
    base_dir,
    indexes,
    control_cols,
    nsc_max_dict=None,
    sd_min_dict=None,
    sampling_depth_dict=None,
    top_disynthon_specs=None,
    top_delta_specs=None,
    disynthon_thresholds=None,
):
    today_date = datetime.today().strftime("%Y%m%d")
    analysis_dir = os.path.join(base_dir, f"{today_date}_analysis")

    # Define directories for different plot categories
    plot_dirs = {
        "disynthon_plots": os.path.join(analysis_dir, "disynthon"),
        "trisynthon_plots": os.path.join(analysis_dir, "trisynthon"),
        "top_disynthons_plots": os.path.join(analysis_dir, "top_disynthons"),
        "top_hits_plots": os.path.join(analysis_dir, "top_hits"),
        "top_delta_plots": os.path.join(analysis_dir, "top_deltas"),
        "ml_fingerprints_to_RF_plots": os.path.join(analysis_dir, "ml_fingerprints_to_RF"),
        "ml_fingerprints_to_clf_plots": os.path.join(analysis_dir, "ml_fingerprints_to_clf"),
        "gnn_classifier_plots": os.path.join(analysis_dir, "gnn"),
        "monosynthon_chemical_space_plots": os.path.join(analysis_dir, "clusters"),
    }

    # Convert plots to base64
    plot_data = {}
    for key, path in plot_dirs.items():
        if key == "monosynthon_chemical_space_plots":
            html_files = get_html_files(
                path,
                ["monosynthon_"],
                include_fn=lambda name: len(name.split("_")) >= 5 and "comparison" not in name,
            )
            # Add a type field to distinguish in template
            for item in html_files:
                item["type"] = "html"
            plot_data[key] = html_files
        else:
            # Handle PNG and SVG files for other plots
            plots = get_plot_files(path, [".svg"])
            if key in {"disynthon_plots", "top_disynthons_plots"}:
                plots = sorted(plots, key=lambda x: _plot_order_key(x["name"]))
            plot_data[key] = plots

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

    # Use the source template directly to ensure we get the latest version
    template_path = os.path.join(
        os.path.dirname(__file__), "..", "templates", "analysis_report.html"
    )
    template = jinja2.Template(open(template_path).read())

    rendered_html = template.render(
        experiment_info=experiment_info,
        sampling_depth_values=sampling_depth_values,
        sampling_depth_only=sampling_depth_values_only,
        top_disynthon_specs=top_disynthon_specs or [],
        top_delta_specs=top_delta_specs or [],
        disynthon_thresholds=disynthon_thresholds,
        now_date=datetime.now().strftime("%Y-%m-%d"),
        now_year=datetime.now().year,
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
