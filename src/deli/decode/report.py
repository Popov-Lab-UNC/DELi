"""handles making html reports for deli decoding"""

import importlib.resources as resources
import os
from datetime import datetime

import jinja2
import plotly.graph_objects as go

from deli.dels import Selection

from .decoder import DecodeStatistics


def _generate_calling_pie_chart(
    decode_stats: DecodeStatistics,
):
    """Generate the main calling results pie chart html"""
    _all_colors = [
        "orange",
        "yellow",
        "firebrick",
        "pink",
        "orchid",
        "crimson",
        "lightpink",
        "royalblue",
    ]
    _all_labels = [
        "Read to Small",
        "Read to Large",
        "Failed Library Call",
        "Library Call Malformed",
        "Failed Barcode Alignment",
        "Failed Building Block Call",
        "UMI Match to Small",
        "Decoded",
    ]
    _all_values: list[int] = [
        decode_stats.num_failed_too_short,
        decode_stats.num_failed_too_long,
        decode_stats.num_failed_library_call,
        decode_stats.num_failed_library_match_too_short,
        decode_stats.num_failed_building_block_call,
        decode_stats.num_failed_alignment,
        decode_stats.num_failed_umi_match_too_short,
        decode_stats.num_seqs_decoded,
    ]

    picked_labels: list[str] = []
    picked_colors: list[str] = []
    picked_values: list[int] = []

    for _label, _color, _val in zip(_all_labels, _all_colors, _all_values):
        if _val > 0:
            picked_values.append(_val)
            picked_labels.append(_label)
            picked_colors.append(_color)

    fig = go.Figure(
        data=[
            go.Pie(
                labels=picked_labels,
                values=picked_values,
                hole=0.3,
                pull=[0] * (len(picked_labels) - 1) + [0.2],
                hoverinfo="skip",
                textinfo="percent+value+label",
                textfont_size=18,
                marker=dict(
                    colors=picked_colors,
                    line=dict(color="#000000", width=2),
                ),
            )
        ]
    )

    fig.update_layout(
        margin=dict(t=0, b=0, l=0, r=0),
        autosize=False,
        showlegend=False,
        width=850,
        height=650,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def _generate_seq_length_hist(decode_stats: DecodeStatistics):
    """Generate the html for the seq length histogram"""
    x_vals = list(
        range(min(decode_stats.seq_lengths.keys()), max(decode_stats.seq_lengths.keys()) + 1)
    )
    y_vals = [decode_stats.seq_lengths.get(key, 0) for key in x_vals]
    fig = go.Figure(data=[go.Bar(x=x_vals, y=y_vals)])
    fig.update_layout(
        xaxis_title="Sequence Read Length",
        yaxis_title="Count",
        margin=dict(t=0, b=0, l=0, r=0),
        autosize=False,
        showlegend=False,
        width=1000,
        height=300,
    )
    fig.update_xaxes(range=[0, 1000])
    return fig.to_html(full_html=False, include_plotlyjs=False)


def _generate_lib_decode_pie_chart(decode_stats: DecodeStatistics):
    """Generate the html for a lib decode pie chart"""
    _labels = []
    _values = []
    for idx, count in decode_stats.num_seqs_decoded_per_lib.items():
        if count == 0:
            continue
        _labels.append(idx)
        _values.append(count)

    fig = go.Figure(
        data=[
            go.Pie(
                labels=_labels,
                values=_values,
                hole=0.3,
                hoverinfo="skip",
                textinfo="percent+value+label",
                textfont_size=18,
                marker=dict(line=dict(color="#000000", width=2)),
            )
        ]
    )

    fig.update_layout(
        margin=dict(t=0, b=0, l=0, r=0),
        autosize=False,
        showlegend=False,
        width=500,
        height=400,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def _generate_lib_degen_pie_chart(decode_stats: DecodeStatistics):
    """Generate the html for a lib degen pie chart"""
    _labels = []
    _values = []
    for idx, count in decode_stats.num_seqs_degen_per_lib.items():
        if count == 0:
            continue
        _labels.append(idx)
        _values.append(count)

    fig = go.Figure(
        data=[
            go.Pie(
                labels=_labels,
                values=_values,
                hole=0.3,
                hoverinfo="skip",
                textinfo="percent+value+label",
                textfont_size=18,
                marker=dict(line=dict(color="#000000", width=2)),
            )
        ]
    )

    fig.update_layout(
        margin=dict(t=0, b=0, l=0, r=0),
        autosize=False,
        showlegend=False,
        width=500,
        height=400,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def build_decoding_report(
    selection: Selection,
    stats: DecodeStatistics,
    out_path: str | os.PathLike,
):
    """
    Generates a deli decoding html report from report stats

    The name of the report will be "{selection_id}_decode_report.html"

    Notes
    -----
    Will render the report using the decode_report.html resource

    Parameters
    ----------
    selection: SequencedSelection
        the selection to build the report for
    stats: DecodeStatistics
        the decode run statistics to build the report for
    out_path: Union[str, os.PathLike]
        the path to save the report to
    """
    _seq_count_data: dict = {
        **{
            "Total": (
                "{:,}".format(selection.library_collection.collection_size),
                "{:,}".format(stats.num_seqs_decoded),
                "{:,}".format(stats.num_seqs_degen),
            )
        },
        **{
            lib.library_id: (
                "{:,}".format(lib.library_size),
                "{:,}".format(stats.num_seqs_decoded_per_lib.get(lib.library_id, 0)),
                "{:,}".format(stats.num_seqs_degen_per_lib.get(lib.library_id, 0)),
            )
            for lib in selection.library_collection.libraries
        },
    }

    jinja_data = {
        "timestamp": datetime.now().strftime("%b %d, %Y %H:%M"),
        "selection": selection.selection_id,
        "run_date": selection.get_run_date_as_str(),
        "target": selection.selection_condition.target_id
        if selection.selection_condition.target_id
        else "NA",
        "selection_cond": selection.selection_condition.selection_condition
        if selection.selection_condition.selection_condition
        else "NA",
        "additional_info": selection.selection_condition.additional_info
        if selection.selection_condition.additional_info
        else "NA",
        "libs": ", ".join(
            sorted([lib.library_id for lib in selection.library_collection.libraries])
        ),
        "num_reads": stats.num_seqs_read,
        "num_decoded": stats.num_seqs_decoded,
        "num_degen": stats.num_seqs_degen,
        "seq_counts": _seq_count_data,
        "read_hist": _generate_seq_length_hist(stats),
        "decoding_pie": _generate_calling_pie_chart(stats),
        "decode_pie": _generate_lib_decode_pie_chart(stats),
        "degen_pie": _generate_lib_degen_pie_chart(stats),
    }

    # write the report as a rendered jinja2 template
    # the template is stored as a resource in the package under templates
    with resources.path("deli.templates", "decode_report.html") as template_path:
        template = jinja2.Template(open(template_path).read())
        rendered_content = template.render(jinja_data)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(rendered_content)
