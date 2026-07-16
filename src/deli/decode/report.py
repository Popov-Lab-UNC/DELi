"""handles making html reports for deli decoding"""

import importlib.resources as resources
import os
from datetime import datetime
from typing import Optional

import jinja2
import plotly.graph_objects as go

from deli.selection import Selection, SequencedSelection

from .stats import DecodeStatistics


def _generate_calling_pie_chart(
    decode_stats: DecodeStatistics,
):
    """Generate the main calling results pie chart HTML"""
    _all_colors = [
        "#F5A623",
        "#F8E71C",
        "#F87171",
        "#EF4444",
        "#DC2626",
        "#F97316",
        "#3B82F6",
    ]
    _all_labels = [
        "Read to Small",
        "Read to Large",
        "Failed to Align Read",
        "Failed Library Call",
        "Failed Building Block Call",
        "Failed to Assign UMI",
        "Decoded",
    ]
    _all_values: list[int] = [
        decode_stats.num_failed_too_short,
        decode_stats.num_failed_too_long,
        decode_stats.num_failed_alignment,
        decode_stats.num_failed_library_call,
        decode_stats.num_failed_building_block_call,
        decode_stats.num_failed_umi,
        decode_stats.num_seqs_decoded,
    ]

    picked_labels: list[str] = []
    picked_colors: list[str] = []
    picked_values: list[int] = []

    for _label, _color, _val in zip(_all_labels, _all_colors, _all_values, strict=False):
        if _val > 0:
            picked_values.append(_val)
            picked_labels.append(_label)
            picked_colors.append(_color)

    # Sort descending by value
    sorted_data = sorted(
        zip(picked_labels, picked_colors, picked_values, strict=False),
        key=lambda x: x[2],
        reverse=True
    )
    if sorted_data:
        picked_labels, picked_colors, picked_values = zip(*sorted_data)
        picked_labels = list(picked_labels)
        picked_colors = list(picked_colors)
        picked_values = list(picked_values)
    else:
        picked_labels, picked_colors, picked_values = [], [], []

    # Pull "Decoded" specifically
    pull = [0] * len(picked_labels)
    try:
        decoded_idx = picked_labels.index("Decoded")
        pull[decoded_idx] = 0.2
    except ValueError:
        pass

    fig = go.Figure(
        data=[
            go.Pie(
                labels=picked_labels,
                values=picked_values,
                hole=0.3,
                pull=pull,
                hoverinfo="label+value+percent",
                textinfo="percent",
                textfont_size=16,
                marker=dict(
                    colors=picked_colors,
                    line=dict(color="#ffffff", width=2),
                ),
            )
        ]
    )

    fig.update_layout(
        margin=dict(t=10, b=10, l=10, r=10),
        autosize=True,
        showlegend=False,
        height=380,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False, config={'responsive': True})


def _generate_calling_bar_chart(
    decode_stats: DecodeStatistics,
):
    """Generate the main calling results bar chart HTML"""
    _all_colors = [
        "#F5A623",
        "#F8E71C",
        "#F87171",
        "#EF4444",
        "#DC2626",
        "#F97316",
        "#3B82F6",
    ]
    _all_labels = [
        "Read to Small",
        "Read to Large",
        "Failed to Align Read",
        "Failed Library Call",
        "Failed Building Block Call",
        "Failed to Assign UMI",
        "Decoded",
    ]
    _all_values: list[int] = [
        decode_stats.num_failed_too_short,
        decode_stats.num_failed_too_long,
        decode_stats.num_failed_alignment,
        decode_stats.num_failed_library_call,
        decode_stats.num_failed_building_block_call,
        decode_stats.num_failed_umi,
        decode_stats.num_seqs_decoded,
    ]

    picked_labels: list[str] = []
    picked_colors: list[str] = []
    picked_values: list[int] = []

    for _label, _color, _val in zip(_all_labels, _all_colors, _all_values, strict=False):
        if _val > 0:
            picked_values.append(_val)
            picked_labels.append(_label)
            picked_colors.append(_color)

    # Sort descending by value (longest bar at the top)
    sorted_data = sorted(
        zip(picked_labels, picked_colors, picked_values, strict=False),
        key=lambda x: x[2],
        reverse=True
    )
    if sorted_data:
        picked_labels, picked_colors, picked_values = zip(*sorted_data)
        picked_labels = list(picked_labels)
        picked_colors = list(picked_colors)
        picked_values = list(picked_values)
    else:
        picked_labels, picked_colors, picked_values = [], [], []

    fig = go.Figure(
        data=[
            go.Bar(
                y=picked_labels,
                x=picked_values,
                orientation="h",
                marker=dict(
                    color=picked_colors,
                    line=dict(color="#ffffff", width=1.5),
                ),
                text=picked_values,
                texttemplate="%{x:,}",
                textposition="outside",
                cliponaxis=False,
                hoverinfo="x+y",
            )
        ]
    )

    fig.update_layout(
        margin=dict(t=20, b=40, l=160, r=40),
        autosize=True,
        showlegend=False,
        height=400,
        xaxis=dict(
            title=dict(text="Sequence Count", font=dict(size=13, family="Inter, sans-serif")),
            tickfont=dict(size=10, family="Inter, sans-serif")
        ),
        yaxis=dict(
            title=dict(text="Calling Category", font=dict(size=13, family="Inter, sans-serif")),
            tickfont=dict(size=10, family="Inter, sans-serif"),
            autorange="reversed"
        ),
    )

    return fig.to_html(full_html=False, include_plotlyjs=False, config={'responsive': True})


def _default_to_NA(value: Optional[str]) -> str:
    """Convert None or empty strings to 'NA'"""
    if value is None or value.strip() == "":
        return "NA"
    return value


def build_decoding_report(
    stats: DecodeStatistics,
    out_path: Optional[os.PathLike] = None,
    selection: Optional[Selection] = None,
):
    """
    Generates a deli decoding HTML report from report stats

    The name of the report will be "{selection_id}_decode_report.html"

    Notes
    -----
    Will render the report using the decode_report.html resource

    Parameters
    ----------
    stats: DecodeStatistics
        the decode run statistics to build the report for
    out_path: os.PathLike
        the path to save the report to
    selection: Selection, optional
        the selection to build the report for.
        if no selection is provided, some fields will be marked as "NA"
        and some libraries info will be missing (like size of the library)
    """
    from pathlib import Path

    if out_path is None:
        _out_path = Path("./decode_report.html")
    else:
        _out_path = Path(out_path)

    if selection is not None:
        libraries_ = [(lib.library_id, str(lib.library_size)) for lib in selection.library_collection.libraries]
    else:
        libraries_ = [(library_id, "NA") for library_id in stats.num_seqs_decoded_per_lib.keys()]

    _seq_count_data: dict = {
        **{
            library_id: (
                "{:,}".format(int(library_size)),
                "{:,}".format(stats.num_seqs_decoded_per_lib.get(library_id, 0)),
            )
            for (library_id, library_size) in libraries_
        },
    }

    _all_colors = [
        "#F5A623",
        "#F8E71C",
        "#F87171",
        "#EF4444",
        "#DC2626",
        "#F97316",
        "#3B82F6",
    ]
    _all_labels = [
        "Read to Small",
        "Read to Large",
        "Failed to Align Read",
        "Failed Library Call",
        "Failed Building Block Call",
        "Failed to Assign UMI",
        "Decoded",
    ]
    _all_values: list[int] = [
        stats.num_failed_too_short,
        stats.num_failed_too_long,
        stats.num_failed_alignment,
        stats.num_failed_library_call,
        stats.num_failed_building_block_call,
        stats.num_failed_umi,
        stats.num_seqs_decoded,
    ]

    picked_legend = []
    for _label, _color, _val in zip(_all_labels, _all_colors, _all_values, strict=False):
        if _val > 0:
            picked_legend.append({"label": _label, "color": _color})

    jinja_data = {
        "selection": "NA",
        "run_date": "NA",
        "target": "NA",
        "sequence_files": "NA",
        "timestamp": datetime.now().strftime("%b %d, %Y %H:%M"),
        "libs": ", ".join(sorted(lib[0] for lib in libraries_)),
        "tool_compounds": [],
        "num_reads": stats.num_seqs_read,
        "num_decoded": stats.num_seqs_decoded,
        "seq_counts": _seq_count_data,
        "decoding_pie": _generate_calling_pie_chart(stats),
        "decoding_bar": _generate_calling_bar_chart(stats),
        "decoding_legend": picked_legend,
    }

    if selection is not None:
        jinja_data["selection"] = selection.selection_id
        if isinstance(selection, SequencedSelection):
            jinja_data["sequence_files"] = "- " + "<br>- ".join([str(f) for f in selection.sequence_files])

        if len(selection.tool_compounds) > 0:
            tool_compound_ids = [comp.compound_id for comp in selection.tool_compounds]
            jinja_data["tool_compounds"] = "- " + "<br>- ".join(tool_compound_ids)

    # write the report as a rendered jinja2 template
    # the template is stored as a resource in the package under templates
    # with resources.path("deli.templates", "decode_report.html") as template_path:
    with resources.as_file(resources.files("deli.templates") / "decode_report.html") as template_path:
        template = jinja2.Template(open(template_path).read())
        rendered_content = template.render(jinja_data)
        with open(_out_path, "w", encoding="utf-8") as f:
            f.write(rendered_content)
