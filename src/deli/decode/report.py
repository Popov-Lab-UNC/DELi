"""handles making html reports for deli decoding"""

import importlib.resources as resources
import os
from datetime import datetime
from typing import Optional

import jinja2
import plotly.graph_objects as go

from deli.selection import DELSelection, Selection, SequencedSelection

from .decoder import DecodeStatistics


def _generate_calling_pie_chart(
    decode_stats: DecodeStatistics,
):
    """Generate the main calling results pie chart HTML"""
    _all_colors = [
        "orange",
        "yellow",
        "firebrick",
        "crimson",
        "royalblue",
    ]
    _all_labels = [
        "Read to Small",
        "Read to Large",
        "Failed Library Call",
        "Failed Building Block Call",
        "Decoded",
    ]
    _all_values: list[int] = [
        decode_stats.num_failed_too_short,
        decode_stats.num_failed_too_long,
        decode_stats.num_failed_library_call,
        decode_stats.num_failed_building_block_call + decode_stats.num_failed_alignment,
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
    selection: DELSelection or SequencedSelection, optional
        the selection to build the report for
        if DEL selection provided will include selection and library info.
        if sequenced selection is provided, will also include sequence file info.
        if no selection is provided, some fields will be marked as "NA"
        and some libraries info will be missing (like size of the library)
    """
    from pathlib import Path

    if out_path is None:
        _out_path = Path("./decode_report.html")
    else:
        _out_path = Path(out_path)

    if isinstance(selection, DELSelection):
        libraries_ = [(lib.library_id, str(lib.library_size)) for lib in selection.library_collection.libraries]
    else:
        libraries_ = [(library_id, "NA") for library_id in stats.num_seqs_decoded_per_lib.keys()]

    _seq_count_data: dict = {
        **{
            library_id: (
                "{:,}".format(library_size) if isinstance(library_size, int) else library_size,
                "{:,}".format(stats.num_seqs_decoded_per_lib.get(library_id, 0)),
            )
            for (library_id, library_size) in libraries_
        },
    }

    jinja_data = {
        "selection": "NA",
        "run_date": "NA",
        "target": "NA",
        "sequence_files": "NA",
        "timestamp": datetime.now().strftime("%b %d, %Y %H:%M"),
        "libs": ", ".join(sorted(l[0] for l in libraries_)),
        "tool_compounds": [],
        "num_reads": stats.num_seqs_read,
        "num_decoded": stats.num_seqs_decoded,
        "seq_counts": _seq_count_data,
        "decoding_pie": _generate_calling_pie_chart(stats),
    }

    if selection is not None:
        jinja_data["selection"] = selection.selection_id
        jinja_data["run_date"] = selection.get_run_date_as_str()
        jinja_data["target"] = _default_to_NA(selection.selection_condition.target_id)
        if isinstance(selection, SequencedSelection):
            jinja_data["sequence_files"] = "- " + "<br>- ".join([str(f) for f in selection.sequence_files])

        if len(selection.tool_compounds) > 0:
            tool_compound_ids = [comp.compound_id for comps in selection.tool_compounds for comp in comps.compounds]
            jinja_data["tool_compounds"] = "- " + "<br>- ".join(tool_compound_ids)

    # write the report as a rendered jinja2 template
    # the template is stored as a resource in the package under templates
    # with resources.path("deli.templates", "decode_report.html") as template_path:
    with resources.as_file(resources.files("deli.templates") / "decode_report.html") as template_path:
        template = jinja2.Template(open(template_path).read())
        rendered_content = template.render(jinja_data)
        with open(_out_path, "w", encoding="utf-8") as f:
            f.write(rendered_content)
