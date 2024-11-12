"""handles making html reports for deli decoding"""

import importlib.resources as resources
import json
import os
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime
from typing import List, Self, Union

import jinja2
import plotly.graph_objects as go


@dataclass
class DecodeReportStats:
    """Dataclass to hold info for decoding statistics"""

    num_reads: int
    read_length: Counter
    num_match_attempts: int
    num_valid_matches: int
    num_match_too_big: int
    num_match_too_small: int
    num_reads_with_match: int
    num_call_attempts: int
    num_valid_calls: int
    num_reads_with_calls: int
    experiment_name: str
    libraries: Counter
    indexes: Counter = field(default_factory=Counter)

    @classmethod
    def load_report_file(cls, path: Union[os.PathLike, str]) -> Self:
        """
        Load a `DecodeReportStats` instance from a json file

        Notes
        -----
        File must have the report stats json format and nothing else

        Parameters
        ----------
        path: Union[os.PathLike, str]
            path to file

        Returns
        -------
        DecodeReportStats
        """
        info = json.load(open(path, "r"))

        return cls(
            num_reads=info["num_reads"],
            read_length=Counter({int(key): val for key, val in info["read_length"].items()}),
            num_match_attempts=info["num_match_attempts"],
            num_valid_matches=info["num_valid_matches"],
            num_match_too_big=info["num_match_too_big"],
            num_match_too_small=info["num_match_too_small"],
            num_reads_with_match=info["num_reads_with_match"],
            num_call_attempts=info["num_call_attempts"],
            num_valid_calls=info["num_valid_calls"],
            num_reads_with_calls=info["num_reads_with_calls"],
            experiment_name=info["experiment_name"],
            libraries=Counter(info["libraries"]),
            indexes=Counter(info["indexes"]),
        )

    def __add__(self, other):
        """Adds to report stats together"""
        if not isinstance(other, DecodeReportStats):
            raise TypeError(
                f"unsupported operand type(s) for +: "
                f"'{self.__class__.__name__}' and '{other.__class__.__name__}'"
            )
        else:
            if self.experiment_name != other.experiment_name:
                raise RuntimeError(
                    f"cannot add two {self.__class__.__name__} object with different "
                    f"'experiment_name': '{self.experiment_name}' and '{other.experiment_name}'"
                )
            return self.__class__(
                num_reads=self.num_reads + other.num_reads,
                read_length=self.read_length + other.read_length,
                num_match_attempts=self.num_match_attempts + other.num_match_attempts,
                num_valid_matches=self.num_valid_matches + other.num_valid_matches,
                num_match_too_big=self.num_match_too_big + other.num_match_too_big,
                num_match_too_small=self.num_match_too_small + other.num_match_too_small,
                num_reads_with_match=self.num_reads_with_match + other.num_reads_with_match,
                num_call_attempts=self.num_call_attempts + other.num_reads_with_calls,
                num_valid_calls=self.num_valid_calls + other.num_valid_calls,
                num_reads_with_calls=self.num_reads_with_calls + other.num_reads_with_calls,
                experiment_name=self.experiment_name,
                libraries=self.libraries + other.libraries,
                indexes=self.indexes + other.indexes,
            )


#
# def build_report_str(
#     num_reads: int,
#     num_match_attempts: int,
#     num_valid_matches: int,
#     num_match_too_big: int,
#     num_match_too_small: int,
#     num_reads_with_match: int,
#     num_call_attempts: int,
#     num_valid_calls: int,
#     num_reads_with_calls: int,
#     experiment_name: str,
#     libraries: Union[Counter, Dict[str, int]],
#     indexes: Union[Counter, Dict[str, int]],
# ):
#     """
#     Given various deli decode stats, generate a report string
#
#     Parameters
#     ----------
#     num_reads: int
#         number of reads from this decoding run
#     num_match_attempts: int
#         number of matches that were attempted (including revcomp)
#     num_valid_matches: int
#         number of matches that were successful
#     num_match_too_big: int
#         number of matches that failed for the sequence being too big
#     num_match_too_small: int
#         number of matches that failed for the sequence being too small
#     num_reads_with_match: int
#         the number of unique reads that had a match in them
#     num_call_attempts: int
#         the number of calls that were attempted
#     num_valid_calls: int
#         the number of valid calls
#     num_reads_with_calls: int
#         the number of unique reads that had a call
#     experiment_name: str
#         the name of the experiment
#     libraries: Union[Counter, Dict[str, int]]
#         a count dict of lib_id -> number of calls with that lib_id
#         can include "FAILED" lib_id for failed lib calls
#     indexes: Union[Counter, Dict[str, int]]
#         a count dict of lib_id -> number of calls with that lib_id
#         can include "FAILED" lib_id for failed lib calls
#
#     Returns
#     -------
#     str
#     """
#     _str = (
#         f"{num_reads}\n"
#         f"{num_match_attempts}\n"
#         f"{num_valid_matches}\n"
#         f"{num_match_too_big}\n"
#         f"{num_match_too_small}\n"
#         f"{num_reads_with_match}\n"
#         f"{num_call_attempts}\n"
#         f"{num_valid_calls}\n"
#         f"{num_reads_with_calls}\n"
#         f"{experiment_name}\n"
#     )
#     _str += "|".join([f"{lib}:{count}" for lib, count in libraries.items()])
#     _str += "|".join([f"{idx}:{count}" for idx, count in indexes.items()])
#     return _str


def build_decoding_report(
    report_stats: Union[List[DecodeReportStats], DecodeReportStats],
    out_path: Union[str, os.PathLike],
):
    """
    Generates a deli decoding html report from report stats

    Notes
    -----
    Will render the report using the decode_report.html resource
    Can take in several `DecodeReportStats` objects as long as
    they all share the same experiment_name

    Parameters
    ----------
    report_stats: Union[List[DecodeReportStats], DecodeReportStats]
        the reporting stats object(s) to generate the report for
    out_path: Union[str, os.PathLike]
        full path (with file name) to save the report to
        will raise exception if directory used does not exist
    """
    if isinstance(report_stats, list):
        _all_report_stats = report_stats[0]
        for _stats in report_stats[1:]:
            _all_report_stats += _stats
    else:
        _all_report_stats = report_stats

    def _generate_calling_pie_chart(
        num_reads,
        num_match_too_big,
        num_match_too_small,
        num_reads_with_match,
        num_reads_with_calls,
    ):
        """Generate the main calling results pie chart html"""
        _total = num_reads
        _too_big = num_match_too_big
        _too_small = num_match_too_small
        _failed_match = num_reads - num_reads_with_match
        _failed_calling = num_reads_with_match - num_reads_with_calls
        _passed = num_reads_with_calls

        fig = go.Figure(
            data=[
                go.Pie(
                    labels=[
                        "Read to large",
                        "Read to small",
                        "Failed matching",
                        "Failed calling",
                        "Called",
                    ],
                    values=[_too_big, _too_small, _failed_match, _failed_calling, _passed],
                    hole=0.3,
                    pull=[0, 0, 0, 0, 0.2],
                    hoverinfo="skip",
                    textinfo="percent+value+label",
                    textfont_size=18,
                    marker=dict(
                        colors=["orange", "yellow", "red", "pink", "royalblue"],
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

    def _generate_seq_length_hist(_seq_lengths):
        """Generate the html for the seq length histogram"""
        x_vals = list(range(min(_seq_lengths.keys()), max(_seq_lengths.keys()) + 1))
        y_vals = [_seq_lengths.get(key, 0) for key in x_vals]
        fig = go.Figure(data=[go.Bar(x=x_vals, y=y_vals)])
        fig.update_layout(
            xaxis_title="Read Length",
            yaxis_title="Count",
            margin=dict(t=0, b=0, l=0, r=0),
            autosize=False,
            showlegend=False,
            width=1000,
            height=400,
        )
        fig.update_xaxes(range=[0, 1000])
        return fig.to_html(full_html=False, include_plotlyjs=False)

    def _generate_sub_call_pie_chart(count_dict):
        """Generate the html for a sub-call (library or index) pie chart"""
        _labels = []
        _values = []
        for idx, count in count_dict.items():
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

    jinja_data = {
        "timestamp": datetime.now().strftime("%b %d, %Y %H:%M"),
        "experiment": _all_report_stats.experiment_name,
        "libs": ", ".join(
            sorted([lib for lib in _all_report_stats.libraries.keys() if lib != "FAILED"])
        ),
        "idx": ", ".join(
            sorted([idx for idx in _all_report_stats.indexes.keys() if idx != "FAILED"])
        )
        if len(_all_report_stats.indexes) > 0
        else "<b>no de-multiplexing occurred during this DELi decoding run</b>",
        "num_reads": _all_report_stats.num_reads,
        "read_hist": _generate_seq_length_hist(_all_report_stats.read_length),
        "pie": _generate_calling_pie_chart(
            num_reads=_all_report_stats.num_reads,
            num_match_too_big=_all_report_stats.num_match_too_big,
            num_match_too_small=_all_report_stats.num_match_too_small,
            num_reads_with_match=_all_report_stats.num_reads_with_match,
            num_reads_with_calls=_all_report_stats.num_reads_with_calls,
        ),
        "indexpie": _generate_sub_call_pie_chart(_all_report_stats.indexes),
        "libpie": _generate_sub_call_pie_chart(_all_report_stats.libraries),
    }

    with resources.path("templates", "decode_report.html") as template_path:
        template = jinja2.Template(open(template_path).read())
        rendered_content = template.render(jinja_data)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(rendered_content)
