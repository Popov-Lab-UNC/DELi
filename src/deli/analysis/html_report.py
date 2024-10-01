import os
from datetime import datetime

import plotly.graph_objects as go
from jinja2 import Template

import numpy as np
import polars as pl

from deli.constants import PACKAGE_DIR


class HTMLReport:
    def __init__(self):
        self.from_file: str = "NA"
        self.error_tolerance: int = 0
        self.pattern: str = "NA"

        self.seq_lengths: list[int] = []

        self.total_reads: int = 0
        self.right_size: int = 0
        self.num_matches_found: int = 0
        self.num_matches_called: int = 0

        self.index_used: list[str] = []
        self.lib_used: list[str] = []
        self.called_csv_path: str = ""
        self.index_counts: list[dict] = []
        self.lib_counts: list[dict] = []

    def get_html_report(self, out_path=None):
        self._process_called_csv()
        if out_path is None:
            out_path = os.path.join(os.getcwd(), f"DELi_report_{datetime.now().strftime('%m-%d-%Y')}.html")

        jinja_data = {
            "file": self.from_file,
            "pattern": self.pattern,
            "error": self.error_tolerance,
            "libs": self.lib_used,
            "idx": self.index_used,
            "numseq": len(self.seq_lengths)
        }
        jinja_data.update(self._generate_calling_pie_chart())
        jinja_data.update(self._generate_seq_length_hist())
        jinja_data.update(self._generate_index_pie_chart())
        jinja_data.update(self._generate_lib_pie_chart())

        jinja_template = Template(open(os.path.join(PACKAGE_DIR, "templates", "report.html")).read())

        with open(out_path, "w", encoding="utf-8") as f:
            f.write(jinja_template.render(jinja_data))

    def _generate_calling_pie_chart(self):
        _total = self.total_reads
        _too_big = _total - self.right_size
        _failed_match = (_total - _too_big) - self.num_matches_found
        _failed_calling = (_total - _too_big - _failed_match) - self.num_matches_called
        _passed = (_total - _too_big - _failed_match - _failed_calling)

        fig = go.Figure(data=[go.Pie(labels=[">1000bp", "Failed matching", "Failed calling", "Called"],
                                     values=[_too_big, _failed_match, _failed_calling, _passed],
                                     hole=0.3, pull=[0, 0, 0, 0.2],
                                     hoverinfo='skip', textinfo='percent+value+label', textfont_size=18,
                                     marker=dict(colors=["orange", "red", "pink", "royalblue"],
                                                 line=dict(color='#000000', width=2)))])

        fig.update_layout(margin=dict(t=0, b=0, l=0, r=0), autosize=False, showlegend=False, width=850, height=650)

        return {"pie": fig.to_html(full_html=False)}

    def _generate_seq_length_hist(self):
        fig = go.Figure(data=[go.Histogram(x=np.clip(self.seq_lengths, 0, 1000))])
        fig.update_layout(margin=dict(t=0, b=0, l=0, r=0), autosize=False, showlegend=False, width=1000, height=400)

        return {"hist": fig.to_html(full_html=False)}

    def _process_called_csv(self):
        self.index_counts = pl.read_csv(self.called_csv_path).group_by("CALLED_INDEX").len().to_dicts()
        self.lib_counts = pl.read_csv(self.called_csv_path).group_by("CALLED_LIB").len().to_dicts()

    def _generate_index_pie_chart(self):

        _labels = []
        _values = []
        for _dict in self.index_counts:
            _labels.append(_dict["CALLED_INDEX"])
            _values.append(_dict["len"])

        fig = go.Figure(data=[go.Pie(labels=_labels, values=_values, hole=0.3, hoverinfo='skip',
                                     textinfo='percent+value+label', textfont_size=18,
                                     marker=dict(line=dict(color='#000000', width=2)))])

        fig.update_layout(margin=dict(t=0, b=0, l=0, r=0), autosize=False, showlegend=False, width=500, height=400)

        return {"indexpie": fig.to_html(full_html=False)}

    def _generate_lib_pie_chart(self):

        _labels = []
        _values = []
        for _dict in self.lib_counts:
            _labels.append(_dict["CALLED_LIB"])
            _values.append(_dict["len"])

        fig = go.Figure(data=[go.Pie(labels=_labels, values=_values, hole=0.3, hoverinfo='skip',
                                     textinfo='percent+value+label', textfont_size=18,
                                     marker=dict(line=dict(color='#000000', width=2)))])

        fig.update_layout(margin=dict(t=0, b=0, l=0, r=0), autosize=False, showlegend=False, width=500, height=400)

        return {"libpie": fig.to_html(full_html=False)}
