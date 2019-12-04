from copy import deepcopy

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from Bio import Restriction
from dna_features_viewer import (
    BiopythonTranslator,
    CircularGraphicRecord,
    GraphicRecord,
)
from dnacauldron.tools import annotate_record

try:
    import pygraphviz
    from networkx.drawing.nx_agraph import graphviz_layout

    GRAPHVIZ_AVAILABLE = True
except ImportError:
    GRAPHVIZ_AVAILABLE = False

from ..tools import record_is_linear


class AssemblyTranslator(BiopythonTranslator):
    """Custom theme for plotting GENBANK assemblies w. Dna Features Viewer."""

    def is_source(self, feature):
        return (feature.type == "misc_feature") and feature.qualifiers.get(
            "source", False
        )

    def compute_feature_color(self, feature):
        if feature:
            if self.is_source(feature):
                return "#ff4c4c"
            else:
                return "#f9edbb"

    def compute_feature_label(self, feature):
        if self.is_source(feature):
            return "".join(feature.qualifiers["source"])
        elif abs(feature.location.end - feature.location.start) > 100:
            label = BiopythonTranslator.compute_feature_label(self, feature)
            return abreviate_string("".join(label), 30)
        else:
            return None


def abreviate_string(string, max_length=30):
    """Truncate and add '...' if the string is too long"""
    return string[:max_length] + ("" if len(string) < max_length else "...")


def plot_cuts(record, enzyme_name, linear="auto", figure_width=5, ax=None):
    """Plot a construct and highlight where an enzyme cuts.

    Parameters
    ----------

    record
      The biopython record to be plotted

    enzyme_name
      Name of the enzyme, e.g. "EcoRI"

    linear
      True for a linear construct, False for a circular construct

    figure_width
      Width of the figure in inches.

    ax
      Matplotlib ax on which to plot the construct. If None is provided, one
      will be created an returned.


    """
    enzyme = Restriction.__dict__[enzyme_name]
    record = deepcopy(record)
    if linear == "auto":
        linear = record_is_linear(record, True)
    cuts = enzyme.search(record.seq, linear=linear)
    for cut in cuts:
        annotate_record(
            record,
            (cut, cut + 1),
            feature_type="misc_feature",
            enzyme=enzyme_name,
        )

    class MyTranslator(BiopythonTranslator):
        def compute_feature_label(self, feature):
            if abs(feature.location.end - feature.location.start) > 100:
                label = BiopythonTranslator.compute_feature_label(
                    self, feature
                )
                return abreviate_string(label, 10)
            else:
                return feature.qualifiers.get("enzyme", None)

        def compute_feature_color(self, feature):
            if feature.qualifiers.get("enzyme", False) and (
                feature.type == "misc_feature"
            ):
                return "#f5eaff"
            else:
                return "#fefefe"

    translator = MyTranslator()
    record_class = GraphicRecord if linear else CircularGraphicRecord
    graphic_record = translator.translate_record(
        record, record_class=record_class
    )
    return graphic_record.plot(ax=ax, figure_width=figure_width)