
from copy import deepcopy
from Bio import Restriction
from ...biotools import record_is_linear, annotate_record
from dna_features_viewer import BiopythonTranslator, GraphicRecord, CircularGraphicRecord


class EnzymesSitesTranslator(BiopythonTranslator):
    def compute_feature_label(self, feature):
        if abs(feature.location.end - feature.location.start) > 100:
            label = BiopythonTranslator.compute_feature_label(
                self, feature
            )
            return self.shorten_string(label, 10)
        else:
            return feature.qualifiers.get("enzyme", None)

    def compute_feature_color(self, feature):
        if feature.qualifiers.get("enzyme", False) and (
            feature.type == "misc_feature"
        ):
            return "#f5eaff"
        else:
            return "#fefefe"
    
    @staticmethod
    def shorten_string(string, max_length=30):
        """Truncate and add '...' if the string is too long"""
        suffix = "" if len(string) < max_length else "..."
        return string[:max_length] + suffix

def plot_cuts(record, enzymes, linear="auto", figure_width=5, ax=None):
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
    record = deepcopy(record)
    if linear == "auto":
        linear = record_is_linear(record, True)
    for enzyme_name in enzymes:
        enzyme = Restriction.__dict__[enzyme_name]
        cuts = enzyme.search(record.seq, linear=linear)
        for cut in cuts:
            annotate_record(
                record,
                (cut, cut + 1),
                feature_type="misc_feature",
                enzyme=enzyme_name,
            )

    translator = EnzymesSitesTranslator()
    record_class = GraphicRecord if linear else CircularGraphicRecord
    graphic_record = translator.translate_record(
        record, record_class=record_class
    )
    return graphic_record.plot(ax=ax, figure_width=figure_width)