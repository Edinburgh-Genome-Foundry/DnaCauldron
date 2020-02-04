import matplotlib.pyplot as plt
from dna_features_viewer import BiopythonTranslator
from ...biotools import write_record, record_is_linear
from .AssemblyPlotTranslator import AssemblyPlotTranslator
from .plot_cuts import plot_cuts


class AssemblyReportPlotsMixin:
    """Mixin for AssemblyReportWriter"""

    def _get_plots_options(self, errors_detected):
        """Compute whether the plots should be generated, depending on assembly
        success or failure."""
        def evaluate(condition):
            return errors_detected if condition == "on_error" else condition

        return {
            "fragment_plots": evaluate(self.include_fragment_plots),
            "mix_graphs_plots": evaluate(self.include_mix_graphs),
            "parts_plots": evaluate(self.include_part_plots),
        }

    def plot_construct(self, construct, directory):
        """Plot schemas of the predicted constructs."""
        target = directory._file(construct.id + ".pdf").open("wb")
        gr_record = AssemblyPlotTranslator().translate_record(construct)
        ax, gr = gr_record.plot(figure_width=16)
        ax.set_title(construct.id)
        ax.set_ylim(top=ax.get_ylim()[1] + 1)
        ax.figure.savefig(target, format="pdf", bbox_inches="tight")
        plt.close(ax.figure)

    def plot_provided_parts(self, report_root, parts_records, enzymes):
        """Plot schemas of the parts provided, with relevant restriction sites.
        """
        provided_parts_dir = report_root._dir("provided_parts_plots")
        for part in parts_records:
            linear = record_is_linear(part, default=False)
            ax, gr = plot_cuts(part, enzymes, linear=linear)
            f = provided_parts_dir._file(part.id + ".pdf").open("wb")
            ax.figure.savefig(f, format="pdf", bbox_inches="tight")
            plt.close(ax.figure)
