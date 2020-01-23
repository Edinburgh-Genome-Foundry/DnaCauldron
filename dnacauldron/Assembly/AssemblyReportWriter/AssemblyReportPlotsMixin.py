import matplotlib.pyplot as plt
from dna_features_viewer import BiopythonTranslator
from ...biotools import write_record, record_is_linear
from .AssemblyPlotTranslator import AssemblyPlotTranslator
from .plot_cuts import plot_cuts


class AssemblyReportPlotsMixin:
    def get_plots_options(self, errors_detected):
        def evaluate(condition):
            return errors_detected if condition == "on_error" else condition

        return {
            "fragments_plots": evaluate(self.include_fragments_plots),
            "mix_graphs_plots": evaluate(self.include_mix_graphs),
            "parts_plots": evaluate(self.include_parts_plots),
        }

    def plot_construct(self, construct, directory):
        target = directory._file(construct.id + ".pdf").open("wb")
        gr_record = AssemblyPlotTranslator().translate_record(construct)
        ax, gr = gr_record.plot(figure_width=16)
        ax.set_title(construct.id)
        ax.set_ylim(top=ax.get_ylim()[1] + 1)
        ax.figure.savefig(target, format="pdf", bbox_inches="tight")
        plt.close(ax.figure)

    def plot_provided_parts(self, report_root, parts_records, enzymes):
        provided_parts_dir = report_root._dir("provided_parts_plots")
        for part in parts_records:
            linear = record_is_linear(part, default=False)
            ax, gr = plot_cuts(part, enzymes, linear=linear)
            f = provided_parts_dir._file(part.id + ".pdf").open("wb")
            ax.figure.savefig(f, format="pdf", bbox_inches="tight")
            plt.close(ax.figure)

    # def plot_slots_graph(self, mix, report_root, highlighted_parts):
    #     ax = mix.plot_slots_graph(
    #         with_overhangs=self.show_overhangs_in_graph,
    #         show_missing=True,
    #         highlighted_parts=highlighted_parts,
    #     )
    #     f = report_root._file("parts_graph.pdf")
    #     ax.figure.savefig(f.open("wb"), format="pdf", bbox_inches="tight")
    #     plt.close(ax.figure)

    # def plot_connections_graph(self, report_root, mix):
    #     ax = mix.plot_connections_graph()
    #     f = report_root._file("connections_graph.pdf")
    #     ax.figure.savefig(f.open("wb"), format="pdf", bbox_inches="tight")
    #     plt.close(ax.figure)
