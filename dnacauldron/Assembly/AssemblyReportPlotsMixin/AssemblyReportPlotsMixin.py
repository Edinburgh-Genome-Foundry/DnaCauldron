import matplotlib.pyplot as plt
from dna_features_viewer import BiopythonTranslator
from ...biotools import write_record, record_is_linear
from .AssemblyPlotTranslator import AssemblyPlotTranslator
from .plot_cuts import plot_cuts

class AssemblyReportPlotsMixin:

    def plot_construct(self, construct, directory):
        target = directory._file(construct.id + ".pdf").open("wb")
        gr_record = AssemblyPlotTranslator().translate_record(construct)
        ax, gr = gr_record.plot(figure_width=16)
        ax.set_title(construct.id)
        ax.set_ylim(top=ax.get_ylim()[1] + 1)
        ax.figure.savefig(target, format="pdf", bbox_inches="tight")
        plt.close(ax.figure)

    def plot_provided_parts(self, report, parts_records, enzymes):
        provided_parts_dir = report._dir("provided_parts")
        for part in parts_records:
            linear = record_is_linear(part, default=False)
            ax, gr = plot_cuts(part, enzymes, linear=linear)
            f = provided_parts_dir._file(part.id + ".pdf").open("wb")
            ax.figure.savefig(f, format="pdf", bbox_inches="tight")
            plt.close(ax.figure)
            gb_file = provided_parts_dir._file(part.id + ".gb")
            write_record(part, gb_file, "genbank")

    def plot_fragments(self, mix, report):
        fragments_dir = report._dir("fragments")
        seen_fragments = {}
        for fragment in mix.fragments:
            gr = BiopythonTranslator().translate_record(fragment)
            ax, _ = gr.plot(strand_in_label_threshold=7)
            name = self.name_fragment(fragment)
            if name not in seen_fragments:
                seen_fragments[name] = 0
            seen_fragments[name] += 1
            file_name = "%s_%02d.pdf" % (name, seen_fragments[name])
            ax.figure.savefig(
                fragments_dir._file(file_name).open("wb"),
                format="pdf",
                bbox_inches="tight",
            )
            plt.close(ax.figure)

    def plot_slots_graph(self, mix, report, highlighted_parts):
        ax = mix.plot_slots_graph(
            with_overhangs=self.show_overhangs_in_graph,
            show_missing=True,
            highlighted_parts=highlighted_parts,
        )
        f = report._file("parts_graph.pdf")
        ax.figure.savefig(f.open("wb"), format="pdf", bbox_inches="tight")
        plt.close(ax.figure)

    def plot_connections_graph(self, report, mix):
        ax = mix.plot_connections_graph()
        f = report._file("connections_graph.pdf")
        ax.figure.savefig(f.open("wb"), format="pdf", bbox_inches="tight")
        plt.close(ax.figure)