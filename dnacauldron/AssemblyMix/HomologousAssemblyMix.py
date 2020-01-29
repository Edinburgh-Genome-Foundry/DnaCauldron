from copy import deepcopy
import matplotlib.pyplot as plt

from .AssemblyMix import AssemblyMix
from ..Fragment.HomologousFragment import HomologousFragment


class HomologousAssemblyMix(AssemblyMix):
    """In construction. Do not use."""

    def __init__(
        self,
        parts,
        homology_checker,
        name="homology_mix",
        annotate_fragments_with_parts=True,
    ):

        self.parts = deepcopy(parts)
        self.name = name
        self.homology_checker = homology_checker
        self.annotate_fragments_with_parts = annotate_fragments_with_parts
        self.fragment_filters = ()
        self.initialize()

    def compute_fragments(self):
        self.fragments = []
        for part in self.parts:
            fragment = HomologousFragment.from_biopython_record(part)
            self.annotate_fragment_with_part(fragment)
            self.fragments.append(fragment)


    def assemble(
        self, fragments, circularize=False, annotate_homologies=False
    ):
        return HomologousFragment.assemble(
            fragments,
            homology_checker=self.homology_checker,
            circularize=circularize,
            annotate_homologies=annotate_homologies,
        )

    def will_clip_in_this_order(self, fragment1, fragment2):
        """Return True iff f1's right sticky end fits f2's left."""
        return fragment1.will_clip_in_this_order_with(
            fragment2, homology_checker=self.homology_checker
        )

    def plot_graphs(self, report_root, assembly, with_overhangs=True):
        file_prefix = assembly.name + "_"
        ax = self.plot_connections_graph()
        f = report_root._file(file_prefix + "connections_graph.pdf")
        ax.figure.savefig(f.open("wb"), format="pdf", bbox_inches="tight")
        plt.close(ax.figure)
