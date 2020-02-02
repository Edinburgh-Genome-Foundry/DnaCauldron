from Bio.SeqFeature import SeqFeature, FeatureLocation
from ..Fragment import Fragment
from ..biotools import set_record_topology
from .AssemblyMix import AssemblyMix
import matplotlib.pyplot as plt


class LigaseCyclingReactionMix(AssemblyMix):
    """Mix to represent and simulate Ligase Cycling Reactions.

    parts
      List of parts records

    bridging_oligos
      List of bridging_oligos records (in direct sense)

    homology_checker
      An HomologyChecker instance defining which homology sizes and melting
      temperatures are valid between one bridging oligo and one part.
    """

    def __init__(
        self,
        parts,
        bridging_oligos,
        homology_checker="default",
        name="lcr_mix",
        annotate_fragments_with_parts=True,
    ):
        self.parts = parts
        self.homology_checker = homology_checker
        self.bridging_oligos = bridging_oligos
        self.reversed_bridging_oligos = []
        for oligo in bridging_oligos:
            reverse = oligo.reverse_complement()
            reverse.id = oligo.id
            self.reversed_bridging_oligos.append(oligo)
        self.name = name
        self.annotate_fragments_with_parts = annotate_fragments_with_parts
        self.fragment_filters = ()
        self.initialize()

    def compute_fragments(self):
        self.fragments = list(self.parts)
        self.fragments = []
        for part in self.parts:
            fragment = Fragment.from_biopython_record(part)
            fragment.original_part = part
            self.annotate_fragment_with_part(fragment)
            self.fragments.append(fragment)

    def assemble(self, fragments, circularize=False, annotate_homologies=False):
        """Assemble sticky-end fragments into a single one (sticky or not).

        Parameters
        ----------

        fragments
          List of fragments records

        circularize
          If True and if the two ends of the final assembly are compatible,
          circularize the construct, i.e. return a non-sticky record
          representing the circular assembly of the fragments.

        annotate_homologies
          If True, all homology regions that where formerly sticky ends will
          be annotated in the final record.
        """
        result = fragments[0]
        oligos = []
        for fragment in fragments[1:]:
            homology = self.find_oligo_homology(result, fragment)
            if homology is None:
                raise ValueError("Homologies should be there at this stage")
            oligo, (start, end) = homology
            oligos.append(oligo)
            result += fragment
            if annotate_homologies:
                annotation = self.create_homology_annotation(
                    start, end, label=oligo.id, annotation_type="homology"
                )
                result.features.append(annotation)

        set_record_topology(result, "circular" if circularize else "linear")
        result.fragments = fragments + oligos
        return result

    def find_oligo_homology(self, fragment1, fragment2):
        """Return (oligo, (start, end)), or None if no bridging oligo is found.
        """
        sequence = str((fragment1 + fragment2).seq)
        L1 = len(fragment1)
        for oligo in self.reversed_bridging_oligos:
            oligo_str = str(oligo.seq)
            if oligo_str in sequence:
                start = sequence.index(oligo_str)
                end = start + len(oligo_str)
                if not start <= L1 < end:
                    continue
                homology1 = sequence[start:L1]
                hom1_valid = self.homology_checker.check_homology(homology1)
                homology2 = sequence[L1:end]
                hom2_valid = self.homology_checker.check_homology(homology2)
                if hom1_valid and hom2_valid:
                    return oligo, (start, end)
        return None

    def will_clip_in_this_order(self, fragment1, fragment2):
        """Return True iff f1's right sticky end fits f2's left."""
        homology = self.find_oligo_homology(fragment1, fragment2)
        return homology is not None

    def create_homology_annotation(
        self, start, end, label, annotation_type, color="#f7e8f7"
    ):
        qualifiers = {
            "label": label,
            "color": color,
            "ApEinfo_fwdcolor": color,
        }
        return SeqFeature(
            FeatureLocation(start, end),
            type=annotation_type,
            qualifiers=qualifiers,
        )

    def plot_graphs(self, report_root, assembly, with_overhangs=True):
        file_prefix = assembly.name + "_"
        ax = self.plot_connections_graph()
        f = report_root._file(file_prefix + "connections_graph.pdf")
        ax.figure.savefig(f.open("wb"), format="pdf", bbox_inches="tight")
        plt.close(ax.figure)