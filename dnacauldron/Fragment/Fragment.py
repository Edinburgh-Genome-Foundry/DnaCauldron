from copy import deepcopy
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dna_features_viewer import BiopythonTranslator


class Fragment(SeqRecord):
    """Base class to represent a DNA fragment that will assemble with other
    fragments, such as HomologousFragments and StickyEndFragments."""

    @staticmethod
    def from_biopython_record(biopython_record):
        """Convert a biopython record into a HomologousFragment (class change).
        """
        new_record = deepcopy(biopython_record)
        new_record.original_part = biopython_record
        new_record.__class__ = Fragment
        return new_record

    def plot(self, ax=None):
        """Plot the fragment and its features on a Matplotlib ax.

        This creates a new ax if no ax is provided. The ax is returned at the
        end.
        """
        graphic_record = BiopythonTranslator().translate_record(self)
        ax, _ = graphic_record.plot(ax=ax, strand_in_label_threshold=7)
        return ax

    def reverse_complement(self):
        """Reverse-complement the fragment while keeping its class"""
        # Note: if the fragment has a StickyEndSeq, it will be properly
        # reverse-complemented with its ends
        new_record = SeqRecord.reverse_complement(self)
        new_record.__class__ = self.__class__
        return new_record

    def to_standard_string(self):
        """Return a standard string to represent and identify the fragment.

        This method is used to standardize and recognize similar FragmentChain
        instances.
        """

        return str(self.seq)

    def create_homology_annotation(
        self, start, end, label, annotation_type, color="#f7e8f7"
    ):
        qualifiers = {
            "label": label,
            "color": color,
            "ApEinfo_fwdcolor": color,
        }
        return SeqFeature(
            FeatureLocation(start, end), type=annotation_type, qualifiers=qualifiers,
        )

    def text_representation_in_plots(self):
        return r"$\bf{%s}$" % self.original_part.id

    def as_biopython_record(self):
        return self
