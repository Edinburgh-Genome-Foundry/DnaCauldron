from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

try:
    # Biopython <1.78
    from Bio.Alphabet import DNAAlphabet

    has_dna_alphabet = True
except ImportError:
    # Biopython >=1.78
    has_dna_alphabet = False

from ...biotools import (
    set_record_topology,
    crop_record_with_saddling_features,
    sequence_to_biopython_record,
    annotate_record,
)
from ..Fragment import Fragment
from .StickyEnd import StickyEnd
from .StickyEndSeq import StickyEndSeq


class StickyEndFragment(Fragment):
    """Biopython SeqRecord whose sequence has sticky ends."""

    def will_clip_in_this_order_with(self, other):
        """Return True if this record's right sticky end is complementary with
        the other record's left sticky end."""
        right_end = self.seq.right_end
        return (right_end is not None) and right_end.will_clip_directly_with(
            other.seq.left_end
        )

    def circularized(
        self,
        annotate_homology=False,
        annotation_type="homology",
        qualifiers=None,
    ):
        """Return the Biopython record obtained by cirularizing the result.

        Only works if the left and right sticky ends are compatible. The
        return is a simple Biopython record where the sticky end has been
        integrated in the sequence.
        """
        if not self.will_clip_in_this_order_with(self):
            raise ValueError(
                "Only constructs with two compatible sticky ends" " can be circularized"
            )
        connector = SeqRecord(Seq(str(self.seq.left_end)))
        if annotate_homology:
            self.annotate_connector(connector, annotation_type=annotation_type)
        result = connector + self
        result.__class__ = SeqRecord
        set_record_topology(result, "circular")

        return result

    def annotate_connector(self, connector, annotation_type="homology"):
        """Annotate a connector to indicate it used to be a sticky end."""
        if len(connector) > 8:
            label = "homology"
        else:
            label = str(connector.seq)
        feature = self.create_homology_annotation(
            start=0,
            end=len(connector),
            annotation_type=annotation_type,
            label=label,
        )
        connector.features = [feature]

    @staticmethod
    def assemble(fragments, circularize=False, annotate_homologies=False):
        """Return the (sticky end) record obtained by assembling the fragments.

        Parameters
        ----------

        fragments
          List of StickyEndFragments to assemble.

        circularize
          True to also assemble the end flanks of the final construct (results
          in a Biopython Record), false to not do it (the result is then a
          StickyEndFragment).

        annotate_homologies
           If true, homologies will have an annotation in the final, predicted
           construct records.
        """
        result = fragments[0]
        for fragment in fragments[1:]:
            result = result.assemble_with(
                fragment, annotate_homology=annotate_homologies
            )
        if circularize:
            result = result.circularized(annotate_homology=annotate_homologies)

        if has_dna_alphabet:  # Biopython <1.78
            result.seq.alphabet = DNAAlphabet()
        result.annotations["molecule_type"] = "DNA"

        return result

    def assemble_with(self, other, annotate_homology=False, annotation_type="homology"):
        connector_str = str(self.seq.right_end)
        connector = SeqRecord(Seq(connector_str))
        if annotate_homology:
            self.annotate_connector(connector, annotation_type=annotation_type)
        selfc = SeqRecord(
            seq=Seq(str(self.seq)),
            features=self.features,
            annotations=self.annotations,
        )
        new_record = SeqRecord.__add__(selfc, connector).__add__(other)
        new_record.seq = self.seq + other.seq
        new_record.__class__ = StickyEndFragment

        if has_dna_alphabet:  # Biopython <1.78
            new_record.seq.alphabet = DNAAlphabet()
        new_record.annotations["molecule_type"] = "DNA"

        return new_record

    @staticmethod
    def list_from_record_digestion(record, enzyme, linear="auto"):
        if linear == "auto":
            linear = record.annotations.get("topology", "linear") == "linear"
        if isinstance(enzyme, (list, tuple)):
            n_cuts = sum([len(e.search(record.seq, linear=linear)) for e in enzyme])
        else:
            n_cuts = len(enzyme.search(record.seq, linear=linear))
        if n_cuts == 0:
            return [record]
        if not linear:
            record.features = [f for f in record.features if f.location is not None]
            record_fragments = StickyEndFragment.list_from_record_digestion(
                record + record, enzyme=enzyme, linear=True
            )
            return record_fragments[1 : n_cuts + 1]
        fragments = StickyEndSeq.list_from_sequence_digestion(
            record.seq, enzyme, linear=linear
        )
        record_fragments = []
        for fragment in fragments:
            index = record.seq.upper().find(fragment)
            if index == -1:
                continue

            def only_parts_indicators(feature):
                return feature.qualifiers.get("indicates_part", False)

            subrecord = crop_record_with_saddling_features(
                record=record,
                start=index,
                end=index + len(fragment),
                filters=(only_parts_indicators,),
            )
            new_stickyend_record = StickyEndFragment(
                seq=fragment,
                features=subrecord.features,
                annotations=subrecord.annotations,
            )
            record_fragments.append(new_stickyend_record)

        return record_fragments

    def to_standard_string(self):
        """Return a string representation of the fragment, used for quick
        comparison of fragments and fragments chains."""
        return "%s%s%s" % (self.seq.left_end, self.seq, self.seq.right_end)

    def text_representation_in_plots(self):
        """Plot a fragment as left//PART_NAME//right (where // is a new line)"""
        lines = [
            str(self.seq.left_end),
            r"$\bf{%s}$" % self.original_part.id,
            str(self.seq.right_end),
        ]
        return "\n".join(lines)

    def as_biopython_record(self):
        record_left = self.seq.left_end.as_biopython_record()
        record_right = self.seq.right_end.as_biopython_record()
        record = record_left + self + record_right
        record.id = self.id
        return record
