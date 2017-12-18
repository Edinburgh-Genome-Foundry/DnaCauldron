from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import DNAAlphabet

from .StickyEnd import StickyEnd
from .StickyEndsSeq import StickyEndsSeq

class StickyEndsSeqRecord(SeqRecord):
    """Biopython SeqRecord whose sequence has sticky ends."""

    def will_clip_in_this_order_with(self, other):
        """Return True iff this record's right sticky end is complementary with
        the other record's left sticky end."""
        right_end = self.seq.right_end
        return ((right_end is not None) and
                right_end.will_clip_directly_with(other.seq.left_end))

    def circularized(self, annotate_homology=False, annotation_type="homology",
                     qualifiers=None):
        """Return the biopython record obtained by cirularizing the result.

        Only works if the left and right sticky ends are compatible. The
        return is a simple Biopython record where the sticky end has been
        integrated in the sequence.
        """
        if not self.will_clip_in_this_order_with(self):
            raise ValueError("Only constructs with two compatible sticky ends"
                             " can be circularized")
        connector_str = str(self.seq.left_end)
        connector = SeqRecord(Seq(str(self.seq.left_end)))
        if annotate_homology:
            label = "homology" if (len(connector) > 8) else connector_str
            connector.features = [
                SeqFeature(FeatureLocation(0, len(connector), 1),
                           type=annotation_type,
                           qualifiers={"label": label})
            ]
        result = connector + self
        result.linear = False
        return result

    @staticmethod
    def assemble(fragments, circularize=False, annotate_homologies=False):
        """Return the (sticky end) record obtained by assembling the fragments.


        """
        result = fragments[0]
        for fragment in fragments[1:]:
            result = result.assemble_with(
                fragment,
                annotate_homology=annotate_homologies
            )
        if circularize:
            result = result.circularized(annotate_homology=annotate_homologies)
        result.seq.alphabet = DNAAlphabet()
        return result

    def assemble_with(self, other, annotate_homology=False,
                      annotation_type="homology", **qualifiers):
        connector_str = str(self.seq.right_end)
        connector = SeqRecord(Seq(connector_str))
        if len(qualifiers) == 0:
            label = "homology" if (len(connector) > 8) else connector_str
            qualifiers = {"label": label}
        if annotate_homology:
            connector.features = [
                SeqFeature(FeatureLocation(0, len(connector), 1),
                           type=annotation_type,
                           qualifiers=qualifiers)
            ]
        selfc = SeqRecord(seq=Seq(str(self.seq)),
                          features=self.features,
                          annotations=self.annotations)
        new_record = SeqRecord.__add__(selfc, connector).__add__(other)
        new_record.seq = self.seq + other.seq
        new_record.__class__ = StickyEndsSeqRecord
        new_record.seq.alphabet = DNAAlphabet()
        return new_record

    def reverse_complement(self):
        new_record = SeqRecord.reverse_complement(self)
        new_record.__class__ = StickyEndsSeqRecord
        return new_record

    def __add__(self, other):
        return self.assemble_with(other)

    # def __hash__(self):
    #     return hash("StickyEndSeqRecord"+str(self.seq))


    @staticmethod
    def list_from_record_digestion(record, enzyme, linear=True):
        n_cuts = len(enzyme.search(record.seq, linear=linear))
        if n_cuts == 0:
            return [record]
        if not linear:
            record_fragments = StickyEndsSeqRecord.list_from_record_digestion(
                record + record, enzyme=enzyme, linear=True)
            return record_fragments[1:n_cuts + 1]
        fragments = StickyEndsSeq.list_from_sequence_digestion(
            record.seq, enzyme, linear=linear)
        record_fragments = []
        for fragment in fragments:
            index = record.seq.find(fragment)
            if index == -1:
                continue
            subrecord = record[index:index + len(fragment)]
            new_stickyend_record = StickyEndsSeqRecord(
                seq=fragment, features=subrecord.features,
                annotations=subrecord.annotations
            )
            record_fragments.append(new_stickyend_record)
        return record_fragments
