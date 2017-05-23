from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import DNAAlphabet

class StickyEnd(Seq):
    """A class to represent the sticky end of a sequence.

    It is used exactly like a Biopython sequence.

    Parameters
    ----------

    data
      A DNA sequence in ATGC format

    strand
      The strand (+1 or -1) on which the protusion is

    **k
      Optional keyword arguments for the sequence, such as ``alphabet`` etc.
    """

    def __init__(self, data, strand, **k):
        Seq.__init__(self, str(data), **k)
        self.strand = strand

    def reverse_complement(self):
        return StickyEnd(
            str(Seq.reverse_complement(self)),
            strand=-self.strand,
            alphabet=self.alphabet
        )

    def __repr__(self):
        return "%s(%s)" % (Seq.__str__(self),
                           {1: "+", -1: "-"}[self.strand])

    def will_clip_directly_with(self, other):
        return ((other is not None) and
                (len(self) > 0) and
                (self.strand == -other.strand) and
                (str(self) == str(other)))


class StickyEndsSeq(Seq):
    """Represent sequences with sticky ends.

    It is used like a Biopython Seq, but with additional flanking sequences
    ``left_end`` and ``right_end``, which are ``StickyEnd`` objects.

    """

    def __init__(self, data, left_end=None, right_end=None, **k):
        Seq.__init__(self, str(data), **k)
        self.left_end = left_end
        self.right_end = right_end

    def reverse_complement(self):
        return StickyEndsSeq(
            str(Seq.reverse_complement(self)),
            left_end=None if self.right_end is None else
            self.right_end.reverse_complement(),
            right_end=None if self.left_end is None else
            self.left_end.reverse_complement(),
            alphabet=self.alphabet
        )

    def will_clip_in_this_order_with(self, other):
        return ((self.right_end is not None) and
                self.right_end.will_clip_directly_with(other.left_end))

    def circularized(self):
        if not self.will_clip_in_this_order_with(self):
            raise ValueError("Only constructs with two compatible sticky ends"
                             " can be circularized")
        result = Seq(str(self.left_end)) + self
        result.linear = False
        return result

    def __repr__(self):
        content = Seq.__str__(self)
        if len(content) > 15:
            content = (content[:5].lower() +
                       ("(%d)" % len(content)) +
                       content[-5:].lower())
        return "(%s-%s-%s)" % (repr(self.left_end),
                               content,
                               repr(self.right_end))

    def __add__(self, other):
        assert self.will_clip_in_this_order_with(other)
        return StickyEndsSeq(
            str(self) + str(self.right_end) + str(other),
            left_end=self.left_end,
            right_end=other.right_end
        )

    def to_standard_sequence(self, discard_sticky_ends=False):
        if discard_sticky_ends:
            return Seq(str(self))


class StickyEndsSeqRecord(SeqRecord):
    """Biopython SeqRecord whose sequence has sticky ends."""

    def will_clip_in_this_order_with(self, other):
        right_end = self.seq.right_end
        return ((right_end is not None) and
                right_end.will_clip_directly_with(other.seq.left_end))

    def circularized(self, annotate_homology=False, annotation_type="Feature",
                     qualifiers=None):
        if not self.will_clip_in_this_order_with(self):
            raise ValueError("Only constructs with two compatible sticky ends"
                             " can be circularized")
        connector = SeqRecord(Seq(str(self.seq.left_end)))
        if annotate_homology:
            connector.features = [
                SeqFeature(FeatureLocation(0, len(connector), 1),
                           type=annotation_type,
                           qualifiers={"label": "homology"})
            ]
        result = connector + self
        result.linear = False
        return result

    @staticmethod
    def assemble(fragments, circularize=False, annotate_homologies=False):
        result = fragments[0]
        for fragment in fragments[1:]:
            result = result.assemble_with(
                fragment,
                annotate_homology=annotate_homologies
            )
        # result = sum(fragments[1:], fragments[0])
        if circularize:
            result = result.circularized(annotate_homology=annotate_homologies)
        result.seq.alphabet = DNAAlphabet()
        return result

    def assemble_with(self, other, annotate_homology=False,
                      annotation_type="misc_feature", **qualifiers):
        connector = SeqRecord(Seq(str(self.seq.right_end)))
        if len(qualifiers) == 0:
            qualifiers = {"label": "homology"}
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


def digest_sequence_with_sticky_ends(sequence, enzyme, linear=True):
    n_cuts = len(enzyme.search(sequence, linear=linear))
    if n_cuts == 0:
        return sequence
    overhang = abs(enzyme.ovhg)
    right_end_sign = +1 if enzyme.is_3overhang() else -1
    fragments = enzyme.catalyse(sequence, linear=linear)
    if linear:
        fragments = enzyme.catalyse(sequence, linear=True)
    else:
        fragments = enzyme.catalyse(
            sequence + sequence, linear=True)[1:n_cuts + 1]
    if right_end_sign == -1:
        if not linear:
            overhang_bit = fragments[0][:overhang]
            new_fragment_seq = fragments[0][overhang:]
            last_right_end = StickyEnd(overhang_bit, right_end_sign)
            first_left_end = StickyEnd(overhang_bit, -right_end_sign)
            sticky_fragments = [StickyEndsSeq(new_fragment_seq,
                                              left_end=first_left_end)]

        else:
            sticky_fragments = [StickyEndsSeq(fragments[0])]
        for f in fragments[1:]:
            overhang_bit, new_fragment_seq = f[:overhang], f[overhang:]
            sticky_fragments[-1].right_end = StickyEnd(
                overhang_bit, right_end_sign)
            new_fragment = StickyEndsSeq(new_fragment_seq,
                                         left_end=StickyEnd(overhang_bit,
                                                            -right_end_sign))
            sticky_fragments.append(new_fragment)
        if not linear:
            sticky_fragments[-1].right_end = last_right_end
    else:
        left_end = None
        for f in fragments[:-1]:
            new_fragment_seq, overhang_bit = f[:-overhang], f[-overhang:]
            right_end = StickyEnd(overhang_bit, right_end_sign)
            new_fragment = StickyEndsSeqRecord(new_fragment_seq,
                                               left_end=left_end,
                                               right_end=right_end)
            sticky_fragments.append(new_fragment)
            left_end = StickyEnd(overhang_bit, -right_end_sign)
        if not linear:
            new_fragment_seq = fragments[-1][:-overhang]
            overhang_bit = fragments[-1][-overhang:]

            first_left_end = StickyEnd(overhang_bit, -right_end_sign)
            last_right_end = StickyEnd(overhang_bit, right_end_sign)
            sticky_fragments[0].left_end = first_left_end
            sticky_fragments = [StickyEndsSeq(new_fragment_seq,
                                              left_end=left_end,
                                              right_end=last_right_end
                                              )]
        else:
            sticky_fragments.append(StickyEndsSeqRecord(fragments[-1],
                                                        left_end=left_end))
    return sticky_fragments


def digest_seqrecord_with_sticky_ends(seqrecord, enzyme, linear=True):
    n_cuts = len(enzyme.search(seqrecord.seq, linear=linear))
    if n_cuts == 0:
        return [seqrecord]
    if not linear:
        record_fragments = digest_seqrecord_with_sticky_ends(
            seqrecord + seqrecord,
            enzyme,
            linear=True
        )[1:n_cuts + 1]
        return record_fragments
    fragments = digest_sequence_with_sticky_ends(
        seqrecord.seq, enzyme, linear=linear)
    record_fragments = []
    for fragment in fragments:
        index = seqrecord.seq.find(fragment)
        if index == -1:
            continue
        subrecord = seqrecord[index:index + len(fragment)]
        new_stickyend_record = StickyEndsSeqRecord(
            seq=fragment, features=subrecord.features,
            annotations=subrecord.annotations
        )
        record_fragments.append(new_stickyend_record)
    return record_fragments
