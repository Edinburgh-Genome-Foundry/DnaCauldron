from Bio.Seq import Seq

try:
    # Biopython <1.78
    from Bio.Alphabet import DNAAlphabet

    has_dna_alphabet = True
except ImportError:
    # Biopython >=1.78
    has_dna_alphabet = False
from .StickyEnd import StickyEnd


class StickyEndSeq(Seq):
    """Represent sequences with sticky ends.

    It is used like a Biopython Seq, but with additional flanking sequences
    ``left_end`` and ``right_end``, which are ``StickyEnd`` objects.

    """

    def __init__(self, data, left_end=None, right_end=None, **k):
        Seq.__init__(self, str(data), **k)
        self.left_end = left_end
        self.right_end = right_end

    def reverse_complement(self, inplace=False):
        """The reverse is a StickyEndSeq with reversed ends

        left-right versions are interchanged and reverse complemented.
        """

        if has_dna_alphabet:  # Biopython <1.78
            sticky_end_seq = StickyEndSeq(
                str(Seq.reverse_complement(self)),
                left_end=None
                if self.right_end is None
                else self.right_end.reverse_complement(),
                right_end=None
                if self.left_end is None
                else self.left_end.reverse_complement(),
                alphabet=self.alphabet,
            )
        else:
            sticky_end_seq = StickyEndSeq(
                str(Seq(self).reverse_complement()),
                left_end=None
                if self.right_end is None
                else self.right_end.reverse_complement(),
                right_end=None
                if self.left_end is None
                else self.left_end.reverse_complement(),
            )

        return sticky_end_seq

    def will_clip_in_this_order_with(self, other):
        """Return whether this sequence will clip in this order with another.
        """
        return (self.right_end is not None) and self.right_end.will_clip_directly_with(
            other.left_end
        )

    def __repr__(self):
        content = Seq.__str__(self)
        if len(content) > 15:
            content = (
                content[:5].lower() + ("(%d)" % len(content)) + content[-5:].lower()
            )
        return "(%s-%s-%s)" % (repr(self.left_end), content, repr(self.right_end),)

    def __add__(self, other):
        assert self.will_clip_in_this_order_with(other)
        return StickyEndSeq(
            str(self) + str(self.right_end) + str(other),
            left_end=self.left_end,
            right_end=other.right_end,
        )

    def to_standard_sequence(self, discard_sticky_ends=False):
        """Return a string, basically left + middle + right"""
        if discard_sticky_ends:
            return Seq(str(self))
        else:
            left = str(self.left_end) if self.left_end else ""
            middle = self.to_standard_sequence(discard_sticky_ends=True)
            right = str(self.right_end) if self.right_end else ""
            return left + middle + right

    def slice_seq(self, start=None, end=None):
        """Slice the StickyEndSeq and return another instance.

        This function replaces the original sticky[start:end] approach.
        
        Parameters
        ----------

        start
            Start index (zero-based). Default None slices from start.
        
        end
            End index. Default None slices to the end.
        """
        # See Github issue #16 for details.
        std_seq = self.to_standard_sequence(discard_sticky_ends=True)
        sliced_seq = std_seq[start:end]  # biopython's slice works properly
        sliced_sticky = StickyEndSeq(str(sliced_seq), left_end=None, right_end=None)
        return sliced_sticky

    @staticmethod
    def list_from_sequence_digestion(sequence, enzyme, linear=True):
        """Compute the StickyEndSeqs resulting from a sequence digestion.

        This is one of the most important methods in DNA Cauldron. It uses
        Biopython to find restriction sites, and deals with the case where
        the sequence is circular, or the sequence is already sticky-ended.

        That was painful to write but hopefully the tests ensure it works well.

        Parameters
        ----------

        sequence
            Biopython Seq instance.

        enzyme
            Biopython Restriction enzyme instance or list of enzymes.

        linear
            True if the sequence is linear.
        """
        if isinstance(enzyme, (list, tuple)):
            if len(enzyme) == 1:
                enzyme = enzyme[0]
            else:
                enzyme, other_enzymes = enzyme[0], enzyme[1:]
                sticky_fragments = StickyEndSeq.list_from_sequence_digestion(
                    sequence=sequence, enzyme=other_enzymes, linear=linear
                )
                if not (sticky_fragments[0] is sequence):
                    linear = True
                return [
                    fragment
                    for sticky in sticky_fragments
                    for fragment in StickyEndSeq.list_from_sequence_digestion(
                        sequence=sticky, enzyme=enzyme, linear=linear
                    )
                ]
        n_cuts = len(enzyme.search(sequence, linear=linear))
        if n_cuts == 0:
            return [sequence]
        overhang = abs(enzyme.ovhg)
        right_end_sign = +1 if enzyme.is_3overhang() else -1

        if linear:
            fragments = enzyme.catalyse(sequence, linear=True)
        else:
            fragments = enzyme.catalyse(sequence + sequence, linear=True)
            fragments = fragments[1 : n_cuts + 1]

        if right_end_sign == -1:
            if not linear:
                # We check if StickyEndSeq instance was passed in the recursive loop
                # StickyEndSeq is a Seq class, so we use hasttr instead:
                if hasattr(fragments[0], "left_end"):  # Seq doesn't have this attr
                    overhang_bit = fragments[0].slice_seq(end=overhang)
                    new_fragment_seq = fragments[0].slice_seq(start=overhang)
                else:  # Seq class
                    overhang_bit = fragments[0][:overhang]
                    new_fragment_seq = fragments[0][overhang:]
                last_right_end = StickyEnd(overhang_bit, right_end_sign)
                first_left_end = StickyEnd(overhang_bit, -right_end_sign)
                sticky_fragments = [
                    StickyEndSeq(new_fragment_seq, left_end=first_left_end)
                ]

            else:
                sticky_fragments = [StickyEndSeq(fragments[0])]
            for f in fragments[1:]:
                if hasattr(f, "left_end"):
                    overhang_bit, new_fragment_seq = (
                        f.slice_seq(end=overhang),
                        f.slice_seq(start=overhang),
                    )
                else:  # Seq class
                    overhang_bit, new_fragment_seq = f[:overhang], f[overhang:]
                sticky_fragments[-1].right_end = StickyEnd(overhang_bit, right_end_sign)
                new_fragment = StickyEndSeq(
                    new_fragment_seq, left_end=StickyEnd(overhang_bit, -right_end_sign),
                )
                sticky_fragments.append(new_fragment)
            if not linear:
                sticky_fragments[-1].right_end = last_right_end
        else:
            left_end = None
            for f in fragments[:-1]:
                if hasattr(f, "left_end"):
                    overhang_bit, new_fragment_seq = (
                        f.slice_seq(end=-overhang),
                        f.slice_seq(start=-overhang),
                    )
                else:  # Seq class
                    new_fragment_seq, overhang_bit = f[:-overhang], f[-overhang:]
                right_end = StickyEnd(overhang_bit, right_end_sign)
                new_fragment = StickyEndSeq(
                    new_fragment_seq, left_end=left_end, right_end=right_end
                )
                sticky_fragments.append(new_fragment)
                left_end = StickyEnd(overhang_bit, -right_end_sign)
            if not linear:
                if hasattr(fragments[-1], "left_end"):
                    overhang_bit, new_fragment_seq = (
                        fragments[-1].slice_seq(end=-overhang),
                        fragments[-1].slice_seq(start=-overhang),
                    )
                else:  # StickyEndSeq instance passed (see above)
                    new_fragment_seq = fragments[-1][:-overhang]
                    overhang_bit = fragments[-1][-overhang:]
                first_left_end = StickyEnd(overhang_bit, -right_end_sign)
                last_right_end = StickyEnd(overhang_bit, right_end_sign)
                sticky_fragments[0].left_end = first_left_end
                sticky_fragments = [
                    StickyEndSeq(
                        new_fragment_seq, left_end=left_end, right_end=last_right_end,
                    )
                ]
            else:
                sticky_fragments.append(StickyEndSeq(fragments[-1], left_end=left_end))
        if hasattr(sequence, "left_end") and sticky_fragments[0].left_end is None:
            sticky_fragments[0].left_end = sequence.left_end
        if hasattr(sequence, "right_end") and sticky_fragments[-1].right_end is None:
            sticky_fragments[-1].right_end = sequence.right_end
        return sticky_fragments

    def ends_tuple(self):
        return (str(self.left_end), str(self.right_end))
