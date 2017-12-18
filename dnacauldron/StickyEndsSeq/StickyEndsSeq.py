from Bio.Seq import Seq
from .StickyEnd import StickyEnd

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
        else:
            left = str(self.left_end) if self.left_end else ''
            middle = self.to_standard_sequence(discard_sticky_ends=True)
            right = str(self.right_end) if self.right_end else ''
            return left + middle + right

    @staticmethod
    def list_from_sequence_digestion(sequence, enzyme, linear=True):
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
                new_fragment = StickyEndsSeq(new_fragment_seq,
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
                sticky_fragments.append(StickyEndsSeq(fragments[-1],
                                                      left_end=left_end))
        return sticky_fragments
