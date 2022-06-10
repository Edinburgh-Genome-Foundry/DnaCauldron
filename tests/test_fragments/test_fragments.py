import pytest
from Bio.Seq import Seq
from dnacauldron.Fragment import (
    StickyEnd,
    StickyEndSeq,
    StickyEndFragment,
)


def test_StickyEnd():
    sticky_end = StickyEnd(Seq("ATGC"), strand=1)
    assert sticky_end.__repr__() == "ATGC(+)"


def test_StickyEndSeq():
    sticky = StickyEndSeq(
        Seq("AAA"),
        left_end=StickyEnd("ATCG", strand=+1),
        # RC of left end so that it self-anneals:
        right_end=StickyEnd("ATCG", strand=-1),
    )
    assert sticky.__repr__() == "(ATCG(+)-AAA-ATCG(-))"

    # Longer than 15 bp:
    sticky = StickyEndSeq(
        Seq("AAAATTTTCCCCGGGG"),
        left_end=StickyEnd("ATCG", strand=+1),
        right_end=StickyEnd("ATCG", strand=-1),
    )
    assert sticky.__repr__() == "(ATCG(+)-aaaat(16)cgggg-ATCG(-))"

    # Test slice_seq()
    sticky = StickyEndSeq(
        Seq("AAAATTTTCCCCGGGG"),
        left_end=StickyEnd("ATCG", strand=+1),
        right_end=StickyEnd("ATCG", strand=-1),
    )
    sliced_sticky = sticky.slice_seq(start=0, end=4)
    assert sliced_sticky.__repr__() == "(None-AAAA-None)"

    sliced_sticky = sticky.slice_seq(start=12)
    assert sliced_sticky.__repr__() == "(None-GGGG-None)"

    sliced_sticky = sticky.slice_seq(end=12)
    assert sliced_sticky.__repr__() == "(None-AAAATTTTCCCC-None)"


def test_StickyEndSeqFragment():
    sticky = StickyEndSeq(
        Seq("TTT"),
        left_end=StickyEnd("AAAA", strand=+1),
        # Incompatible overhang:
        right_end=StickyEnd("ATCG", strand=-1),
    )
    sticky_fragment = StickyEndFragment(sticky)
    with pytest.raises(ValueError):
        sticky_fragment.circularized()
