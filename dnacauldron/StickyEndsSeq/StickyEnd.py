from Bio.Seq import Seq

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
        Seq.__init__(self, str(data).upper(), **k)
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
