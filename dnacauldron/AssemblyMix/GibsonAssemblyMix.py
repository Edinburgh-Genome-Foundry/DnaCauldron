from copy import deepcopy


from .AssemblyMix import AssemblyMix
from ..StickyEndsSeq import StickyEndsSeqRecord

class GibsonAssemblyMix(AssemblyMix):
    """In construction. Do not use."""

    def __init__(self, constructs, min_homology=15, max_homology=200):

        self.constructs = deepcopy(constructs)
        self.min_homology = min_homology
        self.max_homology = max_homology
        self.initialize()

    def compute_fragments(self):
        self.fragments = list(self.constructs)

    @staticmethod
    def assemble(fragments, circularize=False, annotate_homologies=False):
        return StickyEndsSeqRecord.assemble(
            fragments,
            circularize=False,
            annotate_homologies=annotate_homologies
        )
