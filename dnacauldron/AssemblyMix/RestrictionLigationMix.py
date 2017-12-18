"""
"""

from copy import copy

from Bio import Restriction



from ..tools import annotate_record
from ..StickyEndsSeq import StickyEndsSeqRecord
from .AssemblyMix import AssemblyMix
from .Filter import NoRestrictionSiteFilter

class RestrictionLigationMix(AssemblyMix):
    """Assembly mix for an enzymatic Restriction Ligation assembly.

    This includes modern assembly techniques such as Golden Gate as well as
    classical enzyme-based assembly.

    Parameters
    ----------

    constructs
      List of Biopython Seqrecords. Each seqrecord should have an attribute
      `linear` set to true or false (for circular constructs). It is advised to
      use method `load_genbank(filename, linear=True)` from `dnacauldron.tools`
      to load the constructs.

    enzyme
      Name of the ligation enzyme to use, e.g. 'BsmBI'
    """

    def __init__(self, constructs=None, enzyme=None, fragments=None,
                 fragments_filters='default'):
        # shallow copy seems sufficient and problem-free.
        # deepcopy would be safer but it is a computational bottleneck.
        self.constructs = copy(constructs) if constructs else constructs
        self.fragments = copy(fragments) if fragments else fragments
        self.enzyme = None if enzyme is None else Restriction.__dict__[enzyme]
        if fragments_filters == 'default':
            if enzyme is not None:
                fragments_filters = [NoRestrictionSiteFilter(str(self.enzyme))]
            else:
                fragments_filters = ()
        self.fragments_filters = fragments_filters
        self.initialize()

    def compute_digest(self, construct):
        """Compute the fragments resulting from the digestion"""
        return StickyEndsSeqRecord.list_from_record_digestion(
            construct, self.enzyme, linear=construct.linear)

    def compute_fragments(self):
        """Compute the (sticky-ended) fragments resulting from the digestion of
        the mix's constructs by the mix's enzyme.

        Note that all fragments receive an annotation (feature) of type
        "source" that will show in the genbank of final constructs.
        """
        self.fragments = []
        for construct in self.constructs:

            for fragment in self.compute_digest(construct):
                if not isinstance(fragment, StickyEndsSeqRecord):
                    continue
                fragment.original_construct = construct
                annotate_record(
                    fragment,
                    feature_type="misc_feature",
                    source=construct.name,
                    note="From " + construct.name,

                )
                self.fragments.append(fragment)

    @staticmethod
    def assemble(fragments, circularize=False, annotate_homologies=False):
        """Assemble sticky-end fragments into a single one (sticky or not).

        Parameters
        ----------

        fragments
          List of StickyEndsSeqRecord fragments

        circularize
          If True and if the two ends of the final assembly are compatible,
          circularize the construct, i.e. return a non-sticky record
          representing the circular assembly of the fragments.

        annotate_homologies
          If True, all homology regions that where formerly sticky ends will
          be annotated in the final record.
        """
        return StickyEndsSeqRecord.assemble(
            fragments,
            circularize=circularize,
            annotate_homologies=annotate_homologies
        )

    @staticmethod
    def will_clip_in_this_order(fragment1, fragment2):
        """Return True iff f1's right sticky end fits f2's left."""
        return fragment1.will_clip_in_this_order_with(fragment2)
