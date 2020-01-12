
from copy import copy
from Bio import Restriction
from .RestrictionLigationMix import RestrictionLigationMix
from ..StickyEndsSeq import StickyEndsSeqRecord
from ..biotools import autoselect_enzyme
from .Filter import NoRestrictionSiteFilter

class Type2sRestrictionMix(RestrictionLigationMix):
    """Assembly mix for an enzymatic Restriction Ligation assembly.

    This includes modern assembly techniques such as Golden Gate as well as
    classical enzyme-based assembly.

    Parameters
    ----------

    constructs
      List of Biopython Seqrecords. Each seqrecord should have an attribute
      `linear` set to true or false (for circular constructs). It is advised to
      use method `load_record(filename, linear=True)` from `dnacauldron.tools`
      to load the constructs.

    enzyme
      Name of the ligation enzyme to use, e.g. 'BsmBI'
    """

    def __init__(
        self,
        constructs=None,
        enzyme="auto",
        fragments=None,
        fragments_filters="default",
        name=None
    ):
        # shallow copy seems sufficient and problem-free.
        # deepcopy would be safer but it is a computational bottleneck.
        self.constructs = copy(constructs) if constructs else constructs
        self.fragments = copy(fragments) if fragments else fragments
        if enzyme == "auto":
            enzyme = autoselect_enzyme(constructs)
        self.enzyme_name = enzyme
        self.enzyme = None if enzyme is None else Restriction.__dict__[enzyme]
        self.enzymes = None if self.enzyme is None else [self.enzyme]
        if fragments_filters == "default":
            if enzyme is not None:
                fragments_filters = [NoRestrictionSiteFilter(str(self.enzyme))]
            else:
                fragments_filters = ()
        self.fragments_filters = fragments_filters
        self.name = None
        self.initialize()
    
    def compute_digest(self, construct):
        """Compute the fragments resulting from the digestion"""
        return StickyEndsSeqRecord.list_from_record_digestion(
            construct, self.enzyme
        )