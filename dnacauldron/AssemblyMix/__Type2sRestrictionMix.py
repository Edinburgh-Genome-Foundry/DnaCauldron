
from copy import copy
from Bio import Restriction
from ..Fragments.StickyEndFragment import StickyEndFragment
from ..biotools import autoselect_enzyme
from .RestrictionLigationMix import RestrictionLigationMix
from .Filter import NoRestrictionSiteFilter

class Type2sRestrictionMix(RestrictionLigationMix):
    """Assembly mix for an enzymatic Restriction Ligation assembly.

    This includes modern assembly techniques such as Golden Gate as well as
    classical enzyme-based assembly.

    Parameters
    ----------

    parts
      List of Biopython Seqrecords. Each seqrecord should have an attribute
      `record.annotations['topology]` set to 'circular' or 'linear'.

    enzyme
      Name of the ligation enzyme to use, e.g. 'BsmBI'
    """

    def __init__(
        self,
        parts=None,
        enzyme="auto",
        fragments=None,
        fragments_filters="default",
        name="type2s_mix"
    ):
        # shallow copy seems sufficient and problem-free.
        # deepcopy would be safer but it is a computational bottleneck.
        self.parts = copy(parts) if parts else parts
        self.fragments = copy(fragments) if fragments else fragments
        if enzyme == "auto":
            enzyme = autoselect_enzyme(parts)
        self.enzyme_name = enzyme
        self.enzyme = None if enzyme is None else Restriction.__dict__[enzyme]
        self.enzymes = None if self.enzyme is None else [self.enzyme]
        if fragments_filters == "default":
            if enzyme is not None:
                fragments_filters = [NoRestrictionSiteFilter(str(self.enzyme))]
            else:
                fragments_filters = ()
        self.fragments_filters = fragments_filters
        self.name = name
        self.initialize()