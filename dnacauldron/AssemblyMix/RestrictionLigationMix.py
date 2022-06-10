"""
"""

from copy import copy

from Bio import Restriction

from ..Fragment.StickyEndFragment import StickyEndFragment
from ..Filter import NoRestrictionSiteFilter
from .StickyEndAssemblyMix import StickyEndAssemblyMix
from ..biotools import autoselect_enzyme


class RestrictionLigationMix(StickyEndAssemblyMix):
    """Mix to represent and simulate restriction-ligation reactions.

    For reactions using type-2S enzymes (Golden Gate assembly etc.),
    see the helper function generate_type2s_restriction_mix().

    Parameters
    ----------

    parts
      List of Biopython records representing circular or blunt-end linear DNA
      fragments which will be restricted into sticky-end fragments and
      assembled toghether via sticky-end homologies.

    enzymes
      List of names of the enzyme to be used for the restriction 

    fragments
      A list of StickyEndFragment instances, which will assemble together based
      of sticky ends perfect homologies.

    fragment_filters
      List of functions of the form fragment=>True/False. If a fragment
      generates a "False" by at least one filter, it is taken out of the mix
      (this is used to remove unstable fragments that won't make it to the
      final assembly, for instance fragments with internal restriction sites).

    name
      Name of the mix as it will appear in reports

    annotate_fragments_with_parts
      If True, final constructs will have annotations "From xxx" indicating
      which part each sequence segment comes from.

    """

    def __init__(
        self,
        parts=None,
        enzymes=None,
        fragments=None,
        fragment_filters=(),
        name="restriction_mix",
        annotate_fragments_with_parts=True,
    ):
        # shallow copy seems sufficient and problem-free.
        # deepcopy would be safer but it is a computational bottleneck.
        self.parts = copy(parts) if parts else parts
        self.fragments = copy(fragments) if fragments else fragments
        if enzymes is not None:
            enzymes = [Restriction.__dict__[e] for e in enzymes]
        self.enzymes = enzymes
        self.fragment_filters = fragment_filters
        self.name = name
        self.annotate_fragments_with_parts = annotate_fragments_with_parts
        self.initialize()

    def compute_digest(self, part):
        """Compute the fragments resulting from the digestion"""
        return StickyEndFragment.list_from_record_digestion(part, self.enzymes)

    def compute_fragments(self):
        """Compute the (sticky-ended) fragments resulting from the digestion of
        the mix's parts by the mix's enzyme.

        Note that all fragments receive an annotation (feature) of type
        "source" that will show in the genbank of final constructs.
        """
        self.fragments = []
        for part in self.parts:
            for fragment in self.compute_digest(part):
                # The part is not a fragment if it hasn't been cut at all and
                # therefore doesn't have sticky ends. Exclude from fragments.
                if not hasattr(fragment.seq, "left_end"):
                    continue
                fragment.original_part = part
                self.annotate_fragment_with_part(fragment)
                self.fragments.append(fragment)


def generate_type2s_restriction_mix(parts, enzyme="auto", name="type2s_mix"):
    """Return a Type-2S compatible RestrictionLigationMix.
    
    If enzyme is "auto", it will be autoselected among BsmBI, BsaI, etc.
    """
    if enzyme == "auto":
        enzyme = autoselect_enzyme(parts)
    return RestrictionLigationMix(
        parts=parts,
        enzymes=[enzyme],
        fragment_filters=[NoRestrictionSiteFilter(str(enzyme))],
        name=name,
    )
