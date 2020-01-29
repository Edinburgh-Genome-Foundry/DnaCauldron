from ...Fragment.StickyEndFragment import StickyEndFragment
from ..AssemblyMix import AssemblyMix
from .SlotsMixin import SlotsMixin
from .PlotsMixin import PlotsMixin

class StickyEndAssemblyMix(AssemblyMix, SlotsMixin, PlotsMixin):
    """Class to represent and simulate the assembly of sticky-ended fragments.

    Used notably for BASIC assembly and other places where we are mixing
    sticky-ended fragments, not necessarily from enzyme digestions (see the
    subclass RestrictionLigationMix for this case).

    Parameters
    ----------

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
    

    """
    
    def __init__(self, fragments, fragment_filters=(), name='sticky_ends_mix'):
        self.fragments = fragments
        self.fragment_filters = fragment_filters
        self.name = name
        self.initialize()
    
    @staticmethod
    def assemble(fragments, circularize=False, annotate_homologies=False):
        """Assemble sticky-end fragments into a single one (sticky or not).

        Parameters
        ----------

        fragments
          List of StickyEndFragment fragments

        circularize
          If True and if the two ends of the final assembly are compatible,
          circularize the construct, i.e. return a non-sticky record
          representing the circular assembly of the fragments.

        annotate_homologies
          If True, all homology regions that where formerly sticky ends will
          be annotated in the final record.
        """
        return StickyEndFragment.assemble(
            fragments,
            circularize=circularize,
            annotate_homologies=annotate_homologies,
        )

    @staticmethod
    def will_clip_in_this_order(fragment1, fragment2):
        """Return True iff f1's right sticky end fits f2's left."""
        return fragment1.will_clip_in_this_order_with(fragment2)
