"""
"""

from ..biotools import set_record_topology
from . import mixins


class AssemblyMix(
    mixins.ConnectorsMixin,
    mixins.ConstructsMixin,
    mixins.FragmentsMixin,
    mixins.GraphsMixin,
    mixins.PlotsMixin,
):
    """Base class for assembly mixes.

    The subclasses (Type2sRestrictionMix and Gibson) implement
    their own version of how the original parts are broken into
    fragments, when two fragments will clip together, etc.
    """

    def initialize(self):
        """Precompute the fragments and connections graph of the mix."""
        if hasattr(self, 'parts') and self.parts is not None:
            for part in self.parts:
                set_record_topology(part, topology="default_to_linear")
            self.parts_dict = {cst.id: cst for cst in self.parts}
        if not hasattr(self, "fragments") or self.fragments is None:
            self.compute_fragments()
        self.compute_reverse_fragments()
        self.compute_connections_graph()
