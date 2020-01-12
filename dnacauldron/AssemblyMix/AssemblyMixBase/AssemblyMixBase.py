"""
"""

from ...biotools import set_record_topology
from .AssemblyMixConnectorsMixin import AssemblyMixConnectorsMixin
from .AssemblyMixConstructsMixin import AssemblyMixConstructsMixin
from .AssemblyMixFragmentsMixin import AssemblyMixFragmentsMixin
from .AssemblyMixGraphsMixin import AssemblyMixGraphsMixin
from .AssemblyMixPlotsMixin import AssemblyMixPlotsMixin

class AssemblyMixBase(
    AssemblyMixConnectorsMixin,
    AssemblyMixConstructsMixin,
    AssemblyMixFragmentsMixin,
    AssemblyMixGraphsMixin,
    AssemblyMixPlotsMixin
):
    """Base class for assembly mixes.

    The subclasses (Type2sRestrictionMix and GibsonAssemblyMix) implement
    their own version of how the original parts are broken into
    fragments, when two fragments will clip together, etc.
    """

    def initialize(self):
        """Precompute the fragments and connections graph of the mix."""
        if self.parts is not None:
            for part in self.parts:
                set_record_topology(
                    part, topology="linear", pass_if_already_set=True
                )
            self.parts_dict = {cst.id: cst for cst in self.parts}
        if not hasattr(self, "fragments") or self.fragments is None:
            self.compute_fragments()
        self.compute_reverse_fragments()
        self.compute_connections_graph()
