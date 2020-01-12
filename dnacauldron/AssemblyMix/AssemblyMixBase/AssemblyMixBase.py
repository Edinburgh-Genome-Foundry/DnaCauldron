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
    their own version of how the original constructs are broken into
    fragments, when two fragments will clip together, etc.
    """

    def initialize(self):
        """Precompute the fragments and connections graph of the mix."""
        if self.constructs is not None:
            for construct in self.constructs:
                set_record_topology(
                    construct, topology="linear", pass_if_already_set=True
                )
            self.constructs_dict = {cst.id: cst for cst in self.constructs}
        if not hasattr(self, "fragments") or self.fragments is None:
            self.compute_fragments()
        self.compute_reverse_fragments()
        self.compute_connections_graph()
