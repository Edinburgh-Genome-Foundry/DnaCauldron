import networkx as nx
from ..AssemblyMixError import AssemblyMixError


class ConnectorsMixin:
    """Mixin for AssemblyMix"""
    
    def autoselect_connectors(self, connectors_records):
        """Select connectors necessary for circular assemblie(s) in this mix.

        This method assumes that the parts provided in the mix are the
        main parts of either a single or a combinatorial assembly (with
        well-defined slots), and the provided list of ``connector_records``
        contains "bridging" parts, some of which may be necessary for bridging
        the main parts of the assembly.

        The connectors record list contains records of parts (with or
        without backbone).

        If a solution is found, the method automatically adds the connectors
        to the mix, and returns the list of selected connectors records.

        Else, an exception is raised.
        """
        original_parts = self.parts
        all_part_ids = [c.id for c in original_parts]
        connectors_records = [
            c for c in connectors_records if c.id not in all_part_ids
        ]

        slotted_parts_records = [
            self.parts_dict[list(parts)[0]]
            for parts in self.compute_slots().values()
        ]
        self.parts = slotted_parts_records + connectors_records
        self.compute_fragments()
        self.initialize()
        graph = self.filtered_connections_graph
        components = sorted(
            nx.components.connected_components(graph.to_undirected()),
            key=lambda graph_: -len(graph_),
        )

        for component in components:

            newgraph = graph.copy()  # deepcopy(graph)
            newgraph.remove_nodes_from(
                set(newgraph.nodes()).difference(component)
            )
            all_paths = dict(nx.all_pairs_shortest_path(graph))
            parts_ids = set([rec.id for rec in slotted_parts_records])
            parts_nodes = [
                n
                for n in newgraph.nodes()
                if self.fragments_dict[n].original_part.id in parts_ids
            ]
            parts_graph = nx.DiGraph()
            parts_graph.add_edges_from(
                [
                    (node, other_node)
                    for node in parts_nodes
                    for other_node, path in all_paths[node].items()
                    if (other_node != node)
                    and (other_node in parts_nodes)
                    and len(set(path[1:-1]).intersection(set(parts_nodes)))
                    == 0
                ]
            )
            cycle = []
            if len(parts_graph) != len(original_parts):
                continue
            for cycle in nx.cycles.simple_cycles(parts_graph):

                if len(cycle) == len(parts_graph):
                    break
            if len(cycle) == len(parts_graph):
                break
        else:
            err = AssemblyMixError(
                message="No construct found involving all parts", mix=self
            )
            err.graph = graph
            raise err
        if len(cycle) == 0:
            raise ValueError("No solution found - a connector may be missing.")

        selected_connectors = [
            self.fragments_dict[n].original_part
            for (node1, node2) in zip(cycle, cycle[1:] + [cycle[0]])
            for n in all_paths[node1][node2][1:-1]
        ]

        # initialize the mix with the selected connectors
        self.parts = original_parts + selected_connectors
        self.compute_fragments()
        self.initialize()
        return selected_connectors
