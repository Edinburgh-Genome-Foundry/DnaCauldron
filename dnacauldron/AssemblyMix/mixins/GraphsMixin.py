import itertools as itt
import networkx as nx
from ...biotools import reverse_complement


class GraphsMixin:
    """Mixin for AssemblyMix"""

    def compute_connections_graph(self):
        """Compute a graph where nodes are fragments and edges indicate
        which fragments can clip together.

        The graph (stored in self.connection_graph) is directed, an edge
        f1->f2 indicating that fragments f1 and f2 will clip in this order.
        """

        all_fragments = self.fragments + self.reverse_fragments
        items = list(enumerate(all_fragments))
        self.fragments_dict = dict(items)
        self.connections_graph = nx.DiGraph()
        for (i1, fragment1), (i2, fragment2) in itt.combinations(items, 2):
            if self.will_clip_in_this_order(fragment1, fragment2):
                self.connections_graph.add_edge(i1, i2)
            if self.will_clip_in_this_order(fragment2, fragment1):
                self.connections_graph.add_edge(i2, i1)

    def list_overhangs(self, filtered_fragments_only=False):
        """Return a list of overhangs in the mix.
    
        Warning: only overhangs on non-reversed fragments are returned, not
        their reverse-complement.

        Parameters
        ----------

        filtered_fragments_only
          If true, only overhangs from filtered fragments are returned
        """
        fragments = (
            self.filtered_fragments
            if filtered_fragments_only
            else self.fragments
        )
        overhangs = [
            str(end)
            for fragment in fragments
            if not fragment.is_reversed
            for end in (fragment.seq.left_end, fragment.seq.right_end)
        ]
        return sorted(set(overhangs))

    @property
    def filtered_connections_graph(self):
        """Networkx Graph of the filered fragments and how they assemble.

        The nodes of the graph are numbers such that N represents the fragment
        indexed by self.fragments_dict[N]

        """
        graph = nx.DiGraph(self.connections_graph)
        graph.remove_nodes_from(
            [
                node
                for node in graph.nodes()
                if not all(
                    [
                        fl(self.fragments_dict[node])
                        for fl in self.fragment_filters
                    ]
                )
            ]
        )
        return graph

    @property
    def uniquified_connection_graph(self, filtered=True):
        def transform(n):
            return self.fragments_dict[n].original_part.id

        def graph_hash(g):
            sorted_edges = tuple(sorted([tuple(sorted(e)) for e in g.edges]))
            signature = tuple(tuple(sorted(g.nodes)) + sorted_edges)
            return hash(signature)

        def graph_with_transformed_nodes(graph, transformation):
            transformations = {n: transformation(n) for n in graph}
            new_graph = graph.__class__(
                [
                    (transformations[n1], transformations[n2])
                    for n1, n2 in graph.edges()
                ]
            )
            for n in transformations.values():
                if n not in new_graph:
                    new_graph.add_node(n)
            return new_graph

        if filtered:
            graph = self.filtered_connections_graph
        else:
            graph = self.connections_graph
        components = list(nx.connected_components(graph.to_undirected()))
        nodes_by_hash = {}
        for nodes in components:
            component = graph.subgraph(nodes)
            component = graph_with_transformed_nodes(component, transform)
            component_hash = graph_hash(component)
            if component_hash in nodes_by_hash:
                continue
            nodes_by_hash[component_hash] = nodes
        kept_nodes = set().union(*nodes_by_hash.values())
        return graph.subgraph(kept_nodes)

