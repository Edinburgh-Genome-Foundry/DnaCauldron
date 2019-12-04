import itertools as itt
import networkx as nx
from dnacauldron.tools import reverse_complement


class AssemblyMixGraphsMixin:
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
        return sorted(
            set(
                [
                    str(end)
                    for fragment in fragments
                    if not fragment.is_reverse
                    for end in (fragment.seq.left_end, fragment.seq.right_end)
                ]
            )
        )

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
                        for fl in self.fragments_filters
                    ]
                )
            ]
        )
        return graph

    @property
    def uniquified_connection_graph(self, filtered=True):
        def transform(n):
            return self.fragments_dict[n].original_construct.id

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

    def compute_slots(self):
        """Return a dict {standardized_slot: set(fragments_id)}.

        If a fragment has left and right sticky ends (o1, o2), the
        standardized version will be will be whichever is alphabetically
        smaller between  (o1, o2) and (revcomp(o2), revcomp(o1))

        For instance

        """

        def slot(fragment):
            def std_sticky(s):
                return "" if s is None else str(s)

            dir_slot = l, r = (
                std_sticky(fragment.seq.left_end),
                std_sticky(fragment.seq.right_end),
            )
            rev_slot = (reverse_complement(r), reverse_complement(l))
            return min(dir_slot, rev_slot)

        slots = {}
        for f in self.filtered_fragments:
            f_slot = slot(f)
            if f_slot not in slots:
                slots[f_slot] = set()
            slots[f_slot].add(f.original_construct.id)
        return slots

    def slots_graph(self, with_overhangs=True, directed=True):
        """Compute the slots graph of the graph.

        In this graph, a node represents a slot (i.e. left-right overhangs)
        and edges represent slots sharing one overhang. If with_overhangs is
        True, additional nodes are added to represent the overhangs (this
        makes the graph more informative).
        """

        def std_overhang(o):
            return min(o, reverse_complement(o))

        slots_list = list(self.compute_slots().keys())
        edges = []
        if with_overhangs:
            # THE CODE BLOCK BELOW IS COMPLICATED AND DOESNT BRING MUCH,
            # JUST ARROW DIRECTIONS WITH LITTLE INFORMATIVE VALUE IN THE
            # SLOTS GRAPH
            if directed:

                def slot_to_edges(slot):
                    left, right = slot
                    std_left, std_right = [std_overhang(o) for o in slot]
                    edges = []
                    if left != "":
                        edge = (std_left, slot)
                        if std_left != left:
                            edge = edge[::-1]
                        edges.append(edge)
                    if right != "":
                        edge = (slot, std_right)
                        if std_left != left:
                            edge = edge[::-1]
                        edges.append(edge)
                    return edges

                edges = [
                    edge for slot in slots_list for edge in slot_to_edges(slot)
                ]
            else:
                edges = [
                    (slot, std_overhang(e))
                    for slot in slots_list
                    for e in slot
                    if e != ""
                ]
        else:

            def slots_have_common_overhang(s1, s2):
                """Return True iff the slots have 1+ common overhang."""
                s2 = [std_overhang(o) for o in s2]
                return any([std_overhang(o) in s2 for o in s1])

            edges = [
                (s1, s2)
                for (s1, s2) in itt.combinations(slots_list, 2)
                if slots_have_common_overhang(s1, s2)
            ]
        if with_overhangs and directed:
            return nx.DiGraph(edges)
        else:
            return nx.Graph(edges)
