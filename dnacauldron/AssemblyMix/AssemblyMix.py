"""
"""

import itertools as itt
from copy import deepcopy

import networkx as nx
import numpy as np

from ..tools import reverse_complement
from .FragmentsChain import FragmentsChain

class AssemblyError(Exception):
    pass

class AssemblyMix:
    """General class for assembly mixes.

    The subclasses (RestrictionLigationMix and GibsonAssemblyMix) implement
    their own version of how the original constructs are broken into
    fragments, when two fragments will clip together, etc.
    """

    def initialize(self):
        """Precompute the fragments and connections graph of the mix."""
        if self.constructs is not None:
            for construct in self.constructs:
                if not hasattr(construct, "linear"):
                    construct.linear = True  # assumed linear by default
            self.constructs_dict = {
                cst.id: cst
                for cst in self.constructs
            }
        if not hasattr(self, "fragments") or self.fragments is None:
            self.compute_fragments()
        self.compute_reverse_fragments()
        self.compute_connections_graph()

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
        fragments = (self.filtered_fragments if filtered_fragments_only
                     else self.fragments)
        return sorted(set([
            str(end)
            for fragment in fragments
            if not fragment.is_reverse
            for end in (fragment.seq.left_end, fragment.seq.right_end)
        ]))
    @property
    def filtered_connections_graph(self):
        """Networkx Graph of the filered fragments and how they assemble.

        The nodes of the graph are numbers such that N represents the fragment
        indexed by self.fragments_dict[N]

        """
        graph = nx.DiGraph(self.connections_graph)
        graph.remove_nodes_from([
            node for node in graph.nodes()
            if not all([fl(self.fragments_dict[node])
                        for fl in self.fragments_filters])
        ])
        return graph

    def compute_slots(self):
        """Return a dict {standardized_slot: set(fragments_id)}.

        If a fragment has left and right sticky ends (o1, o2), the
        standardized version will be will be whichever is alphabetically
        smaller between  (o1, o2) and (revcomp(o2), revcomp(o1))

        For instance

        """
        def slot(fragment):
            def std_sticky(s):
                return '' if s is None else str(s)
            dir_slot = l, r = (std_sticky(fragment.seq.left_end),
                               std_sticky(fragment.seq.right_end))
            rev_slot = (reverse_complement(r), reverse_complement(l))
            return min(dir_slot, rev_slot)
        slots = {}
        for f in self.filtered_fragments:
            f_slot = slot(f)
            if f_slot not in slots:
                slots[f_slot] = set()
            slots[f_slot].add(f.original_construct.id)
        return slots

    def slots_graph(self, with_overhangs=True):
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
            edges = [
                (slot, std_overhang(e))
                for slot in slots_list
                for e in slot
                if e != ''
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
        return nx.Graph(edges)


    def compute_reverse_fragments(self):
        """Precompute self.reverse_fragments.

        This method also marks all "direct" fragments in the mix as
        `fragment.is_reverse=True` and all "reverse" fragments as
        `fragment.is_reverse=False`.
        """
        self.reverse_fragments = []
        for fragment in self.fragments:
            fragment.is_reverse = False
            new_fragment = fragment.reverse_complement()
            new_fragment.is_reverse = True
            new_fragment.reverse_fragment = fragment
            fragment.reverse_fragment = new_fragment
            new_fragment.original_construct = fragment.original_construct
            self.reverse_fragments.append(new_fragment)

    def compute_random_circular_fragments_sets(self, staling_cutoff=100,
                                               fragments_sets_filters=()):
        """Return an iterator over all the lists of fragments [f1, f2, f3...fn]
        that can assemble into a circular construct.

        This means that fragment f1 will clip with f2 (in this order),
        f2 with f3... and fn with f1.

        Parameters
        ----------

        fragments_sets_filters
          A list of test functions of the form (set->True/False) where "set" is
          a list of fragments from the mix (which assemble into a circular
          construct). Only circular fragments sets passing all these tests
          will be returned

        randomize
          If set to False, the circular fragments sets will be returned one by
          one until the last one, in an order implemented by
          networkx.simple_cycles.
          True, the circular sets returned will be drawn randomly
          (a circular set will only be returned once). This is very practical
          to obtain a sample out of a combinatorial assembly mix. However this
          feature is a bit experimental, and the iteration will certainly stale
          before all cycles have been found, because the randomizer can't find
          any new cycle.

        randomization_staling_cutoff
          If randomize is True, the randomizer will throw an error if the
          latest C cycles it has drawn had already been seen before, where C
          is the randomization staling cutoff.

        """

        def generator():
            """Return random cycles from the connections graph.

            The randomness is introduced by permuting the nodes names,
            running `networkx.circular_paths` once, permuting the nodes
            names again, etc.
            """

            graph = self.filtered_connections_graph
            seen_hashes = set()
            graph_nodes = list(graph.nodes())
            node_to_index = {node: i for i, node in enumerate(graph_nodes)}
            while True:
                permutation = np.arange(len(graph_nodes))
                np.random.shuffle(permutation)
                antipermutation = np.argsort(permutation)
                new_graph = nx.DiGraph([
                    (permutation[node_to_index[node1]],
                     permutation[node_to_index[node2]])
                    for node1, node2 in graph.edges()

                ])
                counter = 0
                for cycle in nx.simple_cycles(new_graph):
                    cycle = [antipermutation[i] for i in cycle]
                    fragments = [self.fragments_dict[graph_nodes[i]]
                                 for i in cycle]
                    cycle = FragmentsChain(fragments,
                                           is_cycle=True).standardized()
                    cycle_hash = hash(cycle)
                    if cycle_hash in seen_hashes:
                        counter += 1
                        if counter > staling_cutoff:
                            raise ValueError(
                                "Randomization staled. Only randomize when"
                                " the search space is huge."
                            )
                        continue
                    seen_hashes.add(cycle_hash)
                    if all(fl(cycle.fragments)
                           for fl in fragments_sets_filters):
                        yield cycle.fragments
                        break
                else:
                    break

        return generator()

    def compute_circular_fragments_sets(self, fragments_sets_filters=()):
        """Return an iterator over all the lists of fragments [f1, f2, f3...fn]
        that can assemble into a circular construct.

        This means that fragment f1 will clip with f2 (in this order),
        f2 with f3... and fn with f1.

        Parameters
        ----------

        fragments_sets_filters
          A list of test functions of the form (set->True/False) where "set" is
          a list of fragments from the mix (which assemble into a circular
          construct). Only circular fragments sets passing all these tests
          will be returned
        """


        def generator():
            """Iterate over all circular paths in the connexion graph
            using Networkx's `simple_paths`."""
            seen_hashes = set()
            for cycle in nx.simple_cycles(self.filtered_connections_graph):
                cycle = [self.fragments_dict[i] for i in cycle]
                cycle = FragmentsChain(cycle, is_cycle=True).standardized()
                cycle_hash = hash(cycle)
                if cycle_hash in seen_hashes:
                    continue
                seen_hashes.add(cycle_hash)
                if all(fl(cycle.fragments)
                       for fl in fragments_sets_filters):
                    yield cycle.fragments

        return generator()

    def compute_circular_assemblies(self, fragments_sets_filters=(),
                                    seqrecord_filters=(),
                                    annotate_homologies=False,
                                    randomize=False,
                                    randomization_staling_cutoff=100):
        """Return a generator listing the circular assemblies in the graph.

        Parameters
        ----------

        fragments_sets_filters
          A list of test functions of the form (set->True/False) where "set" is
          a list of fragments from the mix (which assemble into a circular
          construct). Only circular fragments sets passing all these tests
          will be returned

        seqrecord_filters
          A list of test functions of the form (record->True/False) where
          "record" is the biopython record of a circular assembly found.
          Only records passing all these tests will be returned

        annotate_homologies
          If True, the junctions between assembled fragments will be annotated
          in the final record with a feature of type 'homology' and label
          equal to the homology (if <8bp), else simply 'homology'.

        randomize
          If set to False, the circular fragments sets will be returned one by
          one until the last one, in an order implemented by
          networkx.simple_cycles.
          True, the circular sets returned will be drawn randomly
          (a circular set will only be returned once). This is very practical
          to obtain a sample out of a combinatorial assembly mix. However this
          feature is a bit experimental, and the iteration will certainly stale
          before all cycles have been found, because the randomizer can't find
          any new cycle.

        randomization_staling_cutoff
          If randomize is True, the randomizer will throw an error if the
          latest C cycles it has drawn had already been seen before, where C
          is the randomization staling cutoff.

        """

        def assemblies_generator():
            if randomize:
                fragments = self.compute_random_circular_fragments_sets(
                    fragments_sets_filters=fragments_sets_filters,
                    staling_cutoff=randomization_staling_cutoff
                )
            else:
                fragments = self.compute_circular_fragments_sets(
                    fragments_sets_filters=fragments_sets_filters)
            for fragments in fragments:
                construct = self.assemble(
                    fragments,
                    circularize=True,
                    annotate_homologies=annotate_homologies
                )
                if all(fl(construct) for fl in seqrecord_filters):
                    construct.fragments = fragments
                    yield construct
        return assemblies_generator()

    def compute_linear_assemblies(self,
                                  fragments_sets_filters=(),
                                  min_parts=2,
                                  seqrecord_filters=(),
                                  annotate_homologies=False):
        """Return a generator listing the possible linear assemblies.

        Parameters
        ----------

        fragments_sets_filters
          A list of test functions of the form (set->True/False) where "set" is
          a list of fragments from the mix (which assemble into a circular
          construct). Only circular fragments sets passing all these tests
          will be returned

        min_parts
          Assemblies with less than this number of parts will be ignored.

        seqrecord_filters
          A list of test functions of the form (record->True/False) where
          "record" is the biopython record of a circular assembly found.
          Only records passing all these tests will be returned

        annotate_homologies
          If True, the junctions between assembled fragments will be annotated
          in the final record with a feature of type 'homology' and label
          equal to the homology (if <8bp), else simply 'homology'.

        randomize
          If set to False, the circular fragments sets will be returned one by
          one until the last one, in an order implemented by
          networkx.simple_cycles.
          True, the circular sets returned will be drawn randomly
          (a circular set will only be returned once). This is very practical
          to obtain a sample out of a combinatorial assembly mix. However this
          feature is a bit experimental, and the iteration will certainly stale
          before all cycles have been found, because the randomizer can't find
          any new cycle.

        randomization_staling_cutoff
          If randomize is True, the randomizer will throw an error if the
          latest C cycles it has drawn had already been seen before, where C
          is the randomization staling cutoff.

        Notes
        ------

        This is a bit undertested as there have been little use cases.

        """
        seen_hashes = set()
        g = self.filtered_connections_graph
        for source, targets in nx.shortest_path(g).items():
            for target, path in targets.items():
                if (len(path) < min_parts):
                    continue
                fragments = [self.fragments_dict[f] for f in path]
                if not all([fl(fragments) for fl in fragments_sets_filters]):
                    continue
                chain = FragmentsChain(fragments).standardized()
                chain_hash = hash(chain)
                if chain_hash in seen_hashes:
                    continue
                seen_hashes.add(chain_hash)
                fragments_assembly = self.assemble(
                    fragments, annotate_homologies=annotate_homologies)
                if all([fl(fragments_assembly) for fl in seqrecord_filters]):
                    yield(fragments_assembly)

    @property
    def filtered_fragments(self):
        """Return the fragments of the mix passing all the tests

        Generally used to remove fragments containing a restriction site used
        in a Type2S assembly.
        """
        return [
            f for f in (self.fragments + self.reverse_fragments)
            if all([fl(f) for fl in self.fragments_filters])
        ]

    def autoselect_connectors(self, connectors_records):
        """Select connectors necessary for circular assemblie(s) in this mix.

        This method assumes that the constructs provided in the mix are the
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
        original_constructs = self.constructs
        all_construct_ids = [c.id for c in original_constructs]
        connectors_records = [c for c in connectors_records
                              if c.id not in all_construct_ids]

        slotted_parts_records = [
             self.constructs_dict[list(parts)[0]]
             for parts in self.compute_slots().values()
        ]
        self.constructs = slotted_parts_records + connectors_records
        self.compute_fragments()
        self.initialize()
        graph = self.filtered_connections_graph
        components = sorted(
            nx.components.connected_component_subgraphs(graph.to_undirected()),
            key=lambda graph_: -len(graph_)
        )

        for component in components:

            newgraph = graph.copy()  # deepcopy(graph)
            newgraph.remove_nodes_from(
                set(newgraph.nodes()).difference(component.nodes())
            )
            all_paths = dict(nx.all_pairs_shortest_path(graph))
            parts_ids = set([rec.id for rec in slotted_parts_records])
            parts_nodes = [
                n for n in newgraph.nodes()
                if self.fragments_dict[n].original_construct.id in parts_ids
            ]
            parts_graph = nx.DiGraph()
            parts_graph.add_edges_from([
                (node, other_node)
                for node in parts_nodes
                for other_node, path in all_paths[node].items()
                if (other_node != node)
                and (other_node in parts_nodes)
                and len(set(path[1: -1]).intersection(set(parts_nodes))) == 0
            ])
            cycle = []
            if (len(parts_graph) != len(original_constructs)):
                continue
            for cycle in nx.cycles.simple_cycles(parts_graph):

                if len(cycle) == len(parts_graph):
                    break
            if (len(cycle) == len(parts_graph)):
                break
        else:
            err = AssemblyError("No construct found involving all parts")
            err.graph = graph
            err.mix = self
            raise err
        if len(cycle) == 0:
            raise ValueError("No solution found - a connector may be missing.")

        selected_connectors = [
            self.fragments_dict[n].original_construct
            for (node1, node2) in zip(cycle, cycle[1:] + [cycle[0]])
            for n in all_paths[node1][node2][1:-1]
        ]

        # initialize the mix with the selected connectors
        self.constructs = original_constructs + selected_connectors
        self.compute_fragments()
        self.initialize()
        return selected_connectors
