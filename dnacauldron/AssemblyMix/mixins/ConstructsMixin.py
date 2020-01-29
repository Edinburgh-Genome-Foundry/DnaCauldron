import networkx as nx
import numpy as np
from ...Fragment.FragmentChain import FragmentChain
from ..AssemblyMixError import AssemblyMixError


class ConstructsMixin:
    """Mixin for AssemblyMix"""

    def compute_random_circular_fragments_sets(
        self, staling_cutoff=100, fragments_sets_filters=()
    ):
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
                new_graph = nx.DiGraph(
                    [
                        (
                            permutation[node_to_index[node1]],
                            permutation[node_to_index[node2]],
                        )
                        for node1, node2 in graph.edges()
                    ]
                )
                counter = 0
                for cycle in nx.simple_cycles(new_graph):
                    cycle = [antipermutation[i] for i in cycle]
                    fragments = [
                        self.fragments_dict[graph_nodes[i]] for i in cycle
                    ]
                    cycle = FragmentChain(
                        fragments, is_cycle=True
                    ).standardized()
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
                    if all(
                        fl(cycle.fragments) for fl in fragments_sets_filters
                    ):
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
                cycle = FragmentChain(cycle, is_cycle=True).standardized()
                cycle_hash = hash(cycle)
                if cycle_hash in seen_hashes:
                    continue
                seen_hashes.add(cycle_hash)
                if all(fl(cycle.fragments) for fl in fragments_sets_filters):
                    yield cycle.fragments

        return generator()

    def compute_circular_assemblies(
        self,
        fragments_sets_filters=(),
        seqrecord_filters=(),
        annotate_parts_homologies=False,
        randomize=False,
        randomization_staling_cutoff=100,
    ):
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

        annotate_parts_homologies
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
                    staling_cutoff=randomization_staling_cutoff,
                )
            else:
                fragments = self.compute_circular_fragments_sets(
                    fragments_sets_filters=fragments_sets_filters
                )
            for fragments in fragments:
                construct = self.assemble(
                    fragments,
                    circularize=True,
                    annotate_homologies=annotate_parts_homologies,
                )
                if all(fl(construct) for fl in seqrecord_filters):
                    construct.fragments = fragments
                    yield construct

        return assemblies_generator()

    def compute_linear_assemblies(
        self,
        fragments_sets_filters=(),
        min_parts=2,
        seqrecord_filters=(),
        annotate_parts_homologies=False,
    ):
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

        annotate_parts_homologies
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
                if len(path) < min_parts:
                    continue
                fragments = [self.fragments_dict[f] for f in path]
                if not all([fl(fragments) for fl in fragments_sets_filters]):
                    continue
                chain = FragmentChain(fragments).standardized()
                chain_hash = hash(chain)
                if chain_hash in seen_hashes:
                    continue
                seen_hashes.add(chain_hash)
                fragments_assembly = self.assemble(
                    fragments, annotate_homologies=annotate_parts_homologies
                )
                if all([fl(fragments_assembly) for fl in seqrecord_filters]):
                    yield (fragments_assembly)
