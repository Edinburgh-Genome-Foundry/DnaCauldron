"""
"""

import itertools as itt
from copy import deepcopy

import networkx as nx
from Bio import Restriction
from Bio.Alphabet import DNAAlphabet

from .StickyEndsSeq import (StickyEndsSeqRecord,
                            digest_seqrecord_with_sticky_ends)


class FragmentsCycle:

    def __init__(self, fragments):
        self.fragments = fragments
        self.hash = hash(self)

    def reverse_complement(self):
        return FragmentsCycle([f.reverse_fragment
                               for f in self.fragments][::-1])

    def is_equivalent_to(self, other, consider_reverse=True):
        if self.hash == other.hash:
            return True
        elif consider_reverse:
            return (self.hash == self.reverse_complement().hash)
        else:
            return False

    def standardized(self):

        reverse_proportion = (sum(len(f)
                                  for f in self.fragments
                                  if f.is_reverse) /
                              float(sum(len(f) for f in self.fragments)))
        if reverse_proportion > 0.5:
            std_fragments = self.reverse_complement().fragments
        else:
            std_fragments = self.fragments

        sequences = ["%s%s%s" % (f.seq.left_end, f.seq, f.seq.right_end)
                     for f in self.fragments]
        index = sequences.index(min(sequences))
        std_fragments = std_fragments[index:] + std_fragments[:index]
        return FragmentsCycle(std_fragments)

    def __hash__(self):
        sequences = ["%s%s%s" % (f.seq.left_end, f.seq, f.seq.right_end)
                     for f in self.fragments]
        index = sequences.index(min(sequences))
        sequences = sequences[index:] + sequences[:index]
        return hash("".join(sequences))


class AssemblyMix:
    """General class for assembly mixes.

    The subclasses (RestrictionLigationMix and GibsonAssemblyMix) implement
    their own version of how the original constructs are broken into
    fragments, when two fragments will clip together, etc.

    """

    def compute_connections_graph(self):
        """Compute a graph where nodes are fragments and edges indicate
        which fragments can clip together.

        The graph (stored in self.connection_graph) is directed, an edge
        f1->f2 indicating that fragments f1 and f2 will clip in this order.
        """

        all_fragments = self.fragments + self.reverse_fragments
        self.connections_graph = nx.DiGraph()
        for fragment1, fragment2 in itt.combinations(all_fragments, 2):
            if self.will_clip_in_this_order(fragment1, fragment2):
                self.connections_graph.add_edge(fragment1, fragment2)
            if self.will_clip_in_this_order(fragment2, fragment1):
                self.connections_graph.add_edge(fragment2, fragment1)

    def compute_circular_fragments_sets(self, fragments_set_filters=()):
        """Return all lists of fragments [f1, f2, f3...fn] that can assemble
        into a circular construct.

        Fragment f1 will clip with f2 (in this order), f2 with f3... and
        fn with f1.

        This comes to finding cycles in the mix's connections graph.

        If all fragments in
        """

        def circular_fragments_generators():
            seen_hashes = set()
            for cycle in nx.simple_cycles(self.connections_graph):
                cycle = FragmentsCycle(cycle).standardized()
                if cycle.hash in seen_hashes:
                    continue
                seen_hashes.add(cycle.hash)
                if all(fl(cycle.fragments) for fl in fragments_set_filters):
                    yield cycle.fragments
        return circular_fragments_generators()

    def compute_circular_assemblies(self, fragments_set_filters=(),
                                    seqrecord_filters=()):

        def assemblies_generator():
            for fragments in self.compute_circular_fragments_sets(
                    fragments_set_filters):
                construct = self.assemble(fragments, circularize=True)
                if all(fl(construct) for fl in seqrecord_filters):
                    yield construct
        return assemblies_generator()

    def compute_reverse_fragments(self):
        self.reverse_fragments = []
        for fragment in self.fragments:
            fragment.is_reverse = False
            new_fragment = fragment.reverse_complement()
            new_fragment.is_reverse = True
            new_fragment.reverse_fragment = fragment
            fragment.reverse_fragment = new_fragment
            new_fragment.original_construct = fragment.original_construct
            self.reverse_fragments.append(new_fragment)

    def initialize(self):
        for construct in self.constructs:
            if not hasattr(construct, "linear"):
                construct.linear = True
        self.compute_fragments()
        self.compute_reverse_fragments()
        self.compute_connections_graph()

class RestrictionLigationMix(AssemblyMix):

    def __init__(self, constructs, enzyme):

        self.constructs = deepcopy(constructs)
        self.enzyme = Restriction.__dict__[enzyme]
        self.initialize()

    def compute_fragments(self):
        self.fragments = []
        for construct in self.constructs:

            digest = digest_seqrecord_with_sticky_ends(
                construct, self.enzyme, linear=construct.linear)
            for fragment in digest:
                fragment.original_construct = construct
                self.fragments.append(fragment)

    @staticmethod
    def assemble(fragments, circularize=False):
        result = sum(fragments[1:], fragments[0])
        if circularize:
            result = result.circularized()
        result.seq.alphabet = DNAAlphabet()
        return result

    @staticmethod
    def will_clip_in_this_order(fragment1, fragment2):
        return fragment1.will_clip_in_this_order_with(fragment2)


class GibsonAssemblyMix(AssemblyMix):

    def __init__(self, constructs, min_homology=15, max_homology=200):

        self.constructs = deepcopy(constructs)
        self.min_homology = min_homology
        self.max_homology = max_homology
        self.initialize()

    def compute_fragments(self):
        self.fragments = list(self.constructs)

    @staticmethod
    def assemble(fragments, circularize=False):
        return StickyEndsSeqRecord.assemble(fragments, circularize=False)
