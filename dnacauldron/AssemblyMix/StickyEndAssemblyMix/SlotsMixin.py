import itertools as itt
import networkx as nx
from ...biotools import reverse_complement

class SlotsMixin:
    """Mixin for StickyEndAssemblyMix"""

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
            slots[f_slot].add(f.original_part.id)
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