import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


class PlotsMixin:
    """Mixin for StickyEndAssemblyMix"""
    def plot_slots_graph(
        self,
        ax=None,
        with_overhangs=False,
        directed=True,
        show_missing=True,
        highlighted_parts=None,
        with_arrows=True,
    ):
        """Plot a map of the different assemblies.

        Parameters
        ----------

        ax
        A matplotlib ax on which to plot. If none is provided, one is created.

        with_overhangs
        If true, the overhangs appear in the graph
        """
        slots = self.compute_slots()
        highlighted_parts = set(
            [] if highlighted_parts is None else highlighted_parts
        )
        graph = self.slots_graph(
            with_overhangs=with_overhangs, directed=directed
        )

        # Positioning - a bit complex to deal with multi-components graphs
        pos = {}
        undirected_graph = graph.to_undirected()
        components = list(nx.components.connected_components(undirected_graph))
        if components == []:
            raise ValueError(
                "Empty connections graph. This probably means your "
                "parts were filtered out, possibly because they do "
                "not contain the right enzyme sites"
            )
        max_len = 1.0 * max(len(c) for c in components)
        for i, g in enumerate(components):
            g = undirected_graph.subgraph(g)
            pos.update(
                nx.layout.kamada_kawai_layout(
                    g, center=(0, -i), scale=len(g) / max_len
                )
            )

        parts = [n for n in graph.nodes() if n in slots]

        fig, ax = plt.subplots(1, figsize=(13, 0.5 * len(parts)))
        nx.draw(
            graph,
            pos=pos,
            node_color="w",
            node_size=100,
            ax=ax,
            edge_color="#eeeeee",
        )

        # Draw a "highlights graph" above the other graph
        def highlight_slot(slot):
            return any(
                [part in highlighted_parts for part in slots.get(slot, [])]
            )

        highlighted_nodes = [n for n in graph.nodes() if highlight_slot(n)]
        if len(highlighted_nodes):
            highlighted_subgraph = graph.subgraph(highlighted_nodes)
            nx.draw(
                highlighted_subgraph,
                pos=pos,
                node_color="#3a3aad",
                node_size=300,
                ax=ax,
                edge_color="#aaaaaa",
            )
        legend = []

        def polar(xy):
            x, y = xy - np.array([0.05, 0.05])
            return (np.arctan2(x, -y), -np.sqrt(x ** 2 + y ** 2))

        sorted_pos = sorted(pos.items(), key=lambda c: polar(c[1]))
        for i, (n, (x, y)) in enumerate(sorted_pos):
            if n in parts:
                slot_parts = list(slots[n])
                color = "w" if highlight_slot(n) else "#3a3aad"
                legend.append("\n     ".join(sorted(slot_parts)))
                fontdict = {}
                if highlight_slot(n) or len(slot_parts) > 1:
                    fontdict = {"weight": "bold"}
                ax.text(
                    x,
                    y,
                    len(legend),
                    ha="center",
                    va="center",
                    color=color,
                    fontdict=fontdict,
                )
            else:
                ax.text(
                    x,
                    y,
                    n,
                    ha="center",
                    va="center",
                    color="#333333",
                    size=9,
                    fontdict=dict(family="Inconsolata"),
                )
        text = "\n".join(
            ["Parts:"]
            + [
                "%2s - %s" % (str(i + 1), name)
                for i, name in enumerate(legend)
            ]
        )
        if show_missing:
            all_mix_parts = set([f.original_part.id for f in self.fragments])
            all_slots_parts = set(
                [p for plist in slots.values() for p in plist]
            )
            missing_parts = all_mix_parts.difference(all_slots_parts)
            if len(missing_parts):
                text += "\n!!! Missing parts: " + ", ".join(missing_parts)

        ax.text(
            1.1,
            0.5,
            text,
            va="center",
            transform=ax.transAxes,
            fontdict=dict(size=12, family="Inconsolata"),
        )
        ax.set_aspect("equal")
        if self.name is not None:
            ax.set_title(self.name, loc="left")
        return ax
    
    def plot_graphs(self, report_root, assembly, with_overhangs=True):
        file_prefix = assembly.name + "_" + self.name + "_"

        # SLOTS GRAPH
        slots_graph = self.slots_graph()
        if len(slots_graph):
            highlighted_parts = []
            if hasattr(assembly, "connectors_collection"):
                if assembly.connectors_collection is not None:
                    highlighted_parts = assembly.parts
            # if hasattr(self, 'plot')
            ax = self.plot_slots_graph(
                with_overhangs=with_overhangs,
                show_missing=True,
                highlighted_parts=highlighted_parts,
            )
            
            f = report_root._file(file_prefix + "parts_graph.pdf")
            ax.figure.savefig(f.open("wb"), format="pdf", bbox_inches="tight")
            plt.close(ax.figure)

        # CONNECTIONS GRAPH
        ax = self.plot_connections_graph()
        f = report_root._file(file_prefix + "connections_graph.pdf")
        ax.figure.savefig(f.open("wb"), format="pdf", bbox_inches="tight")
        plt.close(ax.figure)
    