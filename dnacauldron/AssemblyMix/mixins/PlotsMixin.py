import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


class PlotsMixin:
    """Mixin for AssemblyMix"""


    def plot_connections_graph(self, ax=None, figsize=(20, 20)):
        """Plot the mix's graph of connections between fragments."""

        graph = self.uniquified_connection_graph

        def fragment_label(i):
            fragment = self.fragments_dict[i]
            return "\n".join(
                [
                    fragment.text_representation_in_plots().replace("_", "-")
                ]
            )

        labels = {i: fragment_label(i) for i in graph}
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=figsize)
        if len(graph):
            nx.draw_kamada_kawai(
                graph,
                labels=labels,
                font_color="k",
                edge_color="grey",
                font_size=12,
                node_color="w",
                node_size=3000,
            )
        
        if self.name is not None:
            ax.set_title(self.name, loc="left")
        return ax
    
    def plot_fragments(self, report_root):
        """Plot schemas of the different fragments in the mix."""
        fragments_dir = report_root._dir("fragments_in_" + self.name)
        seen_fragments = {}
        for fragment in self.fragments:
            name = fragment.original_part.id
            if name not in seen_fragments:
                seen_fragments[name] = 0
            seen_fragments[name] += 1
            file_name = "%s_%02d.pdf" % (name, seen_fragments[name])
            ax = fragment.plot()
            ax.figure.savefig(
                fragments_dir._file(file_name).open("wb"),
                format="pdf",
                bbox_inches="tight",
            )
            plt.close(ax.figure)
    
    # def plot_graphs(self, report_root, assembly, with_overhangs=True):
    #     file_prefix = assembly.name + "_"
    #     ax = self.plot_connections_graph()
    #     f = report_root._file(file_prefix + "connections_graph.pdf")
    #     ax.figure.savefig(f.open("wb"), format="pdf", bbox_inches="tight")
    #     plt.close(ax.figure)
