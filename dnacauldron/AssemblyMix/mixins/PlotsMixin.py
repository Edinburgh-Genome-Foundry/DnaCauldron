import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


class PlotsMixin:

    def plot_connections_graph(self, ax=None, figsize=(20, 20)):

        graph = self.uniquified_connection_graph

        def fragment_label(i):
            fragment = self.fragments_dict[i]
            return "\n".join(
                [
                    str(fragment.seq.left_end),
                    r"$\bf{%s}$" % fragment.original_part.id,
                    str(fragment.seq.right_end),
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
