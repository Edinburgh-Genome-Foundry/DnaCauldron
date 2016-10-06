import networkx as nx
import matplotlib.pyplot as plt

def plot_mix_connections_graph(mix=None, graph=None, min_fragment_size=0,
                               other_nodes_conditions=(),
                               labels_fun=None,
                               direct_sense_color="#ffaaaa",
                               reverse_sense_color="#aaaaff",
                               ax=None):
    if mix is not None:
        graph = mix.connections_graph
    nodes_conditions = list(other_nodes_conditions)
    nodes_conditions.append(lambda node: len(node) > min_fragment_size)

    if ax is None:
        fig, ax = plt.subplots(1, figsize=(9, 9))
    subgraph = graph.subgraph([
        node
        for node in graph.nodes()
        if all(condition(node) for condition in nodes_conditions)
    ])
    nodes_color = [
        direct_sense_color if node.is_reverse else reverse_sense_color
        for node in subgraph.nodes()
    ]
    positions = nx.layout.spring_layout(subgraph, iterations=200)

    if labels_fun is None:
        nx.draw(subgraph, positions,
                node_color=nodes_color, node_size=70, ax=ax)
    else:

        labels = {
            node: labels_fun(node)
            for node in subgraph.nodes()
        }
        ax.axis("off")
        plot = nx.draw_networkx_nodes(subgraph, positions,
                                      node_color=nodes_color,
                                      node_size=200, ax=ax)
        plot.set_edgecolor("w")
        nx.draw_networkx_edges(subgraph, positions)
        nx.draw_networkx_labels(
            mix.connections_graph, positions, labels, font_size=9)
