import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path


def plot_leveled_graph(
    levels,
    edges,
    draw_node,
    elements_positions=None,
    ax=None,
    width_factor=2.5,
    height_factor=2,
    scale=1.0,
    edge_left_space=0.015,
    edge_right_space=0.015,
    interlevel_shift=0,
    margin=None,
    **txt_kw
):
    """General function for plotting tree graphs.
    Parameters
    ----------
    levels
      A list of lists of nodes grouped by "level", i.e distance to the in the
      graph to the level 0. levels will be displayed on a same column.
    edges
      List of nodes pairs (source node, target node).
    draw_node
      A function f(x , y , node, ax, **kw) which draws something related to the
      node at the position x,y on the given Matplotlib ax.
    ax
      The matplotlib ax to use. If none is provided, a new ax is generated.
    Returns
    -------
    elements_positions, ax
      Dictionary of elements positions, matplotlib ax.
    Examples:
    ---------
    >>> def draw_node(x,y, node, ax):
            ax.text(x,y, node)
    >>> positions, ax = plot_tree_graph(
            levels=[["A","B","C"], ["D,E"], ["F"]],
            edges=[("A","D"),("B","D"),("C","E"), ("D","F"),("E","F")],
            draw_node = draw_node
        )
    """
    levels_dict = {
        element: level
        for level, elements in enumerate(levels)
        for element in elements
    }
    if elements_positions is None:
        elements_positions = {}
        for lvl, elements in enumerate(levels):
            yy = np.linspace(0, 1, len(elements) + 2)[1:-1]
            yy += interlevel_shift * (1 - 2 * (lvl % 2))
            x = 1.0 * (1 + lvl) / (len(levels) + 1)
            for y, element in zip(yy, elements):
                elements_positions[element] = (x, y)

    if ax is None:
        width = width_factor * len(levels) * scale
        height = height_factor * max([len(lvl) for lvl in levels]) * scale
        fig, ax = plt.subplots(1, figsize=(width, height))

    for element, (x, y) in elements_positions.items():
        draw_node(x, y, element, ax, **txt_kw)

    y_spans = [
        elements_positions[elements[1]][1] - elements_positions[elements[0]][1]
        for elements in levels
        if len(elements) > 1
    ]

    delta_y = 0.5 * min(y_spans) if y_spans != [] else 0
    for el1, el2 in edges:
        x1, y1 = elements_positions[el1]
        x2, y2 = elements_positions[el2]
        x1 += edge_left_space * np.sqrt(scale)
        x2 += -edge_right_space * np.sqrt(scale)
        if ((levels_dict[el2] - levels_dict[el1]) > 1) and (y1 == y2):
            patch = mpatches.PathPatch(
                Path(
                    [
                        (x1, y1),
                        (0.5 * x2 + 0.5 * x1, y1 - delta_y),
                        (0.5 * x2 + 0.5 * x1, y2 - delta_y),
                        (x2, y2),
                    ],
                    [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4],
                ),
                facecolor="none",
                lw=1 * scale,
                zorder=-1000,
            )

        else:
            patch = mpatches.PathPatch(
                Path(
                    [
                        (x1, y1),
                        (0.9 * x2 + 0.1 * x1, y1),
                        (0.1 * x2 + 0.9 * x1, y2),
                        (x2, y2),
                    ],
                    [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4],
                ),
                facecolor="none",
                lw=1 * scale,
                zorder=-1000,
            )

        ax.add_patch(patch)

    ax.axis("off")
    if margin is not None:
        xx, yy = [np.array(e) for e in zip(*elements_positions.values())]
        xmin, xmax = xx.min(), xx.max()
        dx = margin * (xmax - xmin)
        ymin, ymax = yy.min(), yy.max()
        dy = margin * (ymax - ymin)
        d = max(dx, dy)
        ax.set_xlim(xmin - d, xmax + d)
        ax.set_ylim(ymin - d, ymax + d)

    return elements_positions, ax
