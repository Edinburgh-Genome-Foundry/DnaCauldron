from copy import deepcopy

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from Bio import Restriction
from dna_features_viewer import (BiopythonTranslator, CircularGraphicRecord,
                                 GraphicRecord)

from dnacauldron.tools import annotate_record
try:
    import pygraphviz
    from networkx.drawing.nx_agraph import graphviz_layout
    GRAPHVIZ_AVAILABLE = True
except ImportError:
    GRAPHVIZ_AVAILABLE = False

class AssemblyTranslator(BiopythonTranslator):

    @staticmethod
    def is_source(feature):
        return ((feature.type == 'misc_feature') and
                feature.qualifiers.get('source', False))

    @staticmethod
    def compute_feature_color(feature):
        if feature:
            if AssemblyTranslator.is_source(feature):
                return '#ff4c4c'
            else:
                return '#f9edbb'

    @staticmethod
    def compute_feature_label(feature):
        if AssemblyTranslator.is_source(feature):
            return feature.qualifiers['source']
        elif abs(feature.location.end - feature.location.start) > 100:
            label = BiopythonTranslator.compute_feature_label(feature)
            return abreviate_string(label, 30)
        else:
            return None

def abreviate_string(string, max_length=30):
    return string[:max_length] + ('' if len(string) < max_length else '...')

def plot_cuts(record, enzyme_name, linear=True, figure_width=5, ax=None):
    enzyme = Restriction.__dict__[enzyme_name]
    record = deepcopy(record)
    cuts = enzyme.search(record.seq, linear=linear)
    for cut in cuts:
        annotate_record(record, (cut, cut+1),
                        feature_type='misc_feature',
                        enzyme=enzyme_name)
    class MyTranslator(BiopythonTranslator):

        @staticmethod
        def compute_feature_label(feature):
            if abs(feature.location.end - feature.location.start) > 100:
                label = BiopythonTranslator.compute_feature_label(feature)
                return abreviate_string(label, 10)
            else:
                return feature.qualifiers.get("enzyme", None)

        @staticmethod
        def compute_feature_color(feature):
            if (feature.qualifiers.get("enzyme", False) and
                (feature.type == 'misc_feature')):
                return '#f5eaff'
            else:
                return '#ffffff'
    translator = MyTranslator()
    grecord_class = GraphicRecord if linear else CircularGraphicRecord
    graphic_record = translator.translate_record(record,
                                                 grecord_class=grecord_class)
    return graphic_record.plot(ax=ax, figure_width=figure_width)



def name_fragment(fragment):
    return (fragment.original_construct.name +
            ("_r" if fragment.is_reverse else ""))

def plot_assembly_graph(mix, ax=None, fragments_display_lim=3,
                        fragments_filters=None, figure_size=(8, 6)):
    """Plot a map of the different assemblies"""

    def normalized_end(end):
        return min(str(end), str(end.reverse_complement()))


    g = nx.Graph()
    all_fragments = [
        f for f in (mix.fragments + mix.reverse_fragments)
        if all([fl(f) for fl in fragments_filters])
    ]
    def fragments_are_equal(f1, f2):
        '''loose equality for the purpose of this method'''
        return str(f1.seq) == str(f2.seq)
    for fragment in all_fragments:
        # print fragment.seq.left_end, fragment.seq.right_end
        if None in [fragment.seq.left_end, fragment.seq.right_end]:
            continue
        left = normalized_end(fragment.seq.left_end)
        right = normalized_end(fragment.seq.right_end)
        if not g.has_edge(left, right):
            g.add_edge(left, right, fragments=[])
        fragments = g[left][right]['fragments']
        if ((not any([fragments_are_equal(fragment, f) for f in fragments])) and
            (not any([fragments_are_equal(fragment.reverse_fragment, f)
                     for f in fragments]))):
            fragments.append(fragment)


    if ax is None:
        fig, ax = plt.subplots(1, figsize=figure_size)
    ax.axis("off")
    if GRAPHVIZ_AVAILABLE:
        layout = graphviz_layout(g, 'circo')
    else:
        layout = nx.layout.circular_layout(g)
    values = list(layout.values())
    all_x = [p[0] for p in values]
    all_y = [p[1] for p in values]
    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = min(all_y), max(all_y)
    dx = 0.1 * (xmax - xmin)
    dy = 0.1 * (ymax - ymin)
    ax.set_xlim(xmin - dx, xmax + dx)
    ax.set_ylim(ymin - dy, ymax + dy)

    for (end1, end2, data) in list(g.edges(data=True)):

        x1, y1 = p1 = np.array(layout[end1])
        x2, y2 = p2 = np.array(layout[end2])
        center = 0.5 * (p1+p2)
        ax.plot([x1, x2], [y1, y2], color='gray')
        fragments = data['fragments']
        g[end1][end2]['fragments'] = fragments
        if len(fragments) <= fragments_display_lim:
            label = "\n".join([
                name_fragment(fragment)
                for fragment in fragments[:fragments_display_lim]
            ])
        else:
            label = "%d parts" % len(fragments)
        ax.annotate(label, center, verticalalignment="center",
                    horizontalalignment="center",
                    xycoords="data",
                    bbox={'facecolor': 'white', 'linewidth': 0},
                    family='Open Sans', weight='bold')

    for (end, pos) in layout.items():
        ax.annotate(end, pos, verticalalignment="center",
                    horizontalalignment="center",
                    xycoords="data",
                    bbox={'facecolor': 'white', 'linewidth': 0},
                    family='Open Sans')
    return ax, g
