from flametree import file_tree
import matplotlib.pyplot as plt
from Bio import SeqIO
from dna_features_viewer import BiopythonTranslator
import pandas
from collections import defaultdict, Counter

from ..AssemblyMix import RestrictionLigationMix, NoRestrictionSiteFilter
from .plots import (plot_cuts, plot_slots_graph, AssemblyTranslator)

def name_fragment(fragment):
    """Return the name of the fragment, or `r_NAME` if the fragment is the
    reverse of another framgnet."""
    return (fragment.original_construct.name +
            ("_r" if fragment.is_reverse else ""))

def full_assembly_report(parts, target, enzyme="BsmBI", max_assemblies=40,
                         connector_records=(),
                         include_fragments=True,
                         include_parts=True,
                         include_assembly_plots=True,
                         fragments_filters='auto',
                         assemblies_prefix='assembly',
                         show_overhangs_in_graph=True,
                         show_overhangs_in_genbank=False,
                         mix_class="restriction"):
    """Write a full assembly report in a folder or a zip.

    The report contains the final sequence(s) of the assembly in Genbank format
    as well as a .csv report on all assemblies produced and PDF figures
    to allow a quick overview or diagnostic.

    Folder ``assemblies`` contains the final assemblies, ``assembly_graph``
    contains a schematic view of how the parts assemble together, folder
    ``fragments`` contains the details of all fragments produced by the enzyme
    digestion, and folder ``provided_parts`` contains the original input
    (genbanks of all parts provided for the assembly mix).

    Parameters
    ----------

    parts
      List of Biopython records representing the parts, potentially on entry
      vectors. All the parts provided should have different attributes ``name``
      as it is used to name the files.

    target
      Either a path to a folder, or to a zip file, or ``@memory`` to return
      a string representing zip data (the latter is particularly useful for
      website backends).

    enzyme
      Name of the enzyme to be used in the assembly

    max_assemblies
      Maximal number of assemblies to consider. If there are more than this
      the additional ones won't be returned.

    fragments_filters
      Fragments filters to be used to filter out fragments before looking for
      assemblies. If left to auto, fragments containing the enzyme site will
      be filtered out.

    connector_records
      List of connector records (a connector is a part that can bridge a gap
      between two other parts), from which only the essential elements to form
      an assembly will be automatically selected and added to the other parts.

    assemblies_prefix
      Prefix for the file names of all assemblies. They will be named
      ``PRE01.gb``,``PRE02.gb``, ``PRE03.gb`` where ``PRE`` is the prefix.


    """

    if mix_class == "restriction":
        mix_class = RestrictionLigationMix
    part_names = [p.name for p in parts]
    non_unique = [e for (e, count) in Counter(part_names).items() if count > 1]
    non_unique = list(set(non_unique))
    if len(non_unique) > 0:
        raise ValueError("All parts provided should have different names. "
                         "Assembly (%s) contains several times the parts %s " %
                         (" ".join(part_names), ", ".join(non_unique)))
    if fragments_filters == 'auto':
        fragments_filters = [NoRestrictionSiteFilter(enzyme)]

    report = file_tree(target, replace=True)

    assemblies_dir = report._dir("assemblies")

    mix = mix_class(parts, enzyme, fragments_filters=fragments_filters)
    if len(connector_records):
        mix.autoselect_connectors(connector_records)

    # PROVIDED PARTS
    if include_parts:
        provided_parts_dir = report._dir("provided_parts")
        for part in parts:
            linear = part.linear if hasattr(part, 'linear') else False
            ax, gr = plot_cuts(part, enzyme, linear=linear)
            f = provided_parts_dir._file(part.name + ".pdf").open('wb')
            ax.figure.savefig(f, format='pdf', bbox_inches="tight")
            plt.close(ax.figure)
            gb_file = provided_parts_dir._file(part.name + ".gb")
            SeqIO.write(part, gb_file.open('w'), 'genbank')

    # FRAGMENTS
    if include_fragments:
        fragments_dir = report._dir("fragments")
        seenfragments = defaultdict(lambda *a: 0)
        for fragment in mix.fragments:
            gr = BiopythonTranslator().translate_record(fragment)
            ax, pos = gr.plot()
            name = name_fragment(fragment)
            seenfragments[name] += 1
            file_name = "%s_%02d.pdf" % (name, seenfragments[name])
            ax.figure.savefig(fragments_dir._file(file_name).open('wb'),
                              format='pdf', bbox_inches="tight")
            plt.close(ax.figure)

    # GRAPH
    ax = plot_slots_graph(mix, with_overhangs=show_overhangs_in_graph,
                          show_missing=True)
    f = report._file('parts_graph.pdf')
    ax.figure.savefig(f.open('wb'), format='pdf', bbox_inches='tight')
    plt.close(ax.figure)

    # ASSEMBLIES
    assemblies = mix.compute_circular_assemblies(
        annotate_homologies=show_overhangs_in_genbank)
    assemblies = sorted(
        [asm for (i, asm) in zip(range(max_assemblies), assemblies)],
        key=lambda asm: str(asm.seq)
    )
    assemblies_data = []
    i_asm = list(zip(range(max_assemblies), assemblies))
    for i, asm in i_asm:
        if len(i_asm) == 1:
            name = assemblies_prefix
        else:
            name = '%s_%03d' % (assemblies_prefix, (i+1))
        assemblies_data.append(dict(
            assembly_name=name,
            parts=" & ".join([name_fragment(f_) for f_ in asm.fragments]),
            number_of_parts=len(asm.fragments),
            assembly_size=len(asm)
        ))
        SeqIO.write(asm, assemblies_dir._file(name + '.gb').open('w'),
                    'genbank')
        if include_assembly_plots:
            gr_record = AssemblyTranslator().translate_record(asm)
            ax, gr = gr_record.plot(figure_width=16)
            ax.set_title(name)
            ax.figure.savefig(assemblies_dir._file(name + '.pdf').open('wb'),
                              format='pdf', bbox_inches='tight')
            plt.close(ax.figure)
    df = pandas.DataFrame.from_records(
        assemblies_data,
        columns=['assembly_name', 'number_of_parts', 'assembly_size', 'parts']
    )
    df.to_csv(report._file('report.csv'), index=False)
    n_constructs = len(df)
    if target == '@memory':
        return n_constructs, report._close()
    else:
        return n_constructs
