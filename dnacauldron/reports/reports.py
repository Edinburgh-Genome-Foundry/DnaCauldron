from collections import defaultdict, Counter, OrderedDict
from copy import deepcopy

from flametree import file_tree
import matplotlib.pyplot as plt
from dna_features_viewer import BiopythonTranslator
import pandas
from proglog import default_bar_logger


from ..AssemblyMix import (RestrictionLigationMix, NoRestrictionSiteFilter,
                           FragmentSetContainsPartsFilter, AssemblyError)
from .plots import (plot_cuts, plot_slots_graph, AssemblyTranslator)
from ..tools import write_record

def name_fragment(fragment, mark_reverse=False):
    """Return the name of the fragment, or optionally `NAME_r` if the fragment
    is the reverse of another fragment."""
    return (fragment.original_construct.name +
            ("_r" if (fragment.is_reverse and mark_reverse) else ""))

def full_assembly_report(parts, target, enzyme="BsmBI", max_assemblies=40,
                         connector_records=(),
                         include_fragments_plots=True,
                         include_parts_plots=True,
                         include_assembly_plots=True,
                         fragments_filters='auto',
                         assemblies_prefix='assembly',
                         show_overhangs_in_graph=True,
                         show_overhangs_in_genbank=False,
                         no_skipped_parts=False,
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

    include_parts_plots, include_assembly_plots
      These two parameters control the rendering of extra figures which are
      great for troubleshooting, but not strictly necessary, and they slow
      down the report generation considerably.


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
        try:
            mix.autoselect_connectors(connector_records)
        except AssemblyError as err:
            ax = plot_slots_graph(mix, with_overhangs=show_overhangs_in_graph,
                                  show_missing=True,
                                  highlighted_parts=part_names)
            f = report._file('parts_graph.pdf')
            ax.figure.savefig(f.open('wb'), format='pdf', bbox_inches='tight')
            plt.close(ax.figure)
            raise err

    # PROVIDED PARTS
    if include_parts_plots:
        provided_parts_dir = report._dir("provided_parts")
        for part in parts:
            linear = part.linear if hasattr(part, 'linear') else False
            ax, gr = plot_cuts(part, enzyme, linear=linear)
            f = provided_parts_dir._file(part.name + ".pdf").open('wb')
            ax.figure.savefig(f, format='pdf', bbox_inches="tight")
            plt.close(ax.figure)
            gb_file = provided_parts_dir._file(part.name + ".gb")
            write_record(part, gb_file, 'genbank')

    # FRAGMENTS
    if include_fragments_plots:
        fragments_dir = report._dir("fragments")
        seenfragments = defaultdict(lambda *a: 0)
        for fragment in mix.fragments:
            gr = BiopythonTranslator().translate_record(fragment)
            ax, _ = gr.plot()
            name = name_fragment(fragment)
            seenfragments[name] += 1
            file_name = "%s_%02d.pdf" % (name, seenfragments[name])
            ax.figure.savefig(fragments_dir._file(file_name).open('wb'),
                              format='pdf', bbox_inches="tight")
            plt.close(ax.figure)

    # GRAPH
    if len(connector_records):
        highlighted_parts = part_names
    else:
        highlighted_parts = []
    ax = plot_slots_graph(mix, with_overhangs=show_overhangs_in_graph,
                          show_missing=True,
                          highlighted_parts=highlighted_parts)
    f = report._file('parts_graph.pdf')
    ax.figure.savefig(f.open('wb'), format='pdf', bbox_inches='tight')
    plt.close(ax.figure)
    graph = mix.slots_graph(with_overhangs=False)
    slots_dict = {
        s: "|".join(list(pts))
        for s, pts in mix.compute_slots().items()
    }
    non_linear_slots = [
        (slots_dict[n], "|".join([slots_dict[b] for b in graph.neighbors(n)]))
        for n in graph.nodes()
        if graph.degree(n) != 2]
    if len(non_linear_slots):
        report._file('non_linear_nodes.csv').write("\n".join(
            ["part,neighbours"] + [
                "%s,%s" % (part, neighbours)
                for part, neighbours in non_linear_slots
            ]
        ))

    # ASSEMBLIES
    filters = (FragmentSetContainsPartsFilter(part_names),)
    assemblies = mix.compute_circular_assemblies(
        annotate_homologies=show_overhangs_in_genbank,
        fragments_sets_filters=filters if no_skipped_parts else ()
    )
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
            name = '%s_%03d' % (assemblies_prefix, (i + 1))
        assemblies_data.append(dict(
            assembly_name=name,
            parts=" & ".join([name_fragment(f_) for f_ in asm.fragments]),
            number_of_parts=len(asm.fragments),
            assembly_size=len(asm)
        ))
        write_record(asm, assemblies_dir._file(name + '.gb'), 'genbank')
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
    df.to_csv(report._file('report.csv').open('w'), index=False)
    n_constructs = len(df)
    if target == '@memory':
        return n_constructs, report._close()
    else:
        return n_constructs

def full_assembly_plan_report(assembly_plan, target, part_records=None,
                              enzyme="BsmBI", assert_single_assemblies=True,
                              logger='bar', connector_records=(),
                              fail_silently=True,
                              errors_with_traceback=False,
                              **report_kwargs):
    """Makes a full report for a plan (list of single assemblies)

    Parameters
    ----------

    assembly_plan
      A list ``[('name', [parts])...]`` or a dict ``{name: [parts]}`` where
      the parts are either records, or simply part names (in that case you
      must provide the records in ``parts_records``)

    parts_records
      A dict {part_name: part_record}.

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

    **report_kwargs
      Any other parameter of ``full_assembly_report``

    """
    logger = default_bar_logger(logger)
    if isinstance(assembly_plan, list):
        assembly_plan = OrderedDict(assembly_plan)
    if isinstance(list(assembly_plan.values())[0][0], str):
        if not hasattr(part_records, 'items'):
            part_records = {r.name: r for r in part_records}
        for part in list(part_records):
            part_records[part] = deepcopy(part_records[part])
            part_records[part].name = part_records[part].id = part
        assembly_plan = OrderedDict([
            (name, [part_records[p] for p in parts])
            for name, parts in assembly_plan.items()
        ])
    root = file_tree(target, replace=True)
    all_records_folder = root._dir("all_records")
    errored_assemblies = []
    assemblies = list(assembly_plan.items())
    for asm_name, parts in logger.iter_bar(assembly=assemblies):
        asm_folder = root._dir(asm_name)
        try:
            n = full_assembly_report(parts, target=asm_folder,
                                     assemblies_prefix=asm_name,
                                     enzyme=enzyme,
                                     connector_records=connector_records,
                                     **report_kwargs)
            if assert_single_assemblies and (n != 1):
                raise ValueError("%s assemblies found instead of 1 for %s." %
                                 (n, asm_name))
            for f in asm_folder.assemblies._all_files:
                if f._extension == 'gb':
                    f.copy(all_records_folder)
        except Exception as err:
            if fail_silently:
                err_string = str(err)
                if errors_with_traceback:
                    err_string += str(err.__traceback__)
                errored_assemblies.append((asm_name, str(err)))
            else:
                raise err

    if len(errored_assemblies):
        root._file("errored_assemblies.txt").write("\n\n".join([
            "%s: %s" % (name, error)
            for name, error in errored_assemblies
        ]))
    f = root._file('assembly_plan.csv')
    f.write("construct, parts")
    all_parts = []
    for f_ in root._all_files:
        if f_._name_no_extension == 'report':
            first_row = f_.read('r').split("\n")[1].split(",")
            if len(first_row) == 4:
                name, _, _, parts = first_row
                parts = parts.split(" & ")
                all_parts += parts
                f.write("\n" + ",".join([name] + parts))
    all_parts = sorted(set(all_parts))
    root._file('all_parts.csv').write(",\n".join(all_parts))
    return errored_assemblies, root._close()
