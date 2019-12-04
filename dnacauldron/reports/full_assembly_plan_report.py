from collections import OrderedDict
from copy import deepcopy

from flametree import file_tree
from proglog import default_bar_logger

from ..tools import autoselect_enzyme
from .full_assembly_report import full_assembly_report


def full_assembly_plan_report(
    assembly_plan,
    target,
    part_records=None,
    enzyme="autoselect",
    assert_single_assemblies=True,
    logger="bar",
    connector_records=(),
    fail_silently=True,
    errors_with_traceback=False,
    **report_kwargs
):
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
      Any other parameter of ``full_assembly_report``. For instance:
      include_fragments_plots, include_parts_plot, include_assembly_plots

    Returns
    -------

    errored_assemblies,zip_data
      list of errored assemblies with errors, and binary zip data (or None if
      the target is not "@memory")
      
    """
    logger = default_bar_logger(logger)
    if isinstance(assembly_plan, list):
        assembly_plan = OrderedDict(assembly_plan)
    if isinstance(list(assembly_plan.values())[0][0], str):
        if not hasattr(part_records, "items"):
            part_records = {r.name: r for r in part_records}
        for part in list(part_records):
            part_records[part] = deepcopy(part_records[part])
            part_records[part].name = part_records[part].id = part
        assembly_plan = OrderedDict(
            [
                (name, [part_records[p] for p in parts])
                for name, parts in assembly_plan.items()
            ]
        )
    root = file_tree(target, replace=True)
    all_records_folder = root._dir("all_records")
    errored_assemblies = []
    assemblies = list(assembly_plan.items())
    selected_enzymes = []  # Used to keep track of autoselected enzymes
    for asm_name, parts in logger.iter_bar(assembly=assemblies):
        if enzyme == "autoselect":
            selected_enzyme = autoselect_enzyme(parts)
            selected_enzymes.append((asm_name, selected_enzyme))
        else:
            selected_enzyme = enzyme
        asm_folder = root._dir(asm_name)
        try:
            n = full_assembly_report(
                parts,
                target=asm_folder,
                assemblies_prefix=asm_name,
                enzyme=selected_enzyme,
                connector_records=connector_records,
                n_expected_assemblies=1 if assert_single_assemblies else None,
                **report_kwargs
            )
            if assert_single_assemblies and (n != 1):
                raise ValueError(
                    "%s assemblies found instead of 1 for %s." % (n, asm_name)
                )
            for f in asm_folder.assemblies._all_files:
                if f._extension == "gb":
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
        root._file("errored_assemblies.txt").write(
            "\n\n".join(
                [
                    "%s: %s" % (name, error)
                    for name, error in errored_assemblies
                ]
            )
        )
    f = root._file("assembly_plan.csv")
    f.write("construct, parts")
    all_parts = []
    for f_ in root._all_files:
        if f_._name_no_extension == "report":
            first_row = f_.read("r").split("\n")[1].split(",")
            if len(first_row) == 4:
                name, _, _, parts = first_row
                parts = parts.split(" & ")
                all_parts += parts
                f.write("\n" + ",".join([name] + parts))
    all_parts = sorted(set(all_parts))
    root._file("all_parts.csv").write(",\n".join(all_parts))
    if enzyme == "autoselect":
        root._file("selected_enzymes_per_construct.csv").write(
            ",\n".join([",".join(selection) for selection in selected_enzymes])
        )
    return errored_assemblies, root._close()
