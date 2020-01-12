from proglog import default_bar_logger
from flametree import file_tree
import pandas
from ..biotools import write_record
from ..Assembly import AssemblyReportGenerator

class AssemblyPlanReportGenerator:
    def __init__(
        self,
        assembly_report_generator="default",
        assert_single_assemblies=True,
        logger="bar",
        fail_silently=True,
        errors_with_traceback=False,
    ):
        """
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
        """
        if assembly_report_generator == "default":
            assembly_report_generator = AssemblyReportGenerator()
        self.assembly_report_generator = assembly_report_generator
        self.logger = default_bar_logger(logger)
        self.assert_single_assemblies = assert_single_assemblies
        self.fail_silently = fail_silently
        self.errors_with_traceback = errors_with_traceback

    def generate_report(self, assembly_plan, sequence_repository, target):
        logger = self.logger
        root = file_tree(target, replace=True)
        all_records_folder = root._dir("all_records")
        errored_assemblies = []
        logger(message="Processing assemblies...")
        assemblies_data = []
        ordered_assemblies = [
            assembly
            for level in sorted(assembly_plan.levels)
            for assembly in assembly_plan.levels[level]
        ]
        for assembly in logger.iter_bar(assembly=ordered_assemblies):
            assembly_folder = root._dir(assembly.name)
            try:
                report_generator = self.assembly_report_generator
                constructs_data = report_generator.generate_report(
                    target=assembly_folder,
                    assembly=assembly,
                    sequence_repository=sequence_repository,
                )
                assemblies_data.extend(constructs_data)
                if assembly.expected_constructs is not None:
                    n_constructs = len(constructs_data)
                    if assembly.expected_constructs != n_constructs:
                        raise ValueError(
                            "%s assemblies found instead of 1 for %s."
                            % (n_constructs, assembly.name)
                        )
                for construct_data in constructs_data:
                    record = construct_data["record"]
                    sequence_repository.constructs[record.id] = record
                    target = all_records_folder._file(record.id + ".gb")
                    write_record(record, target.open("w"), "genbank")
            except Exception as err:
                if self.fail_silently:
                    err_string = str(err)
                    if self.errors_with_traceback:
                        err_string += str(err.__traceback__)
                    errored_assemblies.append((assembly.name, str(err)))
                else:
                    raise err
        logger(message="Wrapping up...")

        if len(errored_assemblies):
            errors = ["%s: %s" % (name, er) for name, er in errored_assemblies]
            root._file("errored_assemblies.txt").write("\n\n".join(errors))

        columns = self.assembly_report_generator.get_data_columns(assembly)
        data = pandas.DataFrame(assemblies_data, columns=columns)
        assembly_plan_has_single_level = set(data["assembly_level"]) == {1}
        # ASSEMBLY PLAN SPREADSHEETS

        for level, subdata in data.groupby("assembly_level"):
            construct_parts = [
                (row.construct_id, row.parts.split(" & "))
                for i, row in subdata.iterrows()
            ]
            if assembly_plan_has_single_level:
                file_name = "assembly_plan.csv"
            else:
                file_name = "assembly_plan_level_%d.csv" % level
            f = root._file(file_name)
            lines = [",".join([c] + parts) for c, parts in construct_parts]
            f.write("\n".join(["construct, parts"] + lines))
        
        # GATHER ALL LEVEL 1 PARTS

        level_1_data = data[data["assembly_level"] == 1]
        all_parts = [
            part
            for i, row in level_1_data.iterrows()
            for part in row.parts.split(" & ")
        ]
        all_parts = sorted(set(all_parts))
        root._file("all_parts.csv").write(",\n".join(all_parts))

        # PLOT SUMMARY

        if (data.assembly_name == data.construct_id).all():
            # Means only 1 construct per assembly. Remove redundant column.
            data = data[columns[1:]]
        
        if assembly_plan_has_single_level:
            columns = [col for col in data.columns if col != "assembly_level"]
            data = data[columns]
        data.to_csv(root._file("summary.csv").open("w"), index=False)

        return assemblies_data, errored_assemblies, root._close()
