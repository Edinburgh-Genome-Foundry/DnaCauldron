import os
import pandas
import networkx as nx
from proglog import default_bar_logger
from .AssemblyPlanSimulation import AssemblyPlanSimulation
from ..Assembly import ASSEMBLY_CLASS_DICT


class AssemblyPlan:
    def __init__(self, assemblies, name="plan", logger="bar"):
        """Class to represent, analyze and simulate assembly plans
        
        Parameters
        ----------

        assemblies
          List of Assembly instances.
        
        name
          Assembly plan name as it will appear in reports.
        
        logger
          Either "bar" for a progress bar, or None, or any Proglog logger.
        """
        self.assemblies = assemblies
        self._raise_an_error_if_duplicate_assembly_names()
        self.assemblies_dict = {a.name: a for a in self.assemblies}
        self._compute_assemblies_levels()
        self.logger = default_bar_logger(logger)
        self.name = name

    def _raise_an_error_if_duplicate_assembly_names(self):
        names_indices = {}
        for i, assembly in enumerate(self.assemblies):
            if assembly.name not in names_indices:
                names_indices[assembly.name] = []
            names_indices[assembly.name].append(i)
        if any(len(indices) > 1 for indices in names_indices.values()):
            duplicates = ", ".join(
                [
                    "%s (lines %s)"
                    % (name, "-".join([str(i) for i in indices]))
                    for name, indices in sorted(names_indices.items())
                    if len(indices) > 1
                ]
            )
            raise ValueError("Multiple assemblies named " + duplicates)

    def _compute_assemblies_levels(self):
        graph_edges = [
            (part, assembly.name)
            for assembly in self.assemblies
            for part in assembly.parts
        ]
        self.graph = nx.DiGraph(graph_edges)
        if not nx.dag.is_directed_acyclic_graph(self.graph):
            cycle = nx.cycles.find_cycle(self.graph)
            raise ValueError("Circular dependency found involving %s" % cycle)
        level_0_nodes = [
            n for n in self.graph if list(self.graph.predecessors(n)) == []
        ]
        nodes_levels = {node: 0 for node in self.graph}

        def mark_depth(node, depth):
            nodes_levels[node] = max(nodes_levels[node], depth)
            for child in self.graph.successors(node):
                mark_depth(child, depth + 1)

        self.all_parts = [
            n for n in self.graph if len(list(self.graph.predecessors(n))) == 0
        ]
        for node in level_0_nodes:
            mark_depth(node, 0)
        for assembly in self.assemblies:
            assembly.dependencies = dict(
                level=nodes_levels[assembly.name],
                depends_on=[
                    part
                    for part in self.graph.predecessors(assembly.name)
                    if part not in self.all_parts
                ],
                used_in=list(self.graph.successors(assembly.name)),
            )
        levels = sorted(set(a.dependencies["level"] for a in self.assemblies))
        self.levels = {
            level: [
                asm
                for asm in self.assemblies
                if asm.dependencies["level"] == level
            ]
            for level in levels
        }

    @staticmethod
    def from_spreadsheet(
        path=None,
        dataframe=None,
        assembly_class="from_spreadsheet",
        sheet_name="all",
        header=None,
        name="auto_from_filename",
        logger="bar",
        assembly_class_dict="default",
        is_csv="auto_from_filename",
        **assembly_params
    ):
        """Import an assembly plan from a spreadsheet.

        You can either read these docs or browse the examples in the repo.
        Note that this function autoselects the enzyme, based on the sites in
        each part. To explicitly set enzymes, set ``assembly.enzyme`` for each
        assembly in ``AssemblyPlan.assemblies``.


        Parameters
        ----------
        
        path
          Path to a spreadsheet file (a dataframe can be used instead).
        
        dataframe
          A pandas dataframe, possibly obtained from a spreadsheet.
        
        sheet_name
          Name of the spreadsheet's sheet on which the assembly plan is
          defined. Use "all" to load assemblies from all the sheets.

        header
          True or False, indicates whether there is a header in the
          spreadsheet.

        name
          Name of the assembly plan (leave to "auto_from_filename" to use the
          file name as assembly plan name).

        logger
          Logger of the created assembly plan. Either "bar" for a progress bar
          or None for none, or any Proglog logger. 

        assembly_params
          Extra keyword parameters which will be fed to each assembly.
        """
        if name == "auto_from_filename":
            if path is None:
                name = "unnamed"
            else:
                filename = os.path.basename(path)
                name, _ = os.path.splitext(filename)
        if is_csv == "auto_from_filename":
            is_csv = path.lower().endswith(".csv")
        if sheet_name == "all" and not is_csv:
            excel_file = pandas.ExcelFile(path)
            return AssemblyPlan(
                name=name,
                logger=logger,
                assemblies=[
                    assembly
                    for _sheet_name in excel_file.sheet_names
                    for assembly in AssemblyPlan.from_spreadsheet(
                        path=path,
                        sheet_name=_sheet_name,
                        header=header,
                        assembly_class=assembly_class,
                        is_csv=False,
                        name=name,
                        **assembly_params
                    ).assemblies
                ],
            )

        if dataframe is None:
            if is_csv:
                dataframe = pandas.read_csv(path, header=header)
                # with open(path, "r") as f:
                #     dataframe = pandas.DataFrame(
                #         [line.split(",") for line in f.read().split("\n")]
                #     )
            else:
                dataframe = pandas.read_excel(
                    path, sheet_name=sheet_name, header=header
                )
        ignore_list = [
            "nan",
            "construct name",
            "construct",
            "none",
            "assembly",
            "",
        ]
        if assembly_class == "from_spreadsheet":
            if assembly_class_dict == "default":
                assembly_class_dict = ASSEMBLY_CLASS_DICT

            def extract_assembly_class_from_row(row):
                row = list(row)
                for i, value in enumerate(row):
                    if value.startswith("class:"):
                        break
                else:
                    row = ", ".join(row)
                    msg = "Could not find assembly class in row: %s" % (row)
                    raise ValueError(msg)
                row.pop(i)
                return value[6:].strip(), row

            assembly_classes_and_rows = [
                extract_assembly_class_from_row(row)
                for i, row in dataframe.iterrows()
                if str(row[0]).lower() not in ignore_list
            ]
            assemblies = [
                assembly_class_dict[_class].from_dataframe_row(row)
                for _class, row in assembly_classes_and_rows
            ]
        else:
            assemblies = [
                assembly_class.from_dataframe_row(row)
                for i, row in dataframe.iterrows()
                if str(row[0]).lower() not in ignore_list
            ]
        return AssemblyPlan(assemblies, name=name, logger=logger)

    def to_spreadsheet(self, path):
        lines = [",".join([asm.name] + asm.parts) for asm in self.assemblies]
        with open(path, "w") as f:
            f.write("\n".join(["construct,parts"] + lines))

    def to_dataframe(self):
        """Return a dataframe describing the assembly plan."""
        sorted_assemblies = sorted(
            self.assemblies, key=lambda a: (a.dependencies["level"], a.name)
        )
        return pandas.DataFrame.from_records(
            [
                {"assembly": assembly.name, "parts": ", ".join(assembly.parts)}
                for assembly in sorted_assemblies
            ],
            columns=["assembly", "parts"],
            index="assembly",
        )

    def simulate(self, sequence_repository):
        """Simulate the whole assembly plan, return an AssemblyPlanSimulation.
        """

        ordered_assemblies = [
            assembly
            for level in sorted(self.levels)
            for assembly in self.levels[level]
        ]
        self.logger(message="Simulating assembly plan %s..." % self.name)
        simulation_results = []
        cancelled_assemblies = {}  # cancelled because dependencies failed
        for assembly in self.logger.iter_bar(assembly=ordered_assemblies):
            if assembly.name in cancelled_assemblies:
                continue
            simulation_result = assembly.simulate(sequence_repository)
            simulation_results.append(simulation_result)
            for record in simulation_result.construct_records:
                sequence_repository.add_record(record, collection="constructs")
            if len(simulation_result.errors):
                for next_assembly in assembly.dependencies["used_in"]:
                    cancelled_assemblies[next_assembly] = assembly.name

        return AssemblyPlanSimulation(
            assembly_plan=self,
            assembly_simulations=simulation_results,
            sequence_repository=sequence_repository,
            cancelled=[
                AssemblySimulationCancellation(assembly_name, dependency)
                for (assembly_name, dependency) in cancelled_assemblies.items()
            ],
        )


class AssemblySimulationCancellation:
    def __init__(self, assembly_name, failed_dependency):
        self.assembly_name = assembly_name
        self.failed_dependency = failed_dependency
