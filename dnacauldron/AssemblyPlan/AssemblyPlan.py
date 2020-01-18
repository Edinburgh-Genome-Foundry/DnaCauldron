import os
import pandas
import networkx as nx
from proglog import default_bar_logger
from .AssemblyPlanSimulation import AssemblyPlanSimulation
from ..Assembly.AssemblySimulation import AssemblySimulation
from ..Assembly.AssemblySimulationError import AssemblySimulationError


class AssemblyPlan:
    def __init__(
        self, assemblies, assembly_class="auto", name="", logger="bar"
    ):
        self.assemblies = assemblies
        self.assemblies_dict = {a.name: a for a in self.assemblies}
        self.compute_assemblies_levels()
        self.logger = default_bar_logger(logger)
        if assembly_class == "auto":
            if len(set([a.__class__ for a in assemblies])) == 1:
                assembly_class = assemblies[0].__class__
            else:
                assembly_class = None
        self.name = name
        self.assembly_class = assembly_class

    def compute_assemblies_levels(self):
        graph_edges = [
            (part, assembly.name)
            for assembly in self.assemblies
            for part in assembly.parts
        ]
        graph = nx.DiGraph(graph_edges)
        if not nx.dag.is_directed_acyclic_graph(graph):
            cycle = nx.cycles.find_cycle(graph)
            raise ValueError("Circular dependency found involving %s" % cycle)
        level_0_nodes = [n for n in graph if list(graph.predecessors(n)) == []]
        nodes_levels = {node: 0 for node in graph}

        def mark_depth(node, depth):
            nodes_levels[node] = max(nodes_levels[node], depth)
            for child in graph.successors(node):
                mark_depth(child, depth + 1)

        for node in level_0_nodes:
            mark_depth(node, 0)
        for assembly in self.assemblies:
            assembly.dependencies = dict(
                level=nodes_levels[assembly.name],
                depends_on=graph.predecessors(assembly.name),
                used_in=graph.successors(assembly.name),
            )

        self.all_parts = [
            n for n in graph if len(list(graph.predecessors(n))) == 0
        ]
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
        assembly_class,
        path=None,
        dataframe=None,
        sheet_name="all",
        header=None,
        name="auto_from_filename",
    ):
        if name == "auto_from_filename":
            if path is None:
                name = "unnamed"
            else:
                filename = os.path.basename(path)
                name, _ = os.path.splitext(filename)
        is_csv = path.lower().endswith(".csv")
        if sheet_name == "all" and not is_csv:
            excel_file = pandas.ExcelFile(path)
            return AssemblyPlan(
                name=name,
                assembly_class=assembly_class,
                assemblies=[
                    assembly
                    for _sheet_name in excel_file.sheet_names
                    for assembly in AssemblyPlan.from_spreadsheet(
                        assembly_class=assembly_class,
                        path=path,
                        sheet_name=_sheet_name,
                        header=header,
                    ).assemblies
                ],
            )

        if dataframe is None:
            if is_csv:
                with open(path, "r") as f:
                    dataframe = pandas.DataFrame(
                        [line.split(",") for line in f.read().split("\n")]
                    )
            else:
                dataframe = pandas.read_excel(
                    path, sheet_name=sheet_name, header=header
                )
        assemblies = [
            assembly_class.from_dataframe_row(row)
            for i, row in dataframe.iterrows()
            if str(row[0]).lower()
            not in ["nan", "construct name", "construct", "none", ""]
        ]
        return AssemblyPlan(
            assemblies, assembly_class=assembly_class, name=name
        )

    def to_spreadsheet(self, path):
        lines = [
            ",".join([asm] + parts) for asm, parts in self.assemblies.items()
        ]
        with open(path, "w") as f:
            f.write("\n".join(["construct,parts"] + lines))

    def simulate(self, sequence_repository, logger=None):

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
                failed_dependency = cancelled_assemblies[assembly.name]
                error = AssemblySimulationError(
                    assembly=assembly,
                    message="Cancelled: depends on failed %s"
                    % failed_dependency,
                )
                simulation_result = AssemblySimulation(
                    assembly=assembly,
                    sequence_repository=sequence_repository,
                    errors=(error,),
                )
            else:
                simulation_result = assembly.simulate(sequence_repository)
            simulation_results.append(simulation_result)
            for record in simulation_result.construct_records:
                sequence_repository.constructs[record.id] = record
            if len(simulation_result.construct_records) != 1:
                for next_assembly in assembly.dependencies["used_in"]:
                    cancelled_assemblies[next_assembly] = assembly.name

        return AssemblyPlanSimulation(
            assembly_plan=self,
            assembly_simulations=simulation_results,
            sequence_repository=sequence_repository,
        )
