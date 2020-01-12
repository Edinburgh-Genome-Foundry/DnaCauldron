import pandas
import networkx as nx


class AssemblyPlan:
    def __init__(self, assemblies, assembly_class="auto"):
        self.assemblies = assemblies
        self.assemblies_dict = {a.name: a for a in self.assemblies}
        self.compute_assemblies_levels()
        if assembly_class == "auto":
            if len(set([a.__class__ for a in assemblies])) == 1:
                assembly_class = assemblies[0].__class__
            else:
                assembly_class = None
        
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
            assembly.level = nodes_levels[assembly.name]
        self.all_parts = [n for n in graph if nodes_levels[node] == 0]
        levels = sorted(set(a.level for a in self.assemblies))
        self.levels = {
            level: [a for a in self.assemblies if a.level == level]
            for level in levels
        }

    @staticmethod
    def from_spreadsheet(
        assembly_class,
        path=None,
        dataframe=None,
        sheet_name="all",
        header=None,
    ):
        is_csv = path.lower().endswith(".csv")
        if sheet_name == "all" and not is_csv:
            excel_file = pandas.ExcelFile(path)
            return AssemblyPlan(
                [
                    assembly
                    for _sheet_name in excel_file.sheet_names
                    for assembly in AssemblyPlan.from_spreadsheet(
                        assembly_class=assembly_class,
                        path=path,
                        sheet_name=_sheet_name,
                        header=header,
                    ).assemblies
                ]
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
        return AssemblyPlan(assemblies)

    def to_spreadsheet(self, path):
        lines = [
            ",".join([asm] + parts) for asm, parts in self.assemblies.items()
        ]
        with open(path, "w") as f:
            f.write("\n".join(["construct,parts"] + lines))