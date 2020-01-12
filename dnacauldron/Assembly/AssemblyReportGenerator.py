import pandas
from flametree import file_tree
from ..biotools import write_record
from ..AssemblyMix import AssemblyMixError
from .AssemblyReportPlotsMixin import AssemblyReportPlotsMixin


class AssemblyReportGenerator(AssemblyReportPlotsMixin):
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

    include_parts_plots, include_assembly_plots
      These two parameters control the rendering of extra figures which are
      great for troubleshooting, but not strictly necessary, and they slow
      down the report generation considerably. They can be True, False, or
      "on_failure" to be True only if the number of assemblies differs from
      n_expected_assemblies


    """

    def __init__(
        self,
        include_fragments_plots="on_failure",
        include_parts_plots="on_failure",
        include_fragments_connection_graph="on_failure",
        include_assembly_plots=False,
        show_overhangs_in_graph=True,
        annotate_parts_homologies=True,
        mix_class="restriction",
    ):
        self.include_fragments_plots = include_fragments_plots
        self.include_parts_plots = include_parts_plots
        self.include_fragments_connection_graph = (
            include_fragments_connection_graph
        )
        self.include_assembly_plots = include_assembly_plots
        self.show_overhangs_in_graph = show_overhangs_in_graph
        self.annotate_parts_homologies = annotate_parts_homologies

    @staticmethod
    def name_fragment(fragment, mark_reverse=False):
        """Return the name of the fragment, or optionally `NAME_r` if the fragment
        is the reverse of another fragment."""
        return fragment.original_construct.name + (
            "_r" if (fragment.is_reverse and mark_reverse) else ""
        )

    def get_plots_options(self, is_failure):
        def evaluate(v):
            return is_failure if v == "on_failure" else v

        return {
            "fragments_plots": evaluate(self.include_fragments_plots),
            "fragments_connection_graph": evaluate(
                self.include_fragments_connection_graph
            ),
            "parts_plots": evaluate(self.include_parts_plots),
        }

    def extract_simulated_construct_data(self, construct):
        parts = [self.name_fragment(f_) for f_ in construct.fragments]
        return dict(
            construct_name=construct.name,
            parts=" & ".join(parts),
            number_of_parts=len(construct.fragments),
            construct_size=len(construct),
        )

    def create_nonlinear_slots_report(self, mix, report):
        graph = mix.slots_graph(with_overhangs=False)
        slots_dict = {
            s: "|".join(list(pts)) for s, pts in mix.compute_slots().items()
        }
        non_linear_slots = [
            (
                slots_dict[n],
                "|".join([slots_dict[b] for b in graph.neighbors(n)]),
            )
            for n in graph.nodes()
            if graph.degree(n) != 2
        ]

        if len(non_linear_slots):
            report._file("non_linear_nodes.csv").write(
                "\n".join(
                    ["part,neighbours"]
                    + [
                        "%s,%s" % (part, neighbours)
                        for part, neighbours in non_linear_slots
                    ]
                )
            )

    def get_data_columns(self, assembly):
        first_columns = [
            "assembly_name",
            "construct_name",
            "assembly_level",
            "construct_size",
            "number_of_parts",
        ]
        extras = sorted(assembly.get_extra_construct_data().keys())
        return first_columns + extras + ["parts"]

    def create_constructs_data_spreadsheet(
        self, assembly, constructs_data, report
    ):
        dataframe = pandas.DataFrame.from_records(
            constructs_data, columns=self.get_data_columns(assembly),
        )
        target = report._file("assembly_summary.csv")
        dataframe.to_csv(target.open("w"), index=False)

    def generate_report(self, assembly, sequences_repository, target):
        report = file_tree(target, replace=True)

        # SIMULATE, GENERATE CONSTRUCTS RECORDS

        try:
            assembly_mix, constructs_records = assembly.simulate(
                sequences_repository=sequences_repository,
                annotate_parts_homologies=self.annotate_parts_homologies,
            )
        except AssemblyMixError as err:
            parts = assembly.parts
            self.plot_slots_graph(err.mix, report, highlighted_parts=parts)
            self.plot_connections_graph(err.mix, report)
            raise err

        # CREATE THE CONSTRUCTS INFOS

        n_assemblies = len(constructs_records)
        constructs_data = []
        for i, construct_record in enumerate(constructs_records):
            name = assembly.name
            if n_assemblies > 1:
                name += "_%03d" % (i + 1)
            construct_record.name = construct_record.id = name
            data = self.extract_simulated_construct_data(construct_record)
            constructs_data.append(data)
            data.update(assembly.get_extra_construct_data())
            data.update(
                dict(
                    assembly_name=assembly.name,
                    record=construct_record,
                    assembly_level=assembly.level,
                )
            )

        n_assemblies = len(constructs_data)
        if n_assemblies > 1:
            assemblies_dir = report._dir("assemblies")
        else:
            assemblies_dir = report  # use the root

        # WRITE GENBANKs

        for construct_data in constructs_data:
            filename = construct_data["construct_name"] + ".gb"
            target = assemblies_dir._file(filename)
            write_record(construct_data["record"], target, "genbank")

        # PLOTS

        if self.include_assembly_plots:
            for construct_data in constructs_data:
                self.plot_construct(construct_data["record"], assemblies_dir)

        is_failure = n_assemblies == 0
        expected = assembly.expected_constructs
        if expected is not None and (expected != n_assemblies):
            is_failure = True

        def evaluate(v):
            return is_failure if v == "on_failure" else v

        if evaluate(self.include_fragments_plots):
            parts_records = sequences_repository.get_records(assembly.parts)
            enzymes = assembly.enzymes if hasattr(assembly, "enzymes") else []
            self.plot_provided_parts(
                report=report, parts_records=parts_records, enzymes=enzymes,
            )
        if evaluate(self.include_fragments_connection_graph):
            self.plot_fragments(report=report, mix=assembly_mix)
        if evaluate(self.include_parts_plots):
            self.plot_connections_graph(report=report, mix=assembly_mix)
        highlighted_parts = (
            [] if assembly.connectors_collection is None else assembly.parts
        )
        self.plot_slots_graph(
            assembly_mix, report, highlighted_parts=highlighted_parts
        )

        self.create_nonlinear_slots_report(assembly_mix, report)
        self.create_constructs_data_spreadsheet(
            assembly, constructs_data, report
        )

        if target == "@memory":
            return constructs_data, report._close()
        else:
            if isinstance(target, str):
                report._close()
            return constructs_data
