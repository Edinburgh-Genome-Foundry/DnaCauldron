import pandas
from flametree import file_tree
from ...biotools import write_record
from ...AssemblyMix import AssemblyMixError
from .AssemblyReportPlotsMixin import AssemblyReportPlotsMixin


class AssemblyReportWriter(AssemblyReportPlotsMixin):
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
      "on_error" to be True only if the number of assemblies differs from
      n_expected_assemblies


    """

    def __init__(
        self,
        include_fragments_plots="on_error",
        include_parts_plots="on_error",
        include_mix_graphs="on_error",
        include_part_records=True,
        include_assembly_plots=False,
        show_overhangs_in_graph=True,
        annotate_parts_homologies=True,
    ):
        self.include_fragments_plots = include_fragments_plots
        self.include_parts_plots = include_parts_plots
        self.include_mix_graphs = include_mix_graphs
        self.include_assembly_plots = include_assembly_plots
        self.show_overhangs_in_graph = show_overhangs_in_graph
        self.annotate_parts_homologies = annotate_parts_homologies
        self.include_part_records = include_part_records

    def _write_constructs_spreadsheet(self, simulation, report_root):
        dataframe = simulation.compute_summary_dataframe()
        target = report_root._file("%s_summary.csv" % simulation.assembly.name)
        dataframe.to_csv(target.open("w"), index=False)

    def _write_records(self, assembly_simulation, report_root):
        if len(assembly_simulation.construct_records) > 1:
            assemblies_dir = report_root._dir("assemblies")
        else:
            assemblies_dir = report_root

        for construct_record in assembly_simulation.construct_records:
            filename = construct_record.id + ".gb"
            target = assemblies_dir._file(filename)
            write_record(construct_record, target, "genbank")

    def _write_part_records(self, simulation, parts_records, report_root):
        provided_parts_dir = report_root._dir("provided_parts_records")
        for part in parts_records:
            write_record(part, provided_parts_dir._file(part.id + ".gb"))

    def _write_records_plots(self, assembly_simulation, report_root):

        if len(assembly_simulation.construct_records) > 1:
            plots_dir = report_root._dir("assemblies_plots")
        else:
            plots_dir = report_root
        for construct_record in assembly_simulation.construct_records:
            self.plot_construct(construct_record, plots_dir)

    def write_report(self, assembly_simulation, target):
        report_root = file_tree(target, replace=True)
        assembly = assembly_simulation.assembly

        # The 3 next lines cover the case where connectors were added to
        # the assembly, and where no assembly was found.
        used_parts = assembly_simulation.list_all_parts_used()
        parts = sorted(set(assembly.parts + used_parts))
        repository = assembly_simulation.sequence_repository
        part_records = repository.get_records(parts)

        self._write_records(assembly_simulation, report_root)
        if self.include_part_records:
            self._write_part_records(
                assembly_simulation, part_records, report_root
            )
        if self.include_assembly_plots:
            self._write_records_plots(assembly_simulation, report_root)

        errors_detected = len(assembly_simulation.errors) != 0
        plot_options = self.get_plots_options(errors_detected)

        if plot_options["parts_plots"]:
            enzymes = assembly.enzymes if hasattr(assembly, "enzymes") else []
            self.plot_provided_parts(
                report_root=report_root,
                parts_records=part_records,
                enzymes=enzymes,
            )
        if plot_options["fragments_plots"]:
            for mix in assembly_simulation.mixes:
                mix.plot_fragments(report_root=report_root)
        if plot_options["mix_graphs_plots"]:
            for mix in assembly_simulation.mixes:
                mix.plot_graphs(
                    report_root=report_root,
                    assembly=assembly,
                    with_overhangs=self.show_overhangs_in_graph,
                )
        if len(assembly_simulation.construct_records):
            self._write_constructs_spreadsheet(
                assembly_simulation, report_root
            )

        if target == "@memory":
            return report_root._close()
