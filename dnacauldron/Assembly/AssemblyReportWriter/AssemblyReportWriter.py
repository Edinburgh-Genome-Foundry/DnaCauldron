import pandas
from flametree import file_tree
from ...biotools import write_record
from ...AssemblyMix import AssemblyMixError
from .AssemblyReportPlotsMixin import AssemblyReportPlotsMixin


class AssemblyReportWriter(AssemblyReportPlotsMixin):
    """Class to configure assembly simulation report writing.

    Responsible for writing the final sequence(s) of the assembly in Genbank
    format as well as a .csv report on all assemblies produced and PDF figures
    to allow a quick overview or diagnostic.

    Folder ``assemblies`` contains the final assemblies, ``assembly_graph``
    contains a schematic view of how the parts assemble together, folder
    ``fragments`` contains the details of all fragments produced by the enzyme
    digestion, and folder ``provided_parts`` contains the original input
    (genbanks of all parts provided for the assembly mix).

    Parameters
    ----------

    include_fragment_plots
      Either True/False/"on_error" to plot schemas of the fragments used in
      the different AssemblyMix throughout the simulation.

    include_part_plots
      Either True/False/"on_error" to plot schemas of the parts used, possibly
      with restriction sites relevant to the AssemblyMix.

    include_mix_graphs
      Either True/False/"on_error" to plot representations of fragment
      connectivity in the AssemblyMix created during the simulation.

    include_part_records
      True/False to include the parts records in the simulation results (makes
      for larger folders and zips, but is better for traceability).

    include_assembly_plots
      True/False to include assembly schemas in the reports (makes the
      report generation slower, but makes it easier to check assemblies at a
      glance).

    show_overhangs_in_graph
      If true, the AssemblyMix graph representations will display the sequence
      of all fragment overhangs.

    include_errors_spreadsheet
      If true and there are errors, an errors spreadsheet will be added to the
      report.

    include_warnings_spreadsheet
      If true and there are warnings, a warnings spreadsheet will be added to
      the report.

    include_pdf_report
      If true, a PDF report file is also generated.
    """

    def __init__(
        self,
        include_fragment_plots="on_error",
        include_part_plots="on_error",
        include_mix_graphs="on_error",
        include_part_records=True,
        include_assembly_plots=False,
        show_overhangs_in_graph=True,
        annotate_parts_homologies=True,
        include_errors_spreadsheet=True,
        include_warnings_spreadsheet=True,
        include_pdf_report=False,
    ):
        self.include_fragment_plots = include_fragment_plots
        self.include_part_plots = include_part_plots
        self.include_mix_graphs = include_mix_graphs
        self.include_assembly_plots = include_assembly_plots
        self.show_overhangs_in_graph = show_overhangs_in_graph
        self.include_part_records = include_part_records
        self.annotate_parts_homologies = annotate_parts_homologies
        self.include_errors_spreadsheet = include_errors_spreadsheet
        self.include_warnings_spreadsheet = include_warnings_spreadsheet
        self.include_pdf_report = include_pdf_report

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
            if hasattr(construct_record, "as_biopython_record"):
                construct_record = construct_record.as_biopython_record()
            filename = construct_record.id + ".gb"
            target = assemblies_dir._file(filename)
            write_record(construct_record, target, "genbank")

    def _write_part_records(self, simulation, parts_records, report_root):
        provided_parts_dir = report_root._dir("provided_parts_records")
        for part in parts_records:
            if hasattr(part, "as_biopython_record"):
                part = part.as_biopython_record()
            write_record(part, provided_parts_dir._file(part.id + ".gb"))

    def _write_records_plots(self, assembly_simulation, report_root):

        if len(assembly_simulation.construct_records) > 1:
            plots_dir = report_root._dir("assemblies_plots")
        else:
            plots_dir = report_root
        for construct_record in assembly_simulation.construct_records:
            if hasattr(construct_record, "as_biopython_record"):
                construct_record = construct_record.as_biopython_record()
            self.plot_construct(construct_record, plots_dir)

    def _write_errors_spreadsheet(self, simulation, report_root, error_type="error"):
        errors = simulation.errors if error_type == "error" else simulation.warnings
        if len(errors) > 0:
            columns = ";".join(["assembly_name", "message", "suggestion", "data"])
            all_error_rows = [
                ";".join([err.assembly.name, err.message, err.data_as_string(),])
                for err in errors
            ]
            filename = "%s.csv" % error_type
            errors_spreadsheet = report_root._file(filename)
            errors_spreadsheet.write("\n".join([columns] + all_error_rows))

    def write_report(self, assembly_simulation, target):
        """Write a comprehensive report for an AssemblySimulation instance.

        ``target`` can be either a path to a folder, to a zip file, or
        ``"@memory"`` to write into a virtual zip file whose raw data is then
        returned.
        """
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
            self._write_part_records(assembly_simulation, part_records, report_root)
        if self.include_assembly_plots:
            self._write_records_plots(assembly_simulation, report_root)

        errors_detected = len(assembly_simulation.errors) != 0
        plot_options = self._get_plots_options(errors_detected)

        if plot_options["parts_plots"]:
            enzymes = assembly.enzymes if hasattr(assembly, "enzymes") else []
            self.plot_provided_parts(
                report_root=report_root, parts_records=part_records, enzymes=enzymes,
            )
        if plot_options["fragment_plots"]:
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
            self._write_constructs_spreadsheet(assembly_simulation, report_root)
        if self.include_errors_spreadsheet:
            self._write_errors_spreadsheet(
                assembly_simulation, report_root, error_type="error"
            )
        if self.include_warnings_spreadsheet:
            self._write_errors_spreadsheet(
                assembly_simulation, report_root, error_type="warnings"
            )

        if (target == "@memory") or str(target).endswith(".zip"):
            return report_root._close()
