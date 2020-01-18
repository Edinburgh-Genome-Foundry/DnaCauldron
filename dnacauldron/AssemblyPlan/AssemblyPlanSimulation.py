from flametree import file_tree
import proglog
import pandas
from ..tools import (
    format_data_dicts_records_for_spreadsheet,
    format_value_for_spreadsheet,
)
from ..biotools import write_record
from ..Assembly.AssemblyReportWriter import AssemblyReportWriter


class AssemblyPlanSimulation:
    def __init__(
        self, assembly_plan, assembly_simulations, sequence_repository=None
    ):
        self.assembly_plan = assembly_plan
        self.assembly_simulations = assembly_simulations
        self.sequence_repository = sequence_repository

    def compute_all_construct_data_dicts(self):
        return [
            data_dict
            for simulation in self.assembly_simulations
            for data_dict in simulation.compute_all_construct_data_dicts()
        ]

    def compute_summary_dataframe(self):
        first_columns = [
            "assembly_name",
            "construct_id",
            "assembly_level",
            "construct_size",
            "number_of_parts",
        ]
        construct_data_dicts = self.compute_all_construct_data_dicts()
        extra_data_columns = [
            field
            for data_dict in construct_data_dicts
            for field in data_dict
            if field not in (first_columns + ["parts"])
        ]
        extra_data_columns = sorted(set(extra_data_columns))
        columns = first_columns + extra_data_columns + ["parts"]
        data = format_data_dicts_records_for_spreadsheet(construct_data_dicts)
        return pandas.DataFrame(data, columns=columns)

    def write_report(
        self,
        target,
        folder_name="auto",
        assembly_report_writer="default",
        logger="bar",
    ):
        if assembly_report_writer == "default":
            assembly_report_writer = AssemblyReportWriter(
                # We'll write all records into one folder for the whole plan
                write_part_records=False
            )
        logger = proglog.default_bar_logger(logger)
        if folder_name == "auto":
            folder_name = self.assembly_plan.name + "_assembly_report"
        report_root = file_tree(target)._dir(folder_name, replace=True)
        self._write_assembly_reports(
            report_root, assembly_report_writer, logger=logger
        )
        self._write_errors_spreadsheet(report_root)
        self._write_all_required_parts(report_root)
        self._write_summary_spreadsheet(report_root)
        self._write_assembly_plan_spreadsheets(report_root)
        if target == "@memory":
            return report_root._close()

    def _get_file_name(self, filename):
        name = self.assembly_plan.name
        prefix = (name + "_") if (name and len(name)) else ""
        return prefix + filename

    def _write_errors_spreadsheet(self, report_root):
        all_errors = [
            error
            for simulation in self.assembly_simulations
            for error in simulation.errors
        ]
        if len(all_errors) > 0:
            columns = ";".join(
                ["assembly_name", "message", "suggestion", "data", "used_in"]
            )
            all_error_rows = [
                ";".join(
                    [
                        err.assembly.name,
                        err.message,
                        err.data_as_string(),
                        " & ".join(err.assembly.dependencies["used_in"]),
                    ]
                )
                for err in all_errors
            ]
            errors_spreadsheet = report_root._file("assembly_errors.csv")
            errors_spreadsheet.write("\n".join([columns] + all_error_rows))

    def _write_assembly_reports(self, report_root, report_writer, logger):
        all_records_folder = report_root._dir("all_records")
        logger(message="Generating assemblies reports...")
        for simulation in logger.iter_bar(assembly=self.assembly_simulations):
            # TODO: skip cancelled assemblies!
            assembly_folder = report_root._dir(simulation.assembly.name)
            simulation.write_report(
                target=assembly_folder, report_writer=report_writer,
            )
            for record in simulation.construct_records:
                target = all_records_folder._file(record.id + ".gb")
                write_record(record, target.open("w"), "genbank")

    def _write_summary_spreadsheet(self, report_root):
        data = self.compute_summary_dataframe()
        file_name = self._get_file_name("summary.csv")
        data.to_csv(report_root._file(file_name).open("w"), index=False)

    def list_all_original_parts(self):
        all_parts = [
            part
            for simulation in self.assembly_simulations
            for part in simulation.list_all_parts_used()
        ]
        assemblies = [
            simulation.assembly.name
            for simulation in self.assembly_simulations
        ]
        parts_that_arent_assembled = set(all_parts).difference(set(assemblies))
        return sorted(parts_that_arent_assembled)

    def _write_all_required_parts(self, report_root):
        all_parts = self.list_all_original_parts()
        file_name = self._get_file_name("all_required_parts.txt")
        report_root._file(file_name).write("\n".join(all_parts))

    def _write_all_required_parts_records(self, report_root):
        all_parts = self.list_all_original_parts()
        part_records = self.sequence_repository.get_records(all_parts)
        records_dir = self._dir("part_records")
        for part_record in part_records:
            filename = part_record.id + ".gb"
            target = records_dir._file(filename)
            write_record(part_record, target, "genbank")

    def _write_assembly_plan_spreadsheets(self, report_root):
        data = self.compute_summary_dataframe()
        assembly_plan_has_single_level = set(data["assembly_level"]) == {1}
        for level, subdata in data.groupby("assembly_level"):
            construct_parts = [
                (row.construct_id, row.parts.split(" & "))
                for i, row in subdata.iterrows()
            ]
            if assembly_plan_has_single_level:
                file_name = self._get_file_name("assembly_plan.csv")
            else:
                file_name = "assembly_plan_level_%s.csv" % level
                file_name = self._get_file_name(file_name)
            f = report_root._file(file_name)
            lines = [",".join([c] + parts) for c, parts in construct_parts]
            f.write("\n".join(["construct, parts"] + lines))
