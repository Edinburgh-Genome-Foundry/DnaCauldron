from .AssemblyReportWriter import AssemblyReportWriter
from ..tools import format_data_dicts_records_for_spreadsheet
import pandas


class AssemblySimulation:
    """Class to represent and report on the simulation of a single assembly.

    Instances are the result of ``assembly.simulate()``.


    Parameters
    ----------

    assembly
      The Assembly instance from which this is the simulation.

    sequence_repository
      The SequenceRepository used to get records for the simulation.

    construct_records
      List of Biopython records (or, sometimes, StickyEndFragment records)
      of the final constructs predicted by the simulation.

    mixes
      A list of AssemblyMix instances generated during the simulation (they
      can be plotted at report writing time).

    warnings
      List of AssemblyFlaw instances that will be flagged as warnings in
      reports and summaries.

    errors
      List of AssemblyFlaw instances that will be flagged as errors in
      reports and summaries.
    """

    def __init__(
        self,
        assembly,
        sequence_repository,
        construct_records=(),
        mixes=(),
        warnings=(),
        errors=(),
    ):
        self.assembly = assembly
        self.construct_records = construct_records
        self.sequence_repository = sequence_repository
        self.mixes = list(mixes)
        self.errors = list(errors)
        self.warnings = list(warnings)

    @staticmethod
    def fragment_part(fragment, mark_reverse=False):
        """Return the name of the fragment, or optionally `NAME_r` if the
        fragment is the reverse of another fragment."""
        if hasattr(fragment, "original_part"):
            name = fragment.original_part.id
        else:
            name = fragment.id
        return name + ("_r" if (fragment.is_reversed and mark_reverse) else "")

    def list_all_parts_used(self):
        """List all parts involved in at least one of the predicted constructs."""
        parts = [
            self.fragment_part(fragment)
            for construct_record in self.construct_records
            for fragment in construct_record.fragments
        ]
        return sorted(set(parts))

    def compute_construct_data_dict(self, construct_record):
        """Return a dictionary with infos on a single construct.

        fields: construct_id, parts, number_of_parts, construct_size,
        assembly_name, depends_on, used_in, assembly_level.
        """
        return dict(
            construct_id=construct_record.id,
            parts=[self.fragment_part(f) for f in construct_record.fragments],
            number_of_parts=len(construct_record.fragments),
            construct_size=len(construct_record),
            assembly_name=self.assembly.name,
            depends_on=self.assembly.dependencies["depends_on"],
            used_in=self.assembly.dependencies["used_in"],
            assembly_level=self.assembly.dependencies["level"],
            **self.assembly.get_extra_construct_data()
        )

    def compute_all_construct_data_dicts(self):
        """Return a list of dictionnaries with infos on a each construct.

        Fields: construct_id, parts, number_of_parts, construct_size,
        assembly_name, depends_on, used_in, assembly_level.
        """
        return [
            self.compute_construct_data_dict(construct_record)
            for construct_record in self.construct_records
        ]

    def compute_summary_dataframe(self):
        """Return a Pandas dataframe with infos on each construct."""
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
            if field not in first_columns + ["parts"]
        ]
        extra_data_columns = sorted(set(extra_data_columns))
        columns = first_columns + extra_data_columns + ["parts"]
        data = format_data_dicts_records_for_spreadsheet(construct_data_dicts)
        return pandas.DataFrame(data, columns=columns)

    def write_report(self, target, report_writer="default"):
        """Write a comprehensive simulation report in a folder or a zip file.

        Parameters
        ----------

        target
          Either a path to a folder, to a zip file, or ``"@memory"`` to write
          into a virtual zip file whose raw data is then returned.

        report_writer
          Either the "default" or any AssemblyReportWriter instance.

        Returns
        -------

        zip_data
          binary zip data (if target="@memory") else None.
        """
        if report_writer == "default":
            report_writer = AssemblyReportWriter()
        return report_writer.write_report(assembly_simulation=self, target=target)
