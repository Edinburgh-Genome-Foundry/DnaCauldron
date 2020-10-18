import os
import dnacauldron as dc

this_directory = os.path.join("tests", "test_hierarchical_type2s")
parts_folder = os.path.join(this_directory, "parts")


def test_single_assembly(tmpdir):
    repository = dc.SequenceRepository()
    repository.import_records(folder=parts_folder, use_file_names_as_ids=True)
    assembly_plan = dc.AssemblyPlan.from_spreadsheet(
        assembly_class=dc.Type2sRestrictionAssembly,
        path=os.path.join(this_directory, "type2s_two-level.csv"),
    )
    plan_simulation = assembly_plan.simulate(sequence_repository=repository)
    stats = plan_simulation.compute_stats()
    report_writer = dc.AssemblyReportWriter(
        include_fragment_plots=False,
        include_part_plots=False,
        include_mix_graphs=False,
        include_assembly_plots=False,
        show_overhangs_in_graph=False,
        annotate_parts_homologies=False,
        include_pdf_report=True,
    )
    plan_simulation.write_report(target="@memory", assembly_report_writer=report_writer)
