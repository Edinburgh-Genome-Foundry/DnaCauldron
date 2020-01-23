import os

this_directory = os.path.dirname(os.path.realpath(__file__))
parts_folder = os.path.join(this_directory, "parts")
import dnacauldron as dc


def test_type2s_hierarchical():
    repository = dc.SequenceRepository()
    repository.import_records(folder=parts_folder, use_file_names_as_ids=True)
    assembly_plan = dc.AssemblyPlan.from_spreadsheet(
        assembly_class=dc.Type2sRestrictionAssembly,
        path=os.path.join(this_directory, "type2s_two-level.csv"),
    )
    assert sorted(assembly_plan.levels) == [1, 2]
    plan_simulation = assembly_plan.simulate(sequence_repository=repository)
    stats = plan_simulation.compute_stats()
    assert stats["valid_assemblies"] == 4


def test_type2s_hierarchical_flawed():
    repository = dc.SequenceRepository()
    repository.import_records(folder=parts_folder, use_file_names_as_ids=True)
    assembly_plan = dc.AssemblyPlan.from_spreadsheet(
        assembly_class=dc.Type2sRestrictionAssembly,
        path=os.path.join(this_directory, "type2s_two-level_flawed.xls"),
    )
    assert sorted(assembly_plan.levels) == [1, 2]
    plan_simulation = assembly_plan.simulate(sequence_repository=repository)
    stats = plan_simulation.compute_stats()
    assert stats["valid_assemblies"] == 2 # Construct 1 and its dependant
    assert stats["errored_assemblies"] == 1
    assert stats["cancelled_assemblies"] == 1

    report_writer = dc.AssemblyReportWriter(
        include_fragments_plots=True,
        include_parts_plots=True,
        include_mix_graphs=True,
        include_assembly_plots=True,
        show_overhangs_in_graph=True,
        annotate_parts_homologies=True,
    )
    plan_simulation.write_report(
        target="@memory", assembly_report_writer=report_writer
    )