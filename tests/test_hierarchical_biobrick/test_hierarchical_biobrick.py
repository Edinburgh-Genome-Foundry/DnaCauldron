import os
import dnacauldron as dc

this_directory = os.path.dirname(os.path.realpath(__file__))
parts_folder = os.path.join(this_directory, "igem_parts_with_backbone")


def test_hierarchical_biobrick():
    repository = dc.SequenceRepository()
    repository.import_records(folder=parts_folder, use_file_names_as_ids=True)
    assembly_plan = dc.AssemblyPlan.from_spreadsheet(
        assembly_class=dc.BioBrickStandardAssembly,
        path=os.path.join(this_directory, "hierarchical_biobrick.csv"),
    )
    plan_simulation = assembly_plan.simulate(sequence_repository=repository)
    stats = plan_simulation.compute_stats()
    assert stats["valid_assemblies"] == 3
    report_writer = dc.AssemblyReportWriter(include_mix_graphs=True)
    plan_simulation.write_report(
        "@memory", assembly_report_writer=report_writer
    )
