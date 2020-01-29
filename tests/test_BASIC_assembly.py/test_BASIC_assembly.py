import os
import dnacauldron as dc

this_directory = os.path.dirname(os.path.realpath(__file__))
parts_and_oligos_folder = os.path.join(this_directory, "parts_and_oligos")
assembly_plan_path = os.path.join(this_directory, "assembly_plan.csv")
flawed_plan_path = os.path.join(this_directory, "assembly_plan_flawed.csv")


def test_BASIC_assembly():
    repo = dc.SequenceRepository()
    repo.import_records(folder=parts_and_oligos_folder)
    plan = dc.AssemblyPlan.from_spreadsheet(
        path=assembly_plan_path, assembly_class="from_spreadsheet"
    )
    simulation = plan.simulate(repo)
    stats = simulation.compute_stats()
    assert stats["valid_assemblies"] == 10
    assert stats["errored_assemblies"] == 0


def test_BASIC_assembly_flawed():
    repo = dc.SequenceRepository()
    repo.import_records(folder=parts_and_oligos_folder)
    plan = dc.AssemblyPlan.from_spreadsheet(
        path=flawed_plan_path, assembly_class="from_spreadsheet"
    )
    simulation = plan.simulate(repo)
    stats = simulation.compute_stats()
    assert stats["errored_assemblies"] == 3
    assert stats["cancelled_assemblies"] == 1
    assert stats["valid_assemblies"] == 8
    simulation.write_report("@memory")
