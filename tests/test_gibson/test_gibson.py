import os
import dnacauldron as dc
from geneblocks import sequences_are_circularly_equal

this_directory = os.path.dirname(os.path.realpath(__file__))
sequences_fasta = os.path.join(this_directory, "gibson_sequences.fa")


def test_gibson_assembly_plan():
    repository = dc.SequenceRepository()
    repository.import_records(files=[sequences_fasta])
    assembly_plan = dc.AssemblyPlan.from_spreadsheet(
        assembly_class=dc.GibsonAssembly,
        path=os.path.join(this_directory, "gibson_plan.csv"),
    )
    plan_simulation = assembly_plan.simulate(sequence_repository=repository)
    stats = plan_simulation.compute_stats()
    assert stats["valid_assemblies"] == 3
    assert stats["errored_assemblies"] == 2
    report_writer = dc.AssemblyReportWriter(
        include_mix_graphs=True,
        include_assembly_plots=True,
        show_overhangs_in_graph=True,
        annotate_parts_homologies=True,
    )
    plan_simulation.write_report(
        target="@memory", assembly_report_writer=report_writer
    )


def test_single_gibson():
    repository = dc.SequenceRepository()
    repository.import_records(files=[sequences_fasta])
    parts = ["Frag_%d" % i for i in [1, 2, 3, 4, 5]]
    expected_record = repository.get_record("expected_sequence")
    assembly = dc.GibsonAssembly(parts=parts, homology_checker="default")
    simulation = assembly.simulate(sequence_repository=repository)
    assert len(simulation.construct_records) == 1
    simulated_record = simulation.construct_records[0]
    assert sequences_are_circularly_equal([simulated_record, expected_record])
