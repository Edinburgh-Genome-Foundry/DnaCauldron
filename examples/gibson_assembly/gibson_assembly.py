import dnacauldron as dc

repository = dc.SequenceRepository()
repository.import_records(files=["gibson_sequences.fa"])
assembly_plan = dc.AssemblyPlan.from_spreadsheet(
    assembly_class=dc.GibsonAssembly, path="gibson_assembly.csv"
)
plan_simulation = assembly_plan.simulate(sequence_repository=repository)
print("Assembly stats:", plan_simulation.compute_stats())

report_writer = dc.AssemblyReportWriter(
    include_mix_graphs=True,
    include_assembly_plots=True,
    show_overhangs_in_graph=True,
    annotate_parts_homologies=True,
)
plan_simulation.write_report(
    target="output", assembly_report_writer=report_writer
)
