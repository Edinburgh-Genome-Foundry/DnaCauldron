import dnacauldron as dc

repository = dc.SequenceRepository()
repository.import_records(folder="parts", use_file_names_as_ids=True)
assembly_plan = dc.AssemblyPlan.from_spreadsheet(
    assembly_class=dc.Type2sRestrictionAssembly,
    path="golden_gate_two_levels.csv",
)
plan_simulation = assembly_plan.simulate(sequence_repository=repository)
plan_simulation.write_report("output")
