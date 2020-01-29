import dnacauldron as dc

repository = dc.SequenceRepository()
repository.import_records(folder="biobrick_parts", use_file_names_as_ids=True)
assembly_plan = dc.AssemblyPlan.from_spreadsheet(
    assembly_class=dc.BioBrickStandardAssembly,
    path="hierarchical_biobrick.csv",
)
plan_simulation = assembly_plan.simulate(sequence_repository=repository)
plan_simulation.write_report(target="output")
print("Done! See the output/ folder for the results.")
