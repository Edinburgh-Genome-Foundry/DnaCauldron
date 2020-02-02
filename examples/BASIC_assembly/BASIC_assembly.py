import dnacauldron as dc

repo = dc.SequenceRepository()
repo.import_records(folder="parts_and_oligos")
plan = dc.AssemblyPlan.from_spreadsheet(
    path="basic_assembly.csv", assembly_class="from_spreadsheet"
)
simulation = plan.simulate(repo)
simulation.write_report("output")
print ("Done! see output/ folder for the results.")