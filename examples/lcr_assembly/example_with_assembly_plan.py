import os
import dnacauldron as dc

repo = dc.SequenceRepository()
files = ["RFP_GFP_plasmid_parts.fa", "RFP_GFP_plasmid_BOs.fa"]
repo.import_records(files=files)
plan = dc.AssemblyPlan.from_spreadsheet(path="assembly_plan.csv")
simulation = plan.simulate(repo)
stats = simulation.compute_stats()

simulation.write_report("output/")
print ("Done! see output/ folder for the results.")