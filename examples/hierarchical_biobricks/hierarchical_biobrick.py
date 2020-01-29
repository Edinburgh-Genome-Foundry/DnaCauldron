import os
import dnacauldron as dc

this_directory = os.path.dirname(os.path.realpath(__file__))
parts_folder = os.path.join(this_directory, "igem_parts_with_backbone")


repository = dc.SequenceRepository()
repository.import_records(folder="biobrick_parts", use_file_names_as_ids=True)
assembly_plan = dc.AssemblyPlan.from_spreadsheet(
    assembly_class=dc.BioBrickStandardAssembly,
    path="hierarchical_biobrick.csv",
)
plan_simulation = assembly_plan.simulate(sequence_repository=repository)
plan_simulation.write_report(target="output")
print ("Done! See the output/ folder for the results.")
