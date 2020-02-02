import dnacauldron as dc
import os

repository = dc.SequenceRepository()
repository.import_records(folder="parts", use_file_names_as_ids=True)
parts_list = list(repository.collections["parts"])
assembly = dc.Type2sRestrictionAssembly(
    name="randomized_combinatorial_asm",
    parts=parts_list,
    expected_constructs="any_number",
    randomize_constructs=True,
    max_constructs=2,
)

simulation = assembly.simulate(sequence_repository=repository,)
output_path = os.path.join("output", "randomized_combinatorial")
simulation.write_report(target=output_path)
print("Done! see output/randomized_combinatorial folder for the results.")
