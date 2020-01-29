import dnacauldron as dc

repository = dc.SequenceRepository()
repository.import_records(folder="parts", use_file_names_as_ids=True)
parts_list = list(repository.collections["parts"])
assembly = dc.Type2sRestrictionAssembly(
    name="combinatorial_asm",
    parts=parts_list,
    expected_constructs="any_number",
)
simulation = assembly.simulate(sequence_repository=repository)
simulation.write_report(target="output")
