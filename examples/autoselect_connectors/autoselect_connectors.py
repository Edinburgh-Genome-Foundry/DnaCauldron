import dnacauldron as dc

repository = dc.SequenceRepository()
repository.import_records(
    collection="parts",
    folder="genetic_parts",
    use_file_names_as_ids=True,
    topology="circular",
)
repository.import_records(
    collection="emma_connectors",
    folder="emma_connectors",
    use_file_names_as_ids=True,
)
all_parts = list(repository.collections["parts"])  # part names, no connectors
assembly = dc.Type2sRestrictionAssembly(
    name="assembly_with_connectors",
    parts=all_parts,
    connectors_collection="emma_connectors",
)
simulation = assembly.simulate(sequence_repository=repository)

simulation.write_report("output")
