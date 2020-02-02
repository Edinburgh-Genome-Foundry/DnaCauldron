import dnacauldron as dc
import os

repository = dc.SequenceRepository()
repository.import_records(folder="parts", use_file_names_as_ids=True)
parts_list = list(repository.collections["parts"])
assembly = dc.Type2sRestrictionAssembly(
    name="combinatorial_asm",
    parts=parts_list,
    expected_constructs="any_number",
)
simulation = assembly.simulate(sequence_repository=repository)
report_writer = dc.AssemblyReportWriter(
    include_mix_graphs=True, include_part_plots=True
)
simulation.write_report(
    target=os.path.join("output", "combinatorial"),
    report_writer=report_writer,
)
print("Done! see output/combinatorial folder for the results.")