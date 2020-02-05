"""Simulate an assembly plan with many different sorts of assemblies.

Note: this is a very generic script which can be used to simulate any assembly
plan, and is equivalent to running

>>> dnacauldron multi_assembly.xlsx part_sequences output --connectors=connector_sequences

at the command line.
"""

import dnacauldron as dc
import proglog
import os

connectors_path = "connector_sequences"
sequences_path = "part_sequences"
assembly_plan_path = "multi_assembly.xlsx"
target_path = "output"
verbose = True

logger = proglog.default_bar_logger("bar" if verbose else None)

logger(message="Importing the sequences...")
sequence_repository = dc.SequenceRepository()
sequence_repository.import_records(folder=sequences_path)
for subfolder_name in os.listdir(connectors_path):
    subfolder_path = os.path.join(connectors_path, subfolder_name)
    sequence_repository.import_records(
        folder=subfolder_path, collection=subfolder_name
    )

logger(message="Simulating the assembly plan...")
assembly_plan = dc.AssemblyPlan.from_spreadsheet(
    assembly_plan_path, logger=logger
)
simulation = assembly_plan.simulate(sequence_repository=sequence_repository)

logger(message="Simulation finished. Here are the stats:")
stats = simulation.compute_stats()
for key, value in sorted(stats.items())[::-1]:
    print((key + ":").ljust(21), value)

logger(message="Now generating the report...")
simulation.write_report(target_path, logger=logger)

logger(message="All done! See %s for the report" % target_path)
