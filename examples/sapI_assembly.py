"""Assembly with SapI.

This example demonstrates an assembly using a 3bp cutter, SapI, for a change.
"""
import os
from dnacauldron import single_assembly

folder = os.path.join("data", "SapI_assembly_parts")
single_assembly(
    parts=[os.path.join(folder, filename) for filename in os.listdir(folder)],
    enzyme="SapI",
    outfile=os.path.join("output_data", "sapI_assembly.gb"),
)
print("Result written in output_data/sapI_assembly.gb")