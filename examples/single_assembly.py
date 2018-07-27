from dnacauldron import single_assembly
import os

parts_path = [
    os.path.join("data", "assemblies", part)
    for part in ["partA.gb", "partB.gb", "partC.gb", "receptor.gb"]
]
single_assembly(
    parts=parts_path,
    enzyme="BsmBI",
    outfile=os.path.join("output_data", "final_sequence.gb")
)
