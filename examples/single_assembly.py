from dnacauldron import single_assembly
import os
os.chdir("data")
single_assembly(
    parts_filenames=["partA.gb", "partB.gb", "partC.gb"],
    receptor_filename="receptor.gb",
    outfile=os.path.join("..", "final_sequence.gb"),
    enzyme="BsmBI"
)
