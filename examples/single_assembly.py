from dnacauldron import single_assembly
import os
os.chdir("data")
single_assembly(parts=["partA.gb", "partB.gb", "partC.gb"],
                receptor="receptor.gb", enzyme="BsmBI",
                outfile=os.path.join("..", "output_data", "final_sequence.gb"))
