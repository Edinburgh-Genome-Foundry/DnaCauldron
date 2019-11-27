import os
from Bio import SeqIO  # for exporting to Genbank
from dnacauldron import full_assembly_report, load_record

parts_dir = os.path.join("data", "assemblies")
parts = [
    load_record(
        os.path.join(parts_dir, filename),
        topology="circular",
        id=filename.split('.')[0],
    )
    for filename in os.listdir(parts_dir)
]
full_assembly_report(
    parts,
    target=os.path.join("output_data", "report"),
    enzyme="BsmBI",
    max_assemblies=40,
    fragments_filters="auto",
    assemblies_prefix="asm",
)
print("Your report is ready at 'output_data/report/'")
