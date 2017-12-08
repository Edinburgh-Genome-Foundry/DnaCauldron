import flametree # for getting/writing files and folders
from Bio import SeqIO  # for exporting to Genbank
from dnacauldron import full_assembly_report, load_genbank
root = flametree.file_tree(".")
parts = [
    load_genbank(f._path, linear=False, name=f._name_no_extension)
    for f in root.data.assemblies._all_files
]
target = root._dir('output_data')._dir('report')._path
full_assembly_report(parts, target=target, enzyme="BsmBI",
                     max_assemblies=40, fragments_filters='auto',
                     assemblies_prefix='asm')
print ("Your report is ready at 'output_data/report/'")
