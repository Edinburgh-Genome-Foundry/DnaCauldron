import os
from Bio import SeqIO # for exporting to Genbank
from dnacauldron import (RestrictionLigationMix, NoRestrictionSiteFilter,
                         load_genbank)
os.chdir("data")
enzyme = "BsmBI"
filters = [NoRestrictionSiteFilter(enzyme)]
parts_files = ["partA.gb", "partA2.gb", "partB.gb", "partB2.gb", "partC.gb",
               "receptor.gb"]
parts = [load_genbank(filename, linear=False) for filename in parts_files]
mix = RestrictionLigationMix(parts, enzyme)
assemblies = mix.compute_circular_assemblies(seqrecord_filters=filters)
for i, assembly in enumerate(assemblies):
    SeqIO.write(assembly, os.path.join("..", "%03d.gb" % i), "genbank")
