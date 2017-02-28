import dnacauldron as dc
import os
parts = [
    dc.load_genbank(os.path.join('data', partfile), linear=False,
                    name=partfile[:-3])
    for partfile in("partA.gb", "partA2.gb", "partB.gb", "partB2.gb",
                    "partC.gb", "receptor.gb")
]
target = os.path.join('../../', 'my_report')
dc.full_assembly_report(parts, target, enzyme="BsmBI",
                        max_assemblies=40, fragments_filters='auto',
                        assemblies_prefix='asm')
