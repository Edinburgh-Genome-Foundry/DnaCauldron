import os
import dnacauldron as dc
from Bio import SeqIO


def test_single_assembly(tmpdir):
    parts_files = [
        os.path.join('tests', 'data', partfile)
        for partfile in ["partA.gb", "partB.gb", "partC.gb"]
    ]
    receptor_file = os.path.join('tests', 'data', 'receptor.gb')
    dc.single_assembly(parts=parts_files,
                       receptor=receptor_file, enzyme="BsmBI",
                       outfile=os.path.join(str(tmpdir), "final_sequence.gb"))


def test_combinatorial_assembly(tmpdir):

    enzyme = "BsmBI"

    parts_files = [
        os.path.join('tests', 'data', partfile)
        for partfile in("partA.gb", "partA2.gb", "partB.gb", "partB2.gb",
                        "partC.gb", "receptor.gb")
    ]
    parts = [
        dc.load_genbank(filename, linear=False)
        for filename in parts_files
    ]
    mix = dc.RestrictionLigationMix(parts, enzyme)
    filters = [dc.NoRestrictionSiteFilter(enzyme)]
    assemblies = mix.compute_circular_assemblies(seqrecord_filters=filters)
    assemblies = list(assemblies)
    assert len(assemblies) == 4
    for i, assembly in enumerate(assemblies):
        filepath = os.path.join(str(tmpdir), "%03d.gb" % i)
        SeqIO.write(assembly, filepath, "genbank")


def test_full_report(tmpdir):

    parts = [
        dc.load_genbank(os.path.join('tests', 'data', partfile), linear=False,
                        name=partfile[:-3])
        for partfile in("partA.gb", "partA2.gb", "partB.gb", "partB2.gb",
                        "partC.gb", "receptor.gb", "connector_A2C.gb")
    ]
    target = os.path.join(str(tmpdir), 'my_report')
    dc.full_assembly_report(parts, target, enzyme="BsmBI",
                            max_assemblies=40, fragments_filters='auto',
                            assemblies_prefix='asm')
    assert os.path.exists(os.path.join(target, 'assemblies', 'asm_05.gb'))
