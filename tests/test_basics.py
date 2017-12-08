import os
import pytest
import matplotlib
matplotlib.use("Agg")
import dnacauldron as dc
from Bio import SeqIO
import flametree


def test_single_assembly(tmpdir):
    parts_files = [
        os.path.join('tests', 'data', 'assemblies', partfile)
        for partfile in ["partA.gb", "partB.gb", "partC.gb"]
    ]
    receptor_file = os.path.join('tests', 'data', 'assemblies', 'receptor.gb')
    dc.single_assembly(parts=parts_files,
                       receptor=receptor_file, enzyme="BsmBI",
                       outfile=os.path.join(str(tmpdir), "final_sequence.gb"))

def test_autoselect_enzyme():
    parts = [
        dc.load_genbank(os.path.join('tests', 'data', 'assemblies', partfile))
        for partfile in ["partA.gb", "partB.gb", "partC.gb"]
    ]
    selected = dc.autoselect_enzyme(parts, enzymes=["BsaI", "BsmBI", "BbsI"])
    assert selected == "BsmBI"

def test_single_assembly_with_wrong_enzyme(tmpdir):
    parts_files = [
        os.path.join('tests', 'data', 'assemblies', partfile)
        for partfile in ["partA.gb", "partB.gb", "partC.gb"]
    ]
    receptor_file = os.path.join('tests', 'data', 'assemblies', 'receptor.gb')
    with pytest.raises(dc.AssemblyError) as exc:
        dc.single_assembly(
            parts=parts_files, receptor=receptor_file, enzyme="BsaI",
            outfile=os.path.join(str(tmpdir), "final_sequence.gb")
        )
    assert "0 assemblies" in str(exc.value)

def test_combinatorial_assembly(tmpdir):

    enzyme = "BsmBI"
    parts_files = [
        os.path.join('tests', 'data', 'assemblies', partfile)
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

def test_swap_donor_vector_part():
    for part_names in [("partA", "partA2"), ("partB", "partB2")]:
        donor, insert = [
            dc.load_genbank(os.path.join('tests', 'data', 'assemblies',
                                         part_name + '.gb'),
                            linear=False, name=part_name[:-3])
            for part_name in part_names
        ]
        record = dc.swap_donor_vector_part(donor, insert, enzyme='BsmBI')


def test_autoselect_connectors():
    data_root = flametree.file_tree(".").tests.data.select_connectors
    parts = [
        dc.load_genbank(f._path, linear=False, name=f._name_no_extension)
        for f in data_root.parts_missing_connectors._all_files
        if f._extension == "gb"
    ]
    connectors = [
            dc.load_genbank(f._path, linear=False, name=f._name_no_extension)
        for f in data_root.connectors._all_files
        if f._extension == "gb"
    ]
    mix = dc.RestrictionLigationMix(parts, enzyme='BsmBI')
    selected_connectors = mix.autoselect_connectors(connectors)
    assert sorted([c.id for c in selected_connectors]) == sorted([
        "conn J-K", "conn L-N", "conn R-W", "conn W-Z", "conn_a_b", "conn D-F"
    ])


def test_full_report(tmpdir):

    parts = [
        dc.load_genbank(os.path.join('tests', 'data', 'assemblies', partfile),
                        linear=False,
                        name=partfile[:-3])
        for partfile in("partA.gb", "partA2.gb", "partB.gb", "partB2.gb",
                        "partC.gb", "receptor.gb", "connector_A2C.gb")
    ]
    target1 = os.path.join(str(tmpdir), 'my_report')
    target2 = os.path.join(str(tmpdir), 'my_report.zip')
    n, _ = dc.full_assembly_report(parts, '@memory', enzyme="BsmBI",
                                   max_assemblies=40, fragments_filters='auto',
                                   assemblies_prefix='asm')
    assert n == 5
    dc.full_assembly_report(parts, target1, enzyme="BsmBI",
                            max_assemblies=40, fragments_filters='auto',
                            assemblies_prefix='asm')
    dc.full_assembly_report(parts, target2, enzyme="BsmBI",
                            max_assemblies=40, fragments_filters='auto',
                            assemblies_prefix='asm')
    assert os.path.exists(os.path.join(target1, 'assemblies', 'asm_005.gb'))
    assert os.path.exists(target2)
