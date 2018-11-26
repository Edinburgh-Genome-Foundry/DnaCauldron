import os
import pytest
import matplotlib
matplotlib.use("Agg")
import dnacauldron as dc
import flametree


records_dict = {
    name: dc.load_record(
        os.path.join('tests', 'data', 'assemblies', name + '.gb'),
        id=name, linear=False
    )
    for name in ("partA", "partA2", "partB", "partB2", "partC", "receptor",
                 "connector_A2C")
}

def test_filters():
    to_record = dc.sequence_to_biopython_record

    filter1 = dc.NoRestrictionSiteFilter('BsmBI')
    assert filter1(to_record('ATGATGATG'))
    assert not filter1(to_record('ACGTCTCTTG'))
    assert not filter1(to_record('ACGGAGACGG'))

    filter2 = dc.TextSearchFilter('GFP', is_forbidden=True)
    record = to_record(dc.random_dna_sequence(50))
    assert filter2(record)
    dc.annotate_record(record, location=(20, 40), label="Here's some GFP!")
    assert not filter2(record)


def test_single_assembly(tmpdir):
    parts_files = [records_dict[part]
                   for part in ["partA", "partB", "partC", 'receptor']]
    dc.single_assembly(parts=parts_files, enzyme="BsmBI",
                       outfile=os.path.join(str(tmpdir), "final_sequence.gb"))

def test_single_assembly_with_skipped_part(tmpdir):
    parts_files = [records_dict[part]
                   for part in ["partA", "connector_A2C", "partC", 'receptor']]
    outfile = os.path.join(str(tmpdir), "final_sequence.gb")
    with pytest.raises(dc.AssemblyError) as exc:
        dc.single_assembly(parts=parts_files, enzyme="BsmBI", outfile=outfile)
    assert "0 assemblies" in str(exc.value)

def test_autoselect_enzyme():
    parts = [records_dict[name] for name in ["partA", "partB", "partC"]]
    selected = dc.autoselect_enzyme(parts, enzymes=["BsaI", "BsmBI", "BbsI"])
    assert selected == "BsmBI"

def test_single_assembly_with_wrong_enzyme(tmpdir):
    parts_files = [records_dict[part]
                   for part in ["partA", "partB", "partC"]]
    with pytest.raises(dc.AssemblyError) as exc:
        dc.single_assembly(
            parts=parts_files, enzyme="BsaI",
            outfile=os.path.join(str(tmpdir), "final_sequence.gb")
        )
    assert "0 assemblies" in str(exc.value)

def test_combinatorial_assembly(tmpdir):

    enzyme = "BsmBI"
    part_names = ["partA", "partB", "partC", "partA2", "partB2", "receptor"]
    parts = [records_dict[name] for name in part_names]
    mix = dc.RestrictionLigationMix(parts, enzyme)
    filters = [dc.NoRestrictionSiteFilter(enzyme)]
    assemblies = mix.compute_circular_assemblies(seqrecord_filters=filters)
    assemblies = list(assemblies)
    assert len(assemblies) == 4
    for i, assembly in enumerate(assemblies):
        filepath = os.path.join(str(tmpdir), "%03d.gb" % i)
        dc.write_record(assembly, filepath, "genbank")

def test_swap_donor_vector_part():
    for part_names in [("partA", "partA2"), ("partB", "partB2")]:
        donor, insert = [records_dict[name] for name in part_names]
        _ = dc.swap_donor_vector_part(donor, insert, enzyme='BsmBI')


def test_autoselect_connectors():
    data_root = flametree.file_tree(".").tests.data.select_connectors
    parts = [
        dc.load_record(f._path, linear=False, id=f._name_no_extension)
        for f in data_root.parts_missing_connectors._all_files
        if f._extension == "gb"
    ]
    connectors = [
            dc.load_record(f._path, linear=False, id=f._name_no_extension)
        for f in data_root.connectors._all_files
        if f._extension == "gb"
    ]
    mix = dc.RestrictionLigationMix(parts, enzyme='BsmBI')
    selected_connectors = mix.autoselect_connectors(connectors)
    assert sorted([c.id for c in selected_connectors]) == sorted([
        "conn J-K", "conn L-N", "conn R-W", "conn W-Z", "conn_a_b", "conn D-F"
    ])


def test_full_report(tmpdir):
    part_names = ["partA", "partA2", "partB", "partB2", "partC", "receptor",
                  "connector_A2C"]
    parts = [records_dict[name] for name in part_names]
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

def test_full_assembly_plan_report(tmpdir):
    constructs = {
        "C1": ["receptor", "partA", "partB", "partC"],
        "C2": ["receptor", "partA2", "partB2", "partC"],
        "C3": ["receptor", "partA", "partA2", "partB", "partC"]
    }
    part_names = set([p for parts in constructs.values() for p in parts])
    parts = {
        name: records_dict[name]
        for name in part_names
    }
    plan = [
        (construct, [parts[p] for p in construct_parts])
        for construct, construct_parts in constructs.items()
    ]
    errors, zip_data = dc.full_assembly_plan_report(plan, "@memory",
        enzyme="BsmBI", assert_single_assemblies=False, fail_silently=True)
    assert errors == []

def test_random_constructs_generator():
    enzyme = "BsmBI"
    part_names = ["partA", "partA2", "partB", "partB2", "partC", "receptor",
                  "connector_A2C"]
    parts = [records_dict[name] for name in part_names]
    mix = dc.RestrictionLigationMix(parts, enzyme)
    circular_assemblies = mix.compute_circular_assemblies(randomize=True)
    assembly_list = list(zip([1, 2, 3], circular_assemblies))
    assert len(assembly_list) == 3

def test_insert_parts_on_backbones(tmpdir):
    backbone_names = ['partA2', 'partB2', 'partC', 'receptor']
    backbone_records = [records_dict[name] for name in backbone_names]
    part_records = [records_dict[name] for name in ['partA', 'partB']]

    choices = resultA, resultB = dc.insert_parts_on_backbones(
        part_records, backbone_records, process_parts_with_backbone=True)
    assert resultA.already_on_backbone
    assert resultB.already_on_backbone
    assert resultA.backbone_record.id == 'partA2'
    assert resultB.backbone_record.id == 'partB2'
    dataframe = dc.BackboneChoice.list_to_infos_spreadsheet(choices)
    dataframe.to_excel(os.path.join(str(tmpdir), 'summary.xls'), index=False)
    dc.BackboneChoice.write_final_records(choices, str(tmpdir))

def test_list_overhangs():
    record = records_dict['partA']
    assert dc.utils.list_overhangs([record]) == ['ATTG', 'GGCT']

def test_substitute_overhangs():
    record = records_dict['partA']
    assert dc.utils.list_overhangs([record]) == ['ATTG', 'GGCT']
    new_record = dc.utils.substitute_overhangs(record, {'ATTG': 'ATAA'})
    assert dc.utils.list_overhangs([new_record]) == ['ATAA', 'GGCT']
    new_record = dc.utils.substitute_overhangs(record, {'ATTG': 'ATAA'})
    assert dc.utils.list_overhangs([new_record]) == ['ATAA', 'GGCT']
    new_record = dc.utils.substitute_overhangs(record, {'ATTG': 'ATAA'},
                                               return_linear_parts=True)
    assert str(new_record.seq[:12]) == "CGTCTCAATAAT"
    