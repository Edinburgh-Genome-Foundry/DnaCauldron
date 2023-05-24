
import dnacauldron as dc
import os
this_directory = os.path.dirname(os.path.realpath(__file__))
records_folder = os.path.join(this_directory, "records")

repo = dc.SequenceRepository()
repo.import_records(folder=records_folder, use_file_names_as_ids=True)

def test_insert_parts_on_backbones(tmpdir):
    backbones = repo.get_records(["partA2", "partB2", "partC", "receptor"])
    parts = repo.get_records(["partA", "partB"])

    choices = resultA, resultB = dc.insert_parts_on_backbones(
        part_records=parts,
        backbone_records=backbones,
        process_parts_with_backbone=True
    )
    assert resultA.already_on_backbone
    assert resultB.already_on_backbone
    assert resultA.backbone_record.id == "partA2"
    assert resultB.backbone_record.id == "partB2"
    dataframe = dc.BackboneChoice.list_to_infos_spreadsheet(choices)
    dataframe.to_excel(os.path.join(str(tmpdir), "summary.xlsx"), index=False)
    dc.BackboneChoice.write_final_records(choices, str(tmpdir))

def test_swap_donor_vector_part():
    for part_names in [("partA", "partA2"), ("partB", "partB2")]:
        donor, insert = repo.get_records(part_names)
        _ = dc.swap_donor_vector_part(donor, insert, enzyme="BsmBI")

def test_list_digestion_overhangs():
    record = repo.get_record("partA")
    assert dc.utils.list_digestion_overhangs([record]) == ["ATTG", "GGCT"]


def test_substitute_overhangs():
    record = repo.get_record("partA")
    assert dc.utils.list_digestion_overhangs([record]) == ["ATTG", "GGCT"]
    new_record = dc.utils.substitute_overhangs(record, {"ATTG": "ATAA"})
    assert dc.utils.list_digestion_overhangs([new_record]) == ["ATAA", "GGCT"]
    new_record = dc.utils.substitute_overhangs(record, {"ATTG": "ATAA"})
    assert dc.utils.list_digestion_overhangs([new_record]) == ["ATAA", "GGCT"]
    new_record = dc.utils.substitute_overhangs(
        record, {"ATTG": "ATAA"}, return_linear_parts=True
    )
    assert str(new_record.seq[:12]) == "CGTCTCAATAAT"

def test_list_overhangs_from_record_annotations():
    record = repo.get_record("construct_1")
    overhangs = dc.utils.list_overhangs_from_record_annotations(record)
    overhangs = sorted(overhangs)
    assert  overhangs == ['AATG', 'CGCT', 'GCTT', 'GGAG', 'GGTA', 'TTCG']