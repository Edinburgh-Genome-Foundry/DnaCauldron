import os
import pytest
import matplotlib

matplotlib.use("Agg")
import dnacauldron as dc
import flametree

RECORDS_FOLDER = os.path.join("tests", "data", "assemblies")
records_dict = {
    name: dc.load_record(
        os.path.join("tests", "data", "assemblies", name + ".gb"),
        id=name,
        topology="circular",
    )
    for name in (
        "partA",
        "partA2",
        "partB",
        "partB2",
        "partC",
        "receptor",
        "connector_A2C",
    )
}


def test_filters():
    to_record = dc.sequence_to_biopython_record

    filter1 = dc.NoRestrictionSiteFilter("BsmBI")
    assert filter1(to_record("ATGATGATG"))
    assert not filter1(to_record("ACGTCTCTTG"))
    assert not filter1(to_record("ACGGAGACGG"))

    filter2 = dc.TextSearchFilter("GFP", is_forbidden=True)
    record = to_record("ATCGCGTGCGTGCACCACACGT")
    assert filter2(record)
    dc.annotate_record(record, location=(20, 40), label="Here's some GFP!")
    assert not filter2(record)






def test_autoselect_enzyme():
    parts = [records_dict[name] for name in ["partA", "partB", "partC"]]
    selected = dc.autoselect_enzyme(parts, enzymes=["BsaI", "BsmBI", "BbsI"])
    assert selected == "BsmBI"








def test_autoselect_connectors():
    data_root = flametree.file_tree(".").tests.data.select_connectors
    parts = [
        dc.load_record(f._path, topology="circular", id=f._name_no_extension)
        for f in data_root.parts_missing_connectors._all_files
        if f._extension == "gb"
    ]
    connectors = [
        dc.load_record(f._path, topology="circular", id=f._name_no_extension)
        for f in data_root.connectors._all_files
        if f._extension == "gb"
    ]
    mix = dc.Type2sRestrictionMix(parts, enzyme="BsmBI")
    selected_connectors = mix.autoselect_connectors(connectors)
    assert sorted([c.id for c in selected_connectors]) == sorted(
        [
            "conn_J-K",
            "conn_L-N",
            "conn_R-W",
            "conn_W-Z",
            "conn_a_b",
            "conn_D-F",
        ]
    )


def test_full_report(tmpdir):
    part_names = [
        "partA",
        "partA2",
        "partB",
        "partB2",
        "partC",
        "receptor",
        "connector_A2C",
    ]
    parts = [records_dict[name] for name in part_names]
    target1 = os.path.join(str(tmpdir), "my_report")
    target2 = os.path.join(str(tmpdir), "my_report.zip")
    n, _ = dc.full_assembly_report(
        parts,
        "@memory",
        enzyme="BsmBI",
        max_assemblies=40,
        fragment_filters="auto",
        assemblies_prefix="asm",
    )
    assert n == 5
    dc.full_assembly_report(
        parts,
        target1,
        enzyme="BsmBI",
        max_assemblies=40,
        fragment_filters="auto",
        assemblies_prefix="asm",
    )
    dc.full_assembly_report(
        parts,
        target2,
        enzyme="BsmBI",
        max_assemblies=40,
        fragment_filters="auto",
        assemblies_prefix="asm",
    )
    assert os.path.exists(os.path.join(target1, "assemblies", "asm_005.gb"))
    assert os.path.exists(target2)


def test_full_assembly_plan_report(tmpdir):
    constructs = {
        "C1": ["receptor", "partA", "partB", "partC"],
        "C2": ["receptor", "partA2", "partB2", "partC"],
        "C3": ["receptor", "partA", "partA2", "partB", "partC"],
    }
    part_names = set([p for parts in constructs.values() for p in parts])
    parts = {name: records_dict[name] for name in part_names}
    plan = [
        (construct, [parts[p] for p in construct_parts])
        for construct, construct_parts in constructs.items()
    ]
    errors, zip_data = dc.full_assembly_plan_report(
        plan,
        "@memory",
        enzyme="autoselect",
        assert_single_assemblies=False,
        fail_silently=True,
    )
    assert errors == []


def test_random_constructs_generator():
    enzyme = "BsmBI"
    part_names = [
        "partA",
        "partA2",
        "partB",
        "partB2",
        "partC",
        "receptor",
        "connector_A2C",
    ]
    parts = [records_dict[name] for name in part_names]
    mix = dc.Type2sRestrictionMix(parts, enzyme)
    circular_assemblies = mix.compute_circular_assemblies(randomize=True)
    assembly_list = list(zip([1, 2, 3], circular_assemblies))
    assert len(assembly_list) == 3

