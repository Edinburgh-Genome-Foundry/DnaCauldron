import os
import dnacauldron.biotools as biotools

this_directory = os.path.dirname(os.path.realpath(__file__))
records_dir = os.path.join(this_directory, "records")


def test_set_record_topology():
    record_path = os.path.join(records_dir, "circular_record.gb")
    record = biotools.load_record(record_path)
    assert record.annotations["topology"] == "circular"
    assert not biotools.record_is_linear(record, default=True)

    # Default only when no topology is set
    biotools.set_record_topology(record, "default_to_circular")
    assert not biotools.record_is_linear(record, default=True)
    record.annotations.pop("topology")
    assert not biotools.record_is_linear(record, default=False)
    biotools.set_record_topology(record, "default_to_linear")
    assert biotools.record_is_linear(record, default=False)


def test_ids_of_multifile_import_are_set_correctly():
    filepaths = [
        os.path.join(records_dir, "circular_record.gb"),
        os.path.join(records_dir, "example_sequences.fa"),
    ]

    records = biotools.load_records_from_files(
        files=filepaths, use_file_names_as_ids=False
    )
    assert sorted([r.id.lower() for r in records]) == [
        "bba_e0040_gfp",
        "expected_construct",
        "frag_0",
        "frag_1",
        "frag_2",
    ]

    records = biotools.load_records_from_files(
        files=filepaths, use_file_names_as_ids=True
    )
    # bba_e0040_gfp now has the file name circular_record as it was a
    # single-file genbank. The other records from the multi-record fasta
    # still use the in-fasta identifier
    assert sorted([r.id.lower() for r in records]) == [
        "circular_record", 
        "expected_construct",
        "frag_0",
        "frag_1",
        "frag_2",
    ]

def test_load_records_from_zip_file():
    zip_path = os.path.join(this_directory, "records.zip")
    records = biotools.load_records_from_files(files=[zip_path])
    assert len(records) == 5