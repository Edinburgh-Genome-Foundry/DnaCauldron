import dnacauldron as dc
import os
import pytest

repo = dc.SequenceRepository()
RECORDS_FOLDER = os.path.join("tests", "data", "assemblies")
repo.import_records(folder=RECORDS_FOLDER, use_file_names_as_ids=True)


def test_single_assembly(tmpdir):
    parts = ["partA", "partB", "partC", "receptor"]
    assembly = dc.Type2sRestrictionAssembly(parts=parts)
    mix, constructs_records = assembly.simulate(sequence_repository=repo)
    assert len(constructs_records) == 1
    assert len(constructs_records[0]) == 8016


def test_single_assembly_with_and_without_skipped_part(tmpdir):
    parts = ["partA", "connector_A2C", "partC", "receptor"]
    asm = dc.Type2sRestrictionAssembly(parts=parts, no_skipped_parts=True)
    mix, construct_records = asm.simulate(sequence_repository=repo)
    assert construct_records == []

    asm = dc.Type2sRestrictionAssembly(parts=parts, no_skipped_parts=False)
    mix, construct_records = asm.simulate(sequence_repository=repo)
    assert len(construct_records) == 1


def test_single_assembly_with_wrong_enzyme(tmpdir):
    assembly = dc.Type2sRestrictionAssembly(
        parts=["partA", "partB", "partC"], enzyme="BsaI"
    )
    mix, constructs_records = assembly.simulate(sequence_repository=repo)
    assert constructs_records == []


def test_combinatorial_assembly(tmpdir):
    parts = ["partA", "partB", "partC", "partA2", "partB2", "receptor"]
    assembly = dc.Type2sRestrictionAssembly(
        parts=parts, no_skipped_parts=False
    )
    mix, constructs_records = assembly.simulate(sequence_repository=repo)
    assert len(constructs_records) == 4
