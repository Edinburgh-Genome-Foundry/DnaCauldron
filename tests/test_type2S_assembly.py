import dnacauldron as dc
import os
import pytest

repo = dc.SequenceRepository()
RECORDS_FOLDER = os.path.join("tests", "data", "assemblies")
repo.import_records(folder=RECORDS_FOLDER, use_file_names_as_ids=True)


# def test_single_assembly(tmpdir):
#     parts = ["partA", "partB", "partC", "receptor"]
#     assembly = dc.Type2sRestrictionAssembly(parts=parts)
#     simulation = assembly.simulate(sequence_repository=repo)
#     assert len(simulation.construct_records) == 1
#     assert len(simulation.construct_records[0]) == 8016


# def test_single_assembly_with_and_without_skipped_part(tmpdir):
#     parts = ["partA", "connector_A2C", "partC", "receptor"]
#     assembly = dc.Type2sRestrictionAssembly(
#         parts=parts, expect_no_unused_parts=True
#     )
#     simulation = assembly.simulate(sequence_repository=repo)
#     assert len(simulation.construct_records) == 1
#     # We expect two errors: 
#     assert len(simulation.errors) == 2

#     assembly = dc.Type2sRestrictionAssembly(
#         parts=parts, expect_no_unused_parts=False
#     )
#     simulation = assembly.simulate(sequence_repository=repo)
#     assert len(simulation.construct_records) == 1
#     assert len(simulation.errors) == 0


def test_single_assembly_with_wrong_enzyme(tmpdir):
    assembly = dc.Type2sRestrictionAssembly(
        parts=["partA", "partB", "partC"], enzyme="BsaI"
    )
    simulation = assembly.simulate(sequence_repository=repo)
    assert simulation.construct_records == []


# def test_combinatorial_assembly(tmpdir):
#     parts = ["partA", "partB", "partC", "partA2", "partB2", "receptor"]
#     assembly = dc.Type2sRestrictionAssembly(
#         parts=parts, expect_no_unused_parts=False
#     )
#     simulation = assembly.simulate(sequence_repository=repo)
#     assert len(simulation.construct_records) == 4
