import os

this_directory = os.path.dirname(os.path.realpath(__file__))
parts_folder = os.path.join(this_directory, "parts")
import dnacauldron as dc


def test_combinatorial_type2s():
    repository = dc.SequenceRepository()
    repository.import_records(folder=parts_folder, use_file_names_as_ids=True)
    parts_list = list(repository.collections["parts"])

    # EXPECT A SINGLE CONSTRUCT, GET AN ERROR

    assembly = dc.Type2sRestrictionAssembly(parts_list, expected_constructs=1)
    simulation = assembly.simulate(sequence_repository=repository)
    assert len(simulation.errors) == 1
    assert len(simulation.construct_records) == 5

    # DON'T EXPECT A CERTAIN NUMBER, GET NO ERROR

    assembly = dc.Type2sRestrictionAssembly(
        parts_list, expected_constructs="any_number"
    )
    simulation = assembly.simulate(sequence_repository=repository)
    assert len(simulation.errors) == 0
    assert len(simulation.construct_records) == 5

    # LIMIT THE NUMBER OF CONSTRUCTS, GET LESS CONSTRUCTS, AND A WARNING

    assembly = dc.Type2sRestrictionAssembly(
        parts_list,
        max_constructs=3,
        expected_constructs="any_number",
        expect_no_unused_parts=False,
    )
    simulation = assembly.simulate(sequence_repository=repository)
    assert len(simulation.errors) == 0
    assert len(simulation.warnings) == 1
    assert len(simulation.construct_records) == 3
