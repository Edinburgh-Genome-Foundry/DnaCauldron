import os
import dnacauldron as dc

this_directory = os.path.dirname(os.path.realpath(__file__))


def test_autoselect_connectors():
    repo = dc.SequenceRepository()
    repo.import_records(
        collection="parts",
        folder=os.path.join(this_directory, "parts"),
        use_file_names_as_ids=True,
        topology="circular"
    )
    repo.import_records(
        collection="emma_connectors",
        folder=os.path.join(this_directory, "emma_connectors"),
        use_file_names_as_ids=True,
    )
    all_parts = list(repo.collections["parts"].keys())
    assembly = dc.Type2sRestrictionAssembly(
        parts=all_parts, connectors_collection="emma_connectors"
    )
    simulation = assembly.simulate(sequence_repository=repo)
    assert len(simulation.construct_records) == 1
    assert len(simulation.warnings) == 1
    warning = simulation.warnings[0]
    # selected_connectors = mix.autoselect_connectors(connectors)
    assert sorted(warning.data["selected_connectors"]) == [
        "conn_A_B",
        "conn_D-F",
        "conn_J-K",
        "conn_L-N",
        "conn_R-W",
        "conn_W-Z",
    ]
