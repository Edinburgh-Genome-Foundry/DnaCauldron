import dnacauldron as dc
import pytest

def test_repository_duplicate_error():

    repo = dc.SequenceRepository()
    records = [("id_1", "ATGC"), ("id_2", "AAAAAT")]
    repo.add_records(records)
    with pytest.raises(dc.RepositoryDuplicateError):
        repo.add_records(records)

def test_not_in_repository_error():

    repo = dc.SequenceRepository()
    records = [("id_1", "ATGC"), ("id_2", "AAAAAT")]
    repo.add_records(records)
    with pytest.raises(dc.NotInRepositoryError):
        repo.get_records(['id_3', 'id_4'])

def test_get_repo_part_names_by_collection():
    repo = dc.SequenceRepository()
    repo.add_records([("id_1", "ATGC")], collection="parts")
    repo.add_records([("id_2", "ATGC")], collection="other")
    
    result = repo.get_part_names_by_collection(format="string")
    assert result == 'parts\n- id_1\nother\n- id_2'