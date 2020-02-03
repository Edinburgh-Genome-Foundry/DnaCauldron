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
