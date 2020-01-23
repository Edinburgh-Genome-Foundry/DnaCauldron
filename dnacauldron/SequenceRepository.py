from .biotools import (
    load_records_from_files,
    set_record_topology,
    sequence_to_biopython_record,
)


class NotInRepositoryError(Exception):
    def __init__(self, parts, repository):
        self.parts = parts
        self.repository = repository
        parts_list = ", ".join(parts)
        if len(parts_list) > 150:
            parts_list = parts_list[:150] + "..."
        parts = "Part%s %s" % ("s" if len(parts) > 1 else "", parts_list)
        repo_name = (" in " + repository.name) if repository.name else ""
        message = parts + " not found" + repo_name
        super().__init__(message)


class SequenceRepository:
    """Sequence repositories store and provide sequence records.
    
    Sequence records are stored in the repository as Biopython SeqRecords.
    They are separated into attributes ``parts`` (main components of an
    assembly), ``constructs``
    """

    def __init__(self, collections=None, name="repo"):
        self.collections = {}
        self.name = name

    def add_record(self, record, collection="parts"):
        if collection not in self.collections:
            self.collections[collection] = {}
        self.collections[collection][record.id] = record

    def add_records(self, records, collection="parts"):

        if len(records) == 0:
            return
        if isinstance(records[0], (tuple, list)):
            records = [
                sequence_to_biopython_record(_record, id=_id)
                for _record, _id in records
            ]
        for record in records:
            self.add_record(record, collection=collection)

    def contains_record(self, name):
        collections = self.collections.values()
        return any(name in collection for collection in collections)

    def get_record(self, name):
        for collection in self.collections.values():
            if name in collection:
                return collection[name]
        raise NotInRepositoryError([name], self)

    def get_records(self, names):
        records = []
        not_in_repository = []
        for name in names:
            if self.contains_record(name):
                records.append(self.get_record(name))
            else:
                not_in_repository.append(name)
        if len(not_in_repository):
            raise NotInRepositoryError(not_in_repository, self)
        return records

    def import_records(
        self,
        files=None,
        folder=None,
        collection="parts",
        use_file_names_as_ids=True,
        topology="auto",
    ):
        if folder is not None:
            records = load_records_from_files(
                folder=folder, use_file_names_as_ids=use_file_names_as_ids
            )
        elif files is not None:
            records = load_records_from_files(
                files=files, use_file_names_as_ids=use_file_names_as_ids,
            )
        else:
            raise ValueError("Provide either ``files`` or ``folder``")
        if topology in ["circular", "linear"]:
            for r in records:
                set_record_topology(r, topology)

        self.add_records(records, collection=collection)

    def get_part_names_by_collection(self, format="dict"):
        """Format: dict or string"""
        result = {
            collection_name: list(parts.keys())
            for collection_name, parts in self.collections.items()
        }
        if format == "dict":
            return result
        else:
            return "\n".join(
                "\n".join([name] + ["- " + part for part in sorted(parts)])
                for name, parts in result.items()
            )
