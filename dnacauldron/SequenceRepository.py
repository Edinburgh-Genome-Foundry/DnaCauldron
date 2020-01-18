from .biotools import load_records_from_files, set_record_topology


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


    def __init__(
        self, parts=None, connectors=None, constructs=None, name="repo"
    ):
        def process(records):
            if records is None:
                return {}
            elif isinstance(records, list):
                return {r.id: r for r in records}
            else:
                return records

        self.parts = process(parts)
        self.connectors = process(connectors)
        self.constructs = process(constructs)
        self.name = name

    def contains_part(self, name):
        repos = [self.parts, self.constructs] + list(self.connectors.values())
        return any([name in repo for repo in repos])

    def get_record(self, name):
        for repo in [self.parts, self.connectors, self.constructs]:
            if name in repo:
                return repo[name]
        raise NotInRepositoryError(name, self)

    def get_records(self, names):
        records = []
        not_in_repository = []
        for name in names:
            if self.contains_part(name):
                records.append(self.get_record(name))
            else:
                not_in_repository.append(name)
        if len(not_in_repository):
            raise NotInRepositoryError(not_in_repository, self)
        return [self.get_record(name) for name in names]

    def import_records(
        self,
        files=None,
        folder=None,
        as_connector_collection=None,
        use_file_names_as_ids=True,
        topology="auto",
    ):
        if folder is not None:
            records = load_records_from_files(
                folder=folder, use_file_names_as_ids=use_file_names_as_ids
            )
        elif files is not None:
            records = load_records_from_files(
                files=files,
                use_file_names_as_ids=use_file_names_as_ids,
            )
        else:
            raise ValueError("Provide either ``files`` or ``folder``")
        if topology in ["circular", "linear"]:
            for r in records:
                set_record_topology(r, topology)
        records = {r.id: r for r in records}
        if as_connector_collection is not None:
            self.connectors[as_connector_collection] = records
        else:
            self.parts.update(records)
