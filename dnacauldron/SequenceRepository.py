from .biotools import load_records_from_files


class SequenceRepository:
    def __init__(self, parts=None, connectors=None, constructs=None):
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

    def is_in_repository(self, name):
        return any(
            [
                name in repo
                for repo in [self.parts, self.connectors, self.constructs]
            ]
        )

    def get_record(self, name):
        for repo in [self.parts, self.connectors, self.constructs]:
            if name in repo:
                return repo[name]
        raise ValueError("%s not found in repository." % name)

    def get_records(self, names):
        return [self.get_record(name) for name in names]

    def import_records(
        self,
        file_paths=None,
        folder=None,
        as_connector_collection=None,
        use_file_names_as_ids=True,
    ):
        if folder is not None:
            records = load_records_from_files(
                folder=folder, use_file_names_as_ids=use_file_names_as_ids
            )
        elif file_paths is not None:
            records = load_records_from_files(
                file_paths=file_paths, use_file_names_as_ids=use_file_names_as_ids
            )
        else:
            raise ValueError("Provide either ``file_paths`` or ``folder``")
        records = {r.id: r for r in records}
        if as_connector_collection is not None:
            self.connectors[as_connector_collection] = records
        else:
            self.parts.update(records)
