class AssemblyBase:

    spreadsheet_import_parameters = ()

    def __init__(
        self,
        name,
        parts,
        connectors_collection=None,
        level=None,
        expected_constructs=1,
        max_constructs=40,
        no_skipped_parts=True,
    ):
        self.name = name
        self.parts = parts
        self.connectors_collection = connectors_collection
        self.expected_constructs = expected_constructs
        self.no_skipped_parts = no_skipped_parts
        self.max_constructs = max_constructs
        self.level = level

    def get_extra_construct_data(self):
        return dict()

    @staticmethod
    def _row_nonempty_cells(row):
        empty_cells_contain = ["-", "nan", "None", ""]
        return [str(e) for e in row if str(e) not in empty_cells_contain]

    @classmethod
    def from_dataframe_row(cls, row):
        line = cls._row_nonempty_cells(row)
        assembly_name = line[0]
        parts, parameters = [], dict()
        for cell in line[1:]:
            for param in cls.spreadsheet_import_parameters:
                prefix = param + ": "
                if cell.startswith(prefix):
                    value = cell[len(prefix) :]
                    try:
                        value = float(value)
                        if int(value) == value:
                            value = int(value)
                    except ValueError:
                        pass
                    parameters[param] = value
                    break
            else:
                parts.append(cell)
        return cls(name=assembly_name, parts=parts, **parameters)

    def get_connectors_records(self, sequence_repository):
        collection = self.connectors_collection
        if collection is None:
            return []
        if isinstance(collection, (tuple, list)):
            return sequence_repository.get_records(collection)
        # last case: the collection is a string (collection name)
        collections = sequence_repository.connectors_collections
        return list(collections[collection].values())

    
