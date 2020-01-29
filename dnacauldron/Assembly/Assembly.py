import numpy

from .AssemblyFlaw import AssemblyFlaw


class Assembly:

    spreadsheet_import_parameters = ()

    def __init__(
        self,
        parts,
        name="unnamed_assembly",
        max_constructs=40,
        dependencies=None,
        expected_constructs=1,
        connectors_collection=None,
    ):
        self.name = name
        self.parts = parts
        self.max_constructs = max_constructs
        if dependencies is None:
            dependencies = dict(level=1, depends_on=[], used_in=[])
        self.dependencies = dependencies
        self.expected_constructs = expected_constructs
        self.connectors_collection = connectors_collection

    def get_extra_construct_data(self):
        return dict()

    @staticmethod
    def _row_nonempty_cells(row):
        empty_cells_contain = ["-", "nan", "None", "_", ""]
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
                    parameters[param] = format_string(cell[len(prefix) :])
                    break
            else:
                parts.append(cell)
        return cls(name=assembly_name, parts=parts, **parameters)

    def attribute_ids_to_constructs(self, construct_records):
        n_records = len(construct_records)
        if n_records == 0:
            return
        if n_records == 1:
            construct_records[0].id = self.name
        else:
            digits = int(numpy.ceil(numpy.log10(n_records - 1)))
            for i, record in enumerate(construct_records):
                record.id = "{name}_{num:0{digits}}".format(
                    num=i, digits=digits, name=self.name
                )

    def get_connectors_records(self, sequence_repository):
        collection = self.connectors_collection
        if collection is None:
            return []
        else:
            return list(sequence_repository.collections[collection].values())

    def _detect_constructs_number_error(self, found, flaws_list):
        """Add a new flaw to the list if unexpected constructs_number"""
        expected = self.expected_constructs
        if (expected is not None) and (expected != found):
            flaw = AssemblyFlaw(
                assembly=self,
                message="Wrong number of constructs",
                suggestion="Check assembly or parts design",
                data={"expected_": expected, "found": found},
            )
            flaws_list.append(flaw)

    def _detect_max_constructs_reached(self, found, flaws_list):
        """Add a new flaw to the list if max constructs is reached"""
        max_allowed = self.max_constructs
        if (max_allowed is not None) and (max_allowed <= found):
            flaw = AssemblyFlaw(
                assembly=self,
                message="Wrong number of constructs",
                suggestion="Check assembly or parts design",
                data={"max": max_allowed, "found": found},
            )
            flaws_list.append(flaw)


def format_string(value):
    if value.lower() == "false":
        return False
    if value.lower() == "true":
        return True
    try:
        value = float(value)
        if int(value) == value:
            value = int(value)
    except ValueError:
        pass
    return value
