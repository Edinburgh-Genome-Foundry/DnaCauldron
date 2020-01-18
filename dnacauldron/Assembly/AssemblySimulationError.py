from ..tools import format_value_for_spreadsheet


class AssemblySimulationError(Exception):
    def __init__(self, message, assembly, suggestion="", data=None):
        self.message = message
        self.assembly = assembly
        self.suggestion = suggestion
        self.data = data or {}
        super().__init__(message)

    def data_as_string(self):
        data_items = sorted(self.data.items())
        items = [
            "%s: %s" % (k, format_value_for_spreadsheet(v))
            for k, v in data_items
        ]
        return ",".join(items)
