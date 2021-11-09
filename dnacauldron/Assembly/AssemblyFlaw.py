from ..tools import format_value_for_spreadsheet


class AssemblyFlaw(Exception):
    """Represents a flaw (error or warning) in an assembly.

    Note that it subclasses Exception which makes it easy to "raise" a flaw,
    but apart from that it is a class of its own.

    Parameters
    ----------

    message
      String explaining the flaw

    assembly
      The Assembly instance on which the flaw was detected.

    suggestion
      Suggestion on how to fix the flaw.

    data
      A dictionary with more data which users could process.
    """

    def __init__(self, message, assembly, suggestion="", data=None):
        self.message = message
        self.assembly = assembly
        self.suggestion = suggestion
        self.data = data or {}
        super().__init__(message)

    def data_as_string(self):
        """Return a comma-separated string of the data, for reports."""
        data_items = sorted(self.data.items())
        items = ["%s: %s" % (k, format_value_for_spreadsheet(v)) for k, v in data_items]
        return ",".join(items)
