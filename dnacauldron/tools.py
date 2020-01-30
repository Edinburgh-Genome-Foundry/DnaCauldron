def format_value_for_spreadsheet(value):
    """Format a value into a string to be written in a spreadsheet cell."""
    if isinstance(value, (list, tuple)):
        return " & ".join([format_value_for_spreadsheet(v) for v in value])
    return str(value)


def format_data_dicts_records_for_spreadsheet(data_dicts):
    """Apply format_value_for_spreadsheet to every value of a dict"""
    return [
        {
            key: format_value_for_spreadsheet(value)
            for key, value in data_dict.items()
        }
        for data_dict in data_dicts
    ]