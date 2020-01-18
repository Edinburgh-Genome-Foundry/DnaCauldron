def format_value_for_spreadsheet(value):
    if isinstance(value, (list, tuple)):
        return " & ".join([format_value_for_spreadsheet(v) for v in value])
    return str(value)


def format_data_dicts_records_for_spreadsheet(data_dicts):
    return [
        {
            key: format_value_for_spreadsheet(value)
            for key, value in data_dict.items()
        }
        for data_dict in data_dicts
    ]