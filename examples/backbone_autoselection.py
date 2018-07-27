"""

"""
from dnacauldron import load_record, insert_parts_on_backbones, BackboneChoice
import flametree  # for getting/writing files and folders

root = flametree.file_tree('.')
output_dir = root._dir('output_data')._dir('backbone_autoselection')

# load lists of genbanks: one list for parts, one list for potential backbones
part_records = [
    load_record(f._path, linear=False, id=f._name_no_extension)
    for f in root.data.assemblies._all_files
    if f._name_no_extension in ['partA', 'partB']
]
backbone_records = [
    load_record(f._path, linear=False, id=f._name_no_extension)
    for f in root.data.backbones._all_files
]

choices = insert_parts_on_backbones(
    part_records, backbone_records, process_parts_with_backbone=True)
dataframe = BackboneChoice.list_to_infos_spreadsheet(choices)
dataframe.to_excel(output_dir._file('summary.xls')._path, index=False)
BackboneChoice.write_final_records(choices, output_dir._dir("records")._path)
