import flametree
from dnacauldron import Type2sRestrictionMix, load_record

data_root = flametree.file_tree(".").data.select_connectors

parts = [
    load_record(f._path, topology='circular', id=f._name_no_extension[:15])
    for f in data_root.parts_missing_connectors._all_files
    if f._extension == "gb"
]
connectors = [
    load_record(f._path, topology='circular', id=f._name_no_extension[:15])
    for f in data_root.connectors._all_files
    if f._extension == "gb"
]
mix = Type2sRestrictionMix(parts, enzyme='BsmBI')
selected_connectors = mix.autoselect_connectors(connectors)
print ("Selected connectors: ", ", ".join([c.id for c in selected_connectors]))
