import flametree
from dnacauldron import RestrictionLigationMix, load_genbank

data_root = flametree.file_tree(".").data.select_connectors

parts = [
    load_genbank(f._path, linear=False, name=f._name_no_extension[:15])
    for f in data_root.parts_missing_connectors._all_files
    if f._extension == "gb"
]
connectors = [
    load_genbank(f._path, linear=False, name=f._name_no_extension[:15])
    for f in data_root.connectors._all_files
    if f._extension == "gb"
]
mix = RestrictionLigationMix(parts, enzyme='BsmBI')
selected_connectors = mix.autoselect_connectors(connectors)
print ("Selected connectors: ", ", ".join([c.id for c in selected_connectors]))
