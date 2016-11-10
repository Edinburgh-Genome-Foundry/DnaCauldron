from Bio import SeqIO
from .Filter import NoRestrictionSiteFilter
from AssemblyMix import RestrictionLigationMix

def genbank_files_to_assembly(parts_files, outfile, receptor_name_contains,
                              enzyme="BsmBI"):

    parts_records = []
    for filename in parts_files:
        record = SeqIO.read(filename, "genbank")
        record.linear = False
        record.name = filename.split("/")[-1].split(".")[0].lower()
        parts_records.append(record)

    def exactly_one_receptor_vector(fragments):
        receptor_fragments = [
            fragment for fragment in fragments
            if receptor_name_contains in fragment.original_construct.name
        ]
        return len(receptor_fragments) == 1

    fragments_filters = [NoRestrictionSiteFilter(enzyme)]
    mix = RestrictionLigationMix(parts_records, enzyme)
    assemblies = mix.compute_circular_assemblies(
        fragments_sets_filters=[exactly_one_receptor_vector],
        fragments_filters=fragments_filters,
        annotate_homologies=True
    )
    assemblies = list(assemblies)
    assert len(assemblies) == 1
    assembly = assemblies[0]
    SeqIO.write(assembly, outfile, "genbank")
    return assembly
