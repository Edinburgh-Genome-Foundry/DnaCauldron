"""Useful functions to simplify the most common operations."""

from Bio import SeqIO, Restriction
from .Filter import NoRestrictionSiteFilter
from .AssemblyMix import RestrictionLigationMix, AssemblyError


def autoselect_enzyme(parts, enzymes):
    """Finds the enzyme that the parts were probably meant to be assembled with

    Parameters
    ----------

    parts
      A list of SeqRecord files. They should have a "linear" attribute set to
      True or False, otherwise

    This finds the enzyme that

    """
    def enzyme_fit_score(enzyme_name):
        enz = Restriction.__dict__[enzyme_name]
        return sum([
            abs(2-len(enz.search(part.seq,
                                 linear=part.__dict__.get('linear', False))))
            for part in parts
        ])
    return min(enzymes, key=enzyme_fit_score)

def single_assembly(parts, receptor, outfile=None,
                    enzyme="BsmBI", annotate_homologies=True,
                    mix_class="restriction"):
    """Return the single assembly obtained by assembling together different
    parts on a receptor vector.

    Parameters
    ----------

    parts
      A list of either filenames or Biopython records or parts. They are
      assumed circular (i.e. parts on backbones) but linear parts should
      work too as long as the receptor is circular. To make sure that
      a part will be treated as linear DNA, provide a Biopython record
      for that part, with a ``.linear`` attribute set to true.

    outfile
      Name of a genbank file where to output the result

    enzyme
      Name of the enzyme used for the assembly.
    """

    if mix_class == "restriction":
        mix_class = RestrictionLigationMix

    def load_genbank(filename):
        """Specific loag_genbank flavor"""
        record = SeqIO.read(filename, "genbank")
        record.linear = False
        record.name = filename.split("/")[-1].split(".")[0].lower()
        if isinstance(receptor, str) and (filename == receptor):
            record.name += " (RECEPTOR)"
        return record
    parts_records = [
        load_genbank(part) if isinstance(part, str) else part
        for part in parts + [receptor]
    ]
    biopython_enzyme = Restriction.__dict__[enzyme]
    sites_in_receptor = \
        len(biopython_enzyme.search(parts_records[-1].seq, linear=False))

    def exactly_one_receptor_vector(fragments):
        receptor_fragments = [
            fragment for fragment in fragments
            if "(RECEPTOR)" in fragment.original_construct.name
        ]
        return len(receptor_fragments) == sites_in_receptor - 1

    mix = mix_class(parts_records, enzyme)
    assemblies = mix.compute_circular_assemblies(
        fragments_sets_filters=[exactly_one_receptor_vector],
        annotate_homologies=annotate_homologies
    )
    assemblies = list(assemblies)
    N = len(assemblies)
    if N != 1:
        raise AssemblyError('Found %d assemblies instead of 1 expected' % N)
    assembly = assemblies[0]
    if outfile is not None:
        SeqIO.write(assembly, outfile, "genbank")
    return assembly

def swap_donor_vector_part(donor_vector, insert, enzyme):
    """Return the records obtained by cloning inserts into a donor vector.

    Meant for Type-2S assembly standards only (Golden Gate, etc.)

    This method is meant to quickly go from a linearized sequence of a part
    to a circular vector (the part in its donor vector) by starting from
    an existing donor vector (with same overhangs) and swapping this vector's
    part for the insert of interest.

    Parameters
    ----------
    donor_vector
      Biopython record of a donor vector. must have an insert producing a
      restriction-free fragment

    insert
      Biopython record of a plasmid or a linear DNA sequence containing an
      insert (i.e. a fragment that is cut out)

    enzyme
      The name of the enzyme to use e.g. 'BsmBI', 'BsaI', ...

    """

    mix = RestrictionLigationMix([donor_vector], enzyme=enzyme)
    donor_fragments = [
        f for f in mix.fragments
        if len(mix.enzyme.search('A' + f.seq.to_standard_sequence())) > 0
    ]
    for fr in donor_fragments:
        fr.features = [f for f in fr.features
                       if 'source' not in f.qualifiers]
    assert(len(donor_fragments) == len(mix.fragments) - 1)

    mix = RestrictionLigationMix([insert], enzyme=enzyme)
    insert_fragments = [
        f for f in mix.fragments
        if len(mix.enzyme.search('A' + f.seq.to_standard_sequence())) == 0
    ]
    assert(len(insert_fragments) == 1)
    mix = RestrictionLigationMix(
        fragments=[insert_fragments[0]] + donor_fragments,
        enzyme=enzyme,
        fragments_filters=()
    )
    assemblies = list(mix.compute_circular_assemblies())
    assert (len(assemblies) == 1)
    return assemblies[0]
