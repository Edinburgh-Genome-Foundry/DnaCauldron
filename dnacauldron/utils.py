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
        fragments_filters=[] if (sites_in_receptor > 2) else
                          [NoRestrictionSiteFilter(enzyme)],
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
