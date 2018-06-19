"""Useful functions built on top of the DnaCauldron classes to simplify the
most common operations."""

from Bio import SeqIO, Restriction


from ..AssemblyMix import (RestrictionLigationMix, AssemblyError,
                           FragmentSetContainsPartsFilter)
from ..tools import reverse_complement, load_genbank

def autoselect_enzyme(parts, enzymes=('BsmBI', 'BsaI', 'BbsI')):
    """Finds the enzyme that the parts were probably meant to be assembled with

    Parameters
    ----------

    parts
      A list of SeqRecord files. They should have a "linear" attribute set to
      True or False, otherwise

    Returns
    --------
    The enzyme that has as near as possible as exactly 2 sites in the different
    constructs.
    """
    def enzyme_fit_score(enzyme_name):
        enz = Restriction.__dict__[enzyme_name]
        return sum([
            abs(2-len(enz.search(part.seq,
                                 linear=part.__dict__.get('linear', False))))
            for part in parts
        ])
    return min(enzymes, key=enzyme_fit_score)


def single_assembly(parts, outfile=None, enzyme="autoselect",
                    annotate_homologies=True, mix_class="restriction"):
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

    part_records = []
    for part in parts:
        if isinstance(part, str):
            name = part.split("/")[-1].split(".")[0].lower()
            part = load_genbank(part, linear=False, name=name)
        part_records.append(part)
    if enzyme == 'autoselect':
        enzyme = autoselect_enzyme(parts, ['BsmBI', 'BsaI', 'BbsI'])

    mix = mix_class(part_records, enzyme)
    part_names = [p.id for p in part_records]
    assemblies = mix.compute_circular_assemblies(
        annotate_homologies=annotate_homologies,
        fragments_sets_filters=(FragmentSetContainsPartsFilter(part_names),)
    )
    first_assemblies = list(zip(assemblies, [1, 2]))
    N = len(first_assemblies)
    if N != 1:
        raise AssemblyError('Found %d assemblies instead of 1 expected' % N)

    assembly = first_assemblies[0][0]
    if outfile is not None:
        SeqIO.write(assembly, outfile, "genbank")
    return assembly

def complement_parts(parts, candidates_parts, enzyme='autoselect'):
    parts = list(parts)
    complement_parts = list(candidates_parts)
    if enzyme == 'autoselect':
        enzyme = autoselect_enzyme(parts + complement_parts)
    mix = RestrictionLigationMix(parts, enzyme=enzyme)
    return mix.autoselect_connectors(complement_parts)
