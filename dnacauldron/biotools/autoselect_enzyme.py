from Bio import Restriction
from .record_operations import record_is_linear

type2S_enzymes = ("BsmBI", "BsaI", "BbsI", "AarI", "SapI")

def autoselect_enzyme(parts, enzymes=type2S_enzymes):
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
    def number_of_sites(enzyme, part):
        linear = record_is_linear(part, default=False)
        return len(enzyme.search(part.seq, linear=linear))

    def enzyme_fit_score(enzyme_name):
        enzyme = Restriction.__dict__[enzyme_name]
        return sum([abs(2 - number_of_sites(enzyme, part)) for part in parts])

    return min(enzymes, key=enzyme_fit_score)