"""Useful functions built on top of the DnaCauldron classes to simplify the
most common operations."""

from Bio import Restriction


from ..AssemblyMix import (RestrictionLigationMix, AssemblyError,
                           FragmentSetContainsPartsFilter)
from ..StickyEndsSeq import StickyEnd
from ..tools import (sequence_to_biopython_record, reverse_complement,
                     annotate_record, load_record, write_record)

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
            part = load_record(part, linear=False, id=name)
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
        write_record(assembly, outfile, "genbank")
    return assembly

def complement_parts(parts, candidates_parts, enzyme='autoselect'):
    parts = list(parts)
    complement_parts = list(candidates_parts)
    if enzyme == 'autoselect':
        enzyme = autoselect_enzyme(parts + complement_parts)
    mix = RestrictionLigationMix(parts, enzyme=enzyme)
    return mix.autoselect_connectors(complement_parts)

def get_overhangs_from_record(rec, with_locations=False):
    """Return a least of the (probable) overhangs used building the construct
    """
    def is_overhang(h):
        return (len(h) == 4) and (set(h) <= set("ATGC"))
    if isinstance(rec, str):
        rec = load_record(rec)
    rec.seq = rec.seq.upper()
    overhangs = [
        (f.location.start, "".join(f.qualifiers.get("label", "")))
        for f in sorted(rec.features,
            key=lambda f: 0 if (f.location is None) else f.location.start)
        if f.type == "homology"
    ]
    overhangs = [(start, o) for start, o in overhangs if is_overhang(o)]
    if overhangs == []:
        # try to identify overhangs another way
        part_locs = []
        for f in rec.features:
            if f.type == "misc_feature":
                note = "".join(f.qualifiers.get("note", ""))
                if note.startswith("From "):
                    part_locs.append((int(f.location.start),
                                      int(f.location.end)))
        part_locs = [(0, 0)] + sorted(part_locs) + [(len(rec), len(rec))]
        seq = str(rec.seq)
        inter_parts = [
            (f1[1], seq[f1[1]:f2[0]])
            for f1, f2 in zip(part_locs, part_locs[1:])
        ]
        overhangs = [(start, o) for start, o in inter_parts if is_overhang(o)]
    if with_locations:
        return overhangs
    else:
        return [o for start, o in overhangs]

def substitute_overhangs(record, substitutions, enzyme='auto',
                         return_linear_parts=False):
    """Replace the record's subsequence that corresponds to overhangs.

    This is practical to change the position of a part in a Type-2S
    assembly standard

    Examples
    ----------

    >>> new_record = replace_overhangs(record, {'ATGC': 'CTCG'})

    Parameters
    ----------

    record
      A Biopython record whose internal sequence needs to be replaced
    
    substitutions
      A dict {overhang: new_overhang} of which overhangs must be replaced
    
    enzyme
      Either 'BsmBI', 'BsaI', etc. or just "auto" for automatic selection.
    """
    if enzyme == 'auto':
        enzyme = autoselect_enzyme([record])
    mix = RestrictionLigationMix([record], enzyme=enzyme)
    fragments = [f for f in mix.fragments if not f.is_reverse]
    for fragment in fragments:
        left = fragment.seq.left_end
        if str(left).upper() in substitutions:
            end = StickyEnd(substitutions[str(left).upper()], left.strand)
            fragment.seq.left_end = end
        right = fragment.seq.right_end
        if str(right).upper() in substitutions:
            end = StickyEnd(substitutions[str(right).upper()], right.strand)
            fragment.seq.right_end = end
    new_mix = RestrictionLigationMix(
        fragments=fragments, enzyme=enzyme, fragments_filters=())
    if return_linear_parts:
        fragment = [f for f in new_mix.filtered_fragments
                    if not f.is_reverse][0]
        site = sequence_to_biopython_record(mix.enzyme.site)
        annotate_record(site, label='%s' % enzyme)
        rev_site = site.reverse_complement()
        annotate_record(site, label='%s' % enzyme)
        left_end = sequence_to_biopython_record(str(fragment.seq.left_end))
        right_end = sequence_to_biopython_record(str(fragment.seq.right_end))
        for end in left_end, right_end:
            annotate_record(end, (0, len(end), 0), label='overhang')
        return site + 'A' + left_end + fragment + right_end + 'A' + rev_site 
    else:
        return list(new_mix.compute_circular_assemblies())[0]

def list_overhangs(records, enzyme='auto', parts_only=True):
    """List all overhangs created by restriction in the provided records.
    
    Warning: only overhangs on non-reversed fragments are returned, not
    their reverse-complement.

    Parameters
    ----------

    records
      List of records
    
    enzyme
      Either 'BsmBI', 'BsaI', etc. or just "auto" for automatic selection.
    
    parts_only
      If true, overhangs created by restriction which are not on a part
      (so for instance inside a backbone) will be ignored.
    """
    if enzyme == 'auto':
        enzyme = autoselect_enzyme(records)
    mix = RestrictionLigationMix(records, enzyme=enzyme)
    return mix.list_overhangs(filtered_fragments_only=parts_only)