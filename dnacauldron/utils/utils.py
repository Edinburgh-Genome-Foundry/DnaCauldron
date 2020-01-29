"""Useful functions built on top of the DnaCauldron classes to simplify the
most common operations."""
import os

from Bio import Restriction


from ..AssemblyMix import (
    RestrictionLigationMix,
    generate_type2s_restriction_mix
)
from ..Filter import FragmentSetContainsPartsFilter, NoRestrictionSiteFilter
from ..Fragment.StickyEndFragment import StickyEnd
from ..biotools import (
    sequence_to_biopython_record,
    annotate_record,
    load_record,
    autoselect_enzyme,
)

def list_overhangs_from_record_annotations(rec, with_locations=False):
    """Return a least of the (probable) overhangs used building the construct
    """

    def is_overhang(h):
        return (len(h) == 4) and (set(h) <= set("ATGC"))

    if isinstance(rec, str):
        rec = load_record(rec)
    rec.seq = rec.seq.upper()
    overhangs = [
        (f.location.start, "".join(f.qualifiers.get("label", "")))
        for f in sorted(
            rec.features,
            key=lambda f: 0 if (f.location is None) else f.location.start,
        )
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
                    part_locs.append(
                        (int(f.location.start), int(f.location.end))
                    )
        part_locs = [(0, 0)] + sorted(part_locs) + [(len(rec), len(rec))]
        seq = str(rec.seq)
        inter_parts = [
            (f1[1], seq[f1[1] : f2[0]])
            for f1, f2 in zip(part_locs, part_locs[1:])
        ]
        overhangs = [(start, o) for start, o in inter_parts if is_overhang(o)]
    if with_locations:
        return overhangs
    else:
        return [o for start, o in overhangs]


def substitute_overhangs(
    record, substitutions, enzyme="auto", return_linear_parts=False
):
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
    if enzyme == "auto":
        enzyme = autoselect_enzyme([record])
    mix = generate_type2s_restriction_mix(parts=[record], enzyme=enzyme)
    fragments = [f for f in mix.fragments if not f.is_reversed]
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
        fragments=fragments, enzymes=[enzyme], fragment_filters=()
    )
    if return_linear_parts:
        fragment = [f for f in new_mix.filtered_fragments if not f.is_reversed][
            0
        ]
        site = sequence_to_biopython_record(mix.enzymes[0].site)
        annotate_record(site, label="%s" % enzyme)
        rev_site = site.reverse_complement()
        annotate_record(site, label="%s" % enzyme)
        left_end = sequence_to_biopython_record(str(fragment.seq.left_end))
        right_end = sequence_to_biopython_record(str(fragment.seq.right_end))
        for end in left_end, right_end:
            annotate_record(end, (0, len(end), 0), label="overhang")
        return site + "A" + left_end + fragment + right_end + "A" + rev_site
    else:
        return list(new_mix.compute_circular_assemblies())[0]


def list_digestion_overhangs(records, enzyme="auto", parts_only=True):
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
    if enzyme == "auto":
        enzyme = autoselect_enzyme(records)
    mix = generate_type2s_restriction_mix(parts=records, enzyme=enzyme)
    return mix.list_overhangs(filtered_fragments_only=parts_only)
