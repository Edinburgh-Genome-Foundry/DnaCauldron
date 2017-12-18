"""
"""

from Bio.Alphabet import DNAAlphabet


from ..tools import annotate_record
from ..StickyEndsSeq import StickyEndsSeqRecord, StickyEndsSeq, StickyEnd
from .Filter import NoRestrictionSiteFilter, TextSearchFilter
from .RestrictionLigationMix import RestrictionLigationMix


class BASICLigationMix(RestrictionLigationMix):

    @staticmethod
    def find_adapter(record):
        for feature in record.features:
            label = feature.qualifiers.get("label", "")
            if isinstance(label, list):
                label = label[0]
            if label == "adapter":
                return (
                    int(feature.location.start),
                    int(feature.location.end),
                    feature.location.strand
                )
        return None

    def fragments_filters(self):
        enzyme_filter = NoRestrictionSiteFilter(str(self.enzyme))
        return [
            lambda frag: (self.find_adapter(frag) or enzyme_filter(frag))
        ]

    def compute_digest(self, construct):

        adapter = self.find_adapter(construct)
        if adapter:
            start, end, strand = adapter
            left_end = StickyEnd(str(construct[:start].seq), strand=1)
            right_end = StickyEnd(str(construct[end:].seq), strand=1)
            sequence = StickyEndsSeq(str(construct[start:end].seq),
                                     left_end=left_end,
                                     right_end=right_end)
            sequence.alphabet = DNAAlphabet()
            record = StickyEndsSeqRecord(seq=sequence)
            annotate_record(record, location=(0, len(sequence), 1),
                            label="adapter")
            return [record]
        else:
            # No feature shows that this is an adapter: use simple restriction
            return RestrictionLigationMix.compute_digest(self, construct)

    @staticmethod
    def assemble_constructs_and_linkers(records_list, enzyme="BsaI"):
        fragments = []
        for linker_left, part, linker_right in records_list:
            linker_left.linear = True
            linker_right.linear = True
            if not isinstance(part, list):
                part = [part]

            for p in part:
                mix = BASICLigationMix([linker_left, p, linker_right],
                                       enzyme="BsaI")
                mix.compute_linear_assemblies
                new_fragment = list(mix.compute_linear_assemblies(
                    fragments_sets_filters=(),
                    min_parts=3,
                    seqrecord_filters=[TextSearchFilter("adapter")],
                    annotate_homologies=False
                ))
                if len(new_fragment) != 1:
                    part_names = str([linker_left.name, p.name,
                                      linker_right.name])
                    raise ValueError(
                        "Something weird happened when trying to assemble "
                        "%s. %d assemblies found" % (
                            part_names, len(new_fragment)))
                new_fragment = new_fragment[0]
                new_fragment.original_construct = p
                fragments.append(new_fragment)
        final_mix = BASICLigationMix(fragments=fragments)
        final_mix.compute_reverse_fragments()
        return final_mix.compute_circular_assemblies()
