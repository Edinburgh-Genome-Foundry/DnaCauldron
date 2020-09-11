from copy import deepcopy
from Bio.SeqRecord import SeqRecord
from ..Fragment import Fragment

try:
    # Biopython <1.78
    from Bio.Alphabet import DNAAlphabet

    has_dna_alphabet = True
except ImportError:
    # Biopython >=1.78
    has_dna_alphabet = False

from ...biotools import set_record_topology, crop_record_with_saddling_features


class HomologousFragment(Fragment):
    @staticmethod
    def from_biopython_record(biopython_record):
        """Convert a Biopython record into a HomologousFragment (class change).
        """
        new_record = deepcopy(biopython_record)
        new_record.original_part = biopython_record
        new_record.__class__ = HomologousFragment
        return new_record

    def circularized(
        self, homology_checker, annotate_homology=False, annotation_type="homology",
    ):
        """Return the Biopython record obtained by cirularizing the result.

        Only works if the left and right sticky ends are compatible. The
        return is a simple Biopython record where the sticky end has been
        integrated in the sequence.
        """
        double_self = self.assemble_with(
            self,
            homology_checker=homology_checker,
            annotate_homology=True,
            annotation_type="homology",
        )

        def only_parts_indicators(feature):
            return feature.qualifiers.get("indicates_part", False)

        result = crop_record_with_saddling_features(
            record=double_self,
            start=len(self),
            end=len(double_self),
            filters=(only_parts_indicators,),
        )
        result.__class__ = SeqRecord
        set_record_topology(result, "circular")
        return result

    def _push_source_features(self, homology_size, side="left"):
        """Relocate one end of a feature so it won't be cut out when the
        fragment is cropped."""

        if side == "left":
            maximum = len(self) - homology_size
            for f in self.features:
                if f.qualifiers.get("is_source", False):
                    f.location.start = min(f.location.start, maximum)
                    f.location.end = min(f.location.end, maximum)
        if side == "right":
            for f in self.features:
                if f.qualifiers.get("is_source", False):
                    f.location.start = max(f.location.start, homology_size)
                    f.location.end = max(f.location.end, homology_size)

    def will_clip_in_this_order_with(self, other_fragment, homology_checker):
        """Return whether the fragment will assemble with anoter via homology
        recombination.

        homology_checker should be an HomologyChecker instance definining the
        homology conditions.
        """
        homology_size = homology_checker.find_end_homologies(self, other_fragment)
        return homology_size > 0

    def assemble_with(
        self,
        fragment,
        homology_checker,
        annotate_homology=True,
        annotation_type="homology",
    ):
        """Return the fragment resulting from the assembly of this fragment
        with another, in that order.

        Parameters
        ----------

        fragment
          The other parameter to assemble with.

        homology_checker
          An HomologyChecker instance definining the homology conditions.

        annotate_homology
          If true, homologies will have an annotation in the final, predicted
          construct records.
        """
        homology_size = homology_checker.find_end_homologies(self, fragment)
        if homology_size == 0:
            raise ValueError(
                "Only fragments with end homologies ends can be assembled."
            )
        new_self = deepcopy(self)
        new_fragment = deepcopy(fragment)
        new_self._push_source_features(homology_size, side="left")
        new_fragment._push_source_features(homology_size, side="right")
        if annotate_homology:
            feature = self.create_homology_annotation(
                start=0,
                end=homology_size,
                annotation_type=annotation_type,
                label="homology",
            )
            new_fragment.features.append(feature)

        def only_parts_indicators(feature):
            return feature.qualifiers.get("indicates_part", False)

        self_subrecord = crop_record_with_saddling_features(
            record=self,
            start=0,
            end=len(self) - homology_size,
            filters=(only_parts_indicators,),
        )
        assert len(self_subrecord) == len(self[:-homology_size])
        result = self_subrecord + new_fragment
        result.__class__ = HomologousFragment
        return result

    @staticmethod
    def assemble(
        fragments, homology_checker, circularize=False, annotate_homologies=False,
    ):
        """Return the record obtained by assembling the fragments.

        Parameters
        ----------

        fragments
          List of HomologousFragments to assemble.

        homology_checker
          An HomologyChecker instance definining the homology conditions.

        circularize
          True to also assemble the end flanks of the final construct.

        annotate_homologies
           If true, homologies will have an annotation in the final, predicted
           construct records.
        """
        result = fragments[0]
        for fragment in fragments[1:]:
            result = result.assemble_with(
                fragment=fragment,
                homology_checker=homology_checker,
                annotate_homology=annotate_homologies,
            )
        if circularize:
            result = result.circularized(
                annotate_homology=annotate_homologies,
                homology_checker=homology_checker,
            )

        if has_dna_alphabet:  # Biopython <1.78
            result.seq.alphabet = DNAAlphabet()
        result.annotations["molecule_type"] = "DNA"

        return result
