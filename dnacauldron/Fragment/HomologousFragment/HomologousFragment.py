from copy import deepcopy
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import DNAAlphabet
from ...biotools import set_record_topology

class HomologousFragment(SeqRecord):

    @staticmethod
    def from_standard_record(biopython_record):
        new_record = deepcopy(biopython_record)
        new_record.__class__ = HomologousFragment
        return new_record

    def circularized(
        self,
        homology_checker,
        annotate_homology=False,
        annotation_type="homology"
    ):
        """Return the biopython record obtained by cirularizing the result.

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
        result = double_self[:len(self)]
        set_record_topology(result, "circular")
        return result

    def push_source_features(self, homology_size, side="left"):

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

    def assemble_with(
        self,
        fragment,
        homology_checker,
        annotate_homology=True,
        annotation_type="homology",
    ):
        homology_size = homology_checker.find_end_homologies(self, self)
        if homology_size == 0:
            raise ValueError(
                "Only fragments with end homologies ends can be assembled."
            )
        new_self = deepcopy(self)
        new_fragment = deepcopy(fragment)
        new_self.push_source_features(homology_size, side="left")
        new_fragment.push_source_features(homology_size, side="right")
        feature = SeqFeature(
            FeatureLocation(0, homology_size, 1),
            type=annotation_type,
            qualifiers={"label": "homology"},
        )
        new_fragment.features.append(feature)
        return new_self[:-homology_size] + new_fragment

    @staticmethod
    def assemble(
        fragments,
        homology_checker,
        circularize=False,
        annotate_homologies=False,
    ):
        """Return the (sticky end) record obtained by assembling the fragments.


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
        result.seq.alphabet = DNAAlphabet()
        return result
