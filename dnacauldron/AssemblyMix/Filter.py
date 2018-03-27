"""
Filters applied in methods such as ``mix.compute_circular_assemblies`` in order
to filter out circular assemblies which would have the wrong marker, or
restriction sites of the digestion enzyme (these are unstable)
"""


from Bio import Restriction
from Bio.Seq import Seq


class NoRestrictionSiteFilter:
    """Filters to ignore fragments and final assemblies containing a given
    restriction site
    """

    def __init__(self, enzyme_name):
        self.enzyme_name = enzyme_name
        self.enzyme = Restriction.__dict__[enzyme_name]

    def __call__(self, seqrecord):
        linear = seqrecord.linear if hasattr(seqrecord, "linear") else True
        if linear:
            # Shameful hack so that enzyme sites of enzymes cutting outside
            # of the sequence (but have their site inside) will be detected
            seq = "AAAAAA" + Seq(str(seqrecord.seq)) + "AAAAAA"
        else:
            seq = seqrecord.seq
        return (self.enzyme.search(seq, linear=linear) == [])

    def __repr__(self):
        return ("NoRestriction(%s)" % self.enzyme_name)

    def __str__(self):
        return ("NoRestriction(%s)" % self.enzyme_name)


class NoPatternFilter:
    """Filters to ignore fragments and final assemblies whose DNA sequence
    contains the given pattern.

    The pattern must be an exact sequence of DNA.
    """
    # TODO: regular expressions

    def __init__(self, pattern):
        self.pattern = pattern

    def __call__(self, seqrecord):
        return seqrecord.seq.find(self.pattern == -1)


class TextSearchFilter:
    """Filters to ignore assemblies containing or not containing some text.

    The text will be looked for in every feature of the construct.
    Constructs which do NOT have the text pattern in at least one feature will
    be filtered out, unless ``is_forbidden`` is set to True, at which case
    constructs which DO have the text pattern will be filtered out.
    """

    def __init__(self, text, is_forbidden=False):
        self.text = text
        self.is_forbidden = is_forbidden

    @staticmethod
    def gather_all_feature_text(feature):
        """Return a single string of all text in the feature (+qualifiers)."""
        return " ".join(
            [feature.type] +
            list(map(str, feature.qualifiers.keys())) +
            list(map(str, feature.qualifiers.values()))
        )

    def gather_all_texts(self, seqrecord):
        """Return a single string of all texts in all record features."""
        return " ".join([self.gather_all_feature_text(feature)
                         for feature in seqrecord.features] +
                        list(map(str, seqrecord.annotations)))

    def __call__(self, seqrecord):
        all_texts = self.gather_all_texts(seqrecord)
        text_found = self.text in all_texts
        if self.is_forbidden:
            return not text_found
        else:
            return text_found

class FragmentSetContainsPartsFilter:

    def __init__(self, part_names):
        self.mandatory_part_names = set(part_names)

    def __call__(self, fragments):
        fragments = set([f.original_construct.name for f in fragments])
        return fragments >= self.mandatory_part_names
