from Bio import Restriction
from Bio.Seq import Seq


class NoRestrictionSiteFilter:

    def __init__(self, enzyme_name):
        self.enzyme_name = enzyme_name
        self.enzyme = Restriction.__dict__[enzyme_name]

    def __call__(self, seqrecord):
        linear = seqrecord.linear if hasattr(seqrecord, "linear") else True
        if linear:
            # Shameful hack so that enzyme sites of ensymes cutting outside
            # of the sequence (but have their site inside) will be detected
            seq = "AAAAAA" + Seq(str(seqrecord.seq)) + "AAAAAA"
        else:
            seq = seqrecord.seq
        return (self.enzyme.search(seq, linear=linear) == [])


class NoPatternFilter:

    def __init__(self, pattern):
        self.pattern = pattern

    def __call__(self, seqrecord):
        return seqrecord.seq.find(self.pattern == -1)


class TextSearchFilter:

    def __init__(self, text, is_forbidden=False):
        self.text = text
        self.is_forbidden = is_forbidden

    @staticmethod
    def gather_all_feature_text(feature):
        return " ".join(
            [feature.type] +
            list(map(str, feature.qualifiers.keys())) +
            list(map(str, feature.qualifiers.values()))
        )

    def gather_all_texts(self, seqrecord):
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
