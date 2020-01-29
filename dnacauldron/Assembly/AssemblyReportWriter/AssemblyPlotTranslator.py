from dna_features_viewer import BiopythonTranslator


class AssemblyPlotTranslator(BiopythonTranslator):
    """Custom theme for plotting GENBANK assemblies with DNA Features Viewer.
    """

    def is_source(self, feature):
        return (feature.type == "misc_feature") and feature.qualifiers.get(
            "source", False
        )

    def compute_feature_color(self, feature):
        if feature:
            if self.is_source(feature):
                return "#ff4c4c"
            else:
                return "#f9edbb"

    def compute_feature_label(self, feature):
        if self.is_source(feature):
            return "".join(feature.qualifiers["source"])
        elif abs(feature.location.end - feature.location.start) > 100:
            label = BiopythonTranslator.compute_feature_label(self, feature)
            return self.shorten_string("".join(label), 30)
        else:
            return None

    @staticmethod
    def shorten_string(string, max_length=30):
        """Truncate and add '...' if the string is too long"""
        suffix = "" if len(string) < max_length else "..."
        return string[:max_length] + suffix
