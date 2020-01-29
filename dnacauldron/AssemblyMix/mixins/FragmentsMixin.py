from ...biotools import annotate_record


class FragmentsMixin:
    """Mixin for AssemblyMix"""


    fragment_annotation_color = "white"

    @property
    def filtered_fragments(self):
        """Return the fragments of the mix passing all the tests

        Generally used to remove fragments containing a restriction site used
        in a Type2S assembly.
        """
        return [
            f
            for f in (self.fragments + self.reverse_fragments)
            if all([fl(f) for fl in self.fragment_filters])
        ]

    def compute_reverse_fragments(self):
        """Precompute self.reverse_fragments.

        This method also marks all "direct" fragments in the mix as
        `fragment.is_reversed=True` and all "reverse" fragments as
        `fragment.is_reversed=False`.
        """
        self.reverse_fragments = []
        for fragment in self.fragments:
            fragment.is_reversed = False
            new_fragment = fragment.reverse_complement()
            new_fragment.is_reversed = True
            new_fragment.reverse_fragment = fragment
            fragment.reverse_fragment = new_fragment
            new_fragment.original_part = fragment.original_part
            self.reverse_fragments.append(new_fragment)

    def annotate_fragment_with_part(self, fragment):
        part = fragment.original_part.id
        if self.annotate_fragments_with_parts:
            annotate_record(
                fragment,
                feature_type="misc_feature",
                source=part,
                indicates_part=True,
                note="From " + part,
                color=self.fragment_annotation_color,
                ApEinfo_fwdcolor=self.fragment_annotation_color,
            )
