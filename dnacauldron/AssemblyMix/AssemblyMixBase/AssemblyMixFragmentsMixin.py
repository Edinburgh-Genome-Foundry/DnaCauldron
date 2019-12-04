class AssemblyMixFragmentsMixin:

    @property
    def filtered_fragments(self):
        """Return the fragments of the mix passing all the tests

        Generally used to remove fragments containing a restriction site used
        in a Type2S assembly.
        """
        return [
            f
            for f in (self.fragments + self.reverse_fragments)
            if all([fl(f) for fl in self.fragments_filters])
        ]

    def compute_reverse_fragments(self):
        """Precompute self.reverse_fragments.

        This method also marks all "direct" fragments in the mix as
        `fragment.is_reverse=True` and all "reverse" fragments as
        `fragment.is_reverse=False`.
        """
        self.reverse_fragments = []
        for fragment in self.fragments:
            fragment.is_reverse = False
            new_fragment = fragment.reverse_complement()
            new_fragment.is_reverse = True
            new_fragment.reverse_fragment = fragment
            fragment.reverse_fragment = new_fragment
            new_fragment.original_construct = fragment.original_construct
            self.reverse_fragments.append(new_fragment)