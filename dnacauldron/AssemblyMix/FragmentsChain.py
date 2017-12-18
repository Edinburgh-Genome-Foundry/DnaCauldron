
class FragmentsChain:
    """Class to represent a set of DNA fragments that can assemble into
    a linear or circular construct.

    Parameters
    ----------

    fragments
      A list of fragments that can be assembled into a linear/cicular construct.

    is_standardized
      Indicates whether the fragment is in standardized form, which saves time
      by avoiding to standardize the fragment more than once.

    is_cycle
      Indicates whether the fragments are expected to assemble circularly

    Note
    ----
    Importantly these objects are not meant to be modified inplace as their
    hashes are cached to accelerate computations
    """

    def __init__(self, fragments, is_standardized=False, is_cycle=False,
                 precomputed_hash=None):
        self.fragments = fragments
        self.is_standardized = is_standardized
        self.is_cycle = is_cycle
        self._hash = precomputed_hash

    def reverse_complement(self):
        """Return a reverse-complemented FragmentsChain.

        The chain is made of the inverted list of the fragments's
        rev-complements, which would therefore assemble as the
        reverse-complement of the current Fragments Chain.
        """
        return FragmentsChain([f.reverse_fragment
                               for f in self.fragments][::-1],
                              is_cycle=self.is_cycle)

    def standardized(self):
        """Return a standardized version of the cycle.

        Useful for spotting cycles that may look different but are just
        two representations of a same circular DNA construct.

        For instance, a cycle with fragments A-B-C represents the same
        construct as C-A-B, or even rev(C)-rev(B)-rev(A) (reverse complement).

        The standardization works as follows:

        - If more than half of the fragments in the cycle are
          "reverse complement" fragments, consider the reverse version of
          the cycle. This way all standardized cycles have less than 50%
          rev-complement fragments
        - If the chaing is a cycle it is "rotated" so that the first fragment
          of the cycle is the largest fragment. If there are several fragments
          of same largest size we choose the first one in alphabetical order of
          the sequence.
        """

        if self.is_standardized:
            # Note: return a copy but don't use deepcopy here
            # it's a computing bottleneck
            return FragmentsChain(self.fragments, self.is_standardized,
                                  is_cycle=self.is_cycle,
                                  precomputed_hash=self._hash)

        # If some backbone is detected in the chain, the standardization
        # is done relatively to this backbone, which will be in direct sense
        # and the first part of the chain if the chain is a cycle
        backbones = [
            (i, fragment)
            for i, fragment in enumerate(self.fragments)
            if fragment.original_construct.__dict__.get("is_backbone", False)
        ]
        if len(backbones) == 1:
            backbone_index, backbone = backbones[0]
            if backbone.is_reverse:
                return self.reverse_complement().standardized()
            elif not self.is_cycle:
                return FragmentsChain(self.fragments,
                                      is_standardized=True,
                                      is_cycle=self.is_cycle,
                                      precomputed_hash=self._hash)
            else:
                std_fragments = (self.fragments[backbone_index:] +
                                 self.fragments[:backbone_index])
                return FragmentsChain(std_fragments,
                                      is_standardized=True,
                                      is_cycle=self.is_cycle)

        # If no backbone is detected in the chain, the standardization
        # is done relatively to this backbone, which will be in direct sense
        # and the first part of the chain if the chain is a cycle

        reverse_proportion = (sum(len(f)
                                  for f in self.fragments
                                  if f.is_reverse) /
                              float(sum(len(f) for f in self.fragments)))
        if reverse_proportion == 0.5:
            f1, f2 = ["%s%s%s" % (f.seq.left_end, f.seq, f.seq.right_end)
                      for f in [self.fragments[0], self.fragments[-1]]]
            if f1 > f2:
                std_fragments = self.reverse_complement().fragments
            else:
                std_fragments = self.fragments
        if (reverse_proportion > 0.5):
            std_fragments = self.reverse_complement().fragments
        else:
            std_fragments = self.fragments

        if self.is_cycle:
            sequences = ["%s%s%s" % (f.seq.left_end, f.seq, f.seq.right_end)
                         for f in std_fragments]
            len_sequences = [len(sequence) for sequence in sequences]
            index = min(range(len(sequences)),
                        key=lambda i: (-len_sequences[i], sequences[i]))
            std_fragments = std_fragments[index:] + std_fragments[:index]

        return FragmentsChain(std_fragments, is_standardized=True,
                              is_cycle=self.is_cycle)

    def __hash__(self):
        """The hash of the cycle is the hash of the concatenation of the
        fragments in the standardized version of the cycle."""

        if self._hash is None:
            self._hash = hash("".join([
                "%s%s%s" % (f.seq.left_end, f.seq, f.seq.right_end)
                for f in self.standardized().fragments
            ]))
        return self._hash
