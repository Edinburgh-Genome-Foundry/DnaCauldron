import Levenshtein
import numpy as np


class HomologyChecker:

    tm_dict = {"A": 2, "T": 2, "G": 4, "C": 4}

    def __init__(
        self, min_size=15, max_size=80, min_tm=0, max_tm=None, max_distance=0
    ):
        """Class to define which homologies are acceptable.

        This class is for instance used to detect and validate end-homologies
        for Gibson Assembly.

        Melting temperatures are currently computed using the naive but
        efficient technique of A/T = 2C, G/C= 4C.

        Parameters
        ==========

        min_size, max_size
          Acceptable size, in nucleotides, for an homology. The upper bound
          helps reduce computational times
        
        min_tm, max_tm
          Acceptable melting temperature range in Celsius for the homology.

        max_distance
          If >0, small 1-, or 2- nucleotide differences (edits, deletion) can
          be accepted in homologies. 
        """
        self.min_size = min_size
        self.max_size = max_size
        self.min_tm = min_tm
        self.max_tm = max_tm
        self.max_distance = max_distance

    def compute_tm(self, sequence):
        """Compute the melting temp. of a sequence"""
        return sum([self.tm_dict[c] for c in sequence], 0)

    def sequence_to_string(self, seq):
        if isinstance(seq, str):
            return seq
        if hasattr(seq, "seq"):
            return str(seq.seq)
        return str(seq)

    def find_end_homologies(self, seq1, seq2):
        """Finds an homology between seq1's end and seq2's start.
        
        Return the size of the homology, or 0 if no homologies found.
        """
        seq1 = self.sequence_to_string(seq1)
        seq2 = self.sequence_to_string(seq2)
        max_size = np.min([self.max_size, len(seq1), len(seq2)])
        for homology_size in range(self.min_size, max_size):
            subseq1 = seq1[-homology_size:]
            subseq2 = seq2[:homology_size]
            if self.check_homology(subseq1, subseq2):
                return homology_size
        return 0

    def check_homology(self, sequence, other_sequence=None):
        """Return whether there is an acceptable full-sequence homology between
        two sequences."""
        sequence = self.sequence_to_string(sequence)

        if other_sequence is not None:
            other_sequence = self.sequence_to_string(other_sequence)
            if self.max_distance == 0:
                if sequence != other_sequence:
                    return False
            else:
                distance = Levenshtein.distance(sequence, other_sequence)
                if distance > self.max_distance:
                    return False
        if len(sequence) < self.min_size:
            return False
        if (self.max_size is not None) and len(sequence) > self.max_size:
            return False
        if self.min_tm > 0 or (self.max_tm is not None):
            tm = self.compute_tm(sequence)
            if tm < self.min_tm:
                return False
            if tm > (self.max_tm or 1e8):
                return False
        return True

    def parameters_as_string(self):
        """Return a string of the parameters for errors and reports."""
        return "%d-%dbp, %.1f-%sC Tm" % (
            self.min_size,
            self.max_size,
            self.min_tm,
            "+" if (self.max_tm is None) else ("%.1f" % self.max_tm),
        )

