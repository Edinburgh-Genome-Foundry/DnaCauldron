import numpy as np
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation


def random_dna_sequence(length, probas=None, seed=None):
    """Return a random DNA sequence ("ATGGCGT...") with the specified length.

    Parameters
    ----------

    length
      Length of the DNA sequence.

    proba
      Frequencies for the different nucleotides, for instance
      ``probas={"A":0.2, "T":0.3, "G":0.3, "C":0.2}``.
      If not specified, all nucleotides are equiprobable (p=0.25).

    seed
      The seed to feed to the random number generator. When a seed is provided
      the random results depend deterministically on the seed, thus enabling
      reproducibility

    """
    if seed is not None:
        np.random.seed(seed)
    if probas is None:
        sequence = np.random.choice(list("ATCG"), length)
    else:
        bases, probas = zip(*probas.items())
        sequence = np.random.choice(bases, length, p=probas)
    return "".join(sequence)


def load_genbank(filename, linear=True, annotation=None):
    record = SeqIO.read(filename, "genbank")
    record.linear = linear
    return record


def annotate_record(seqrecord, location="full", feature_type="source",
                    margin=0, **qualifiers):
    if location == "full":
        location = (margin, len(seqrecord)-margin)

    strand = location[2] if len(location)==3 else 1
    seqrecord.features.append(
        SeqFeature(
            FeatureLocation(location[0], location[1], strand),
            qualifiers=qualifiers,
            type=feature_type
        )
    )
