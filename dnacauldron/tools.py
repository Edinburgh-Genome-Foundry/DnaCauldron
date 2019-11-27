import numpy as np
from copy import deepcopy
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
from Bio.SeqFeature import SeqFeature, FeatureLocation
from snapgene_reader import snapgene_file_to_seqrecord
from Bio import Restriction
from Bio.Seq import Seq

import os


def complement(dna_sequence):
    """Return the complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"TACGGC"``.

    Uses BioPython for speed.
    """
    return str(Seq(dna_sequence).complement())


def set_record_topology(record, topology, pass_if_already_set=False):
    record_topology = record.annotations.get("topology", None)
    do_nothing = pass_if_already_set and (record_topology is not None)
    if not do_nothing:
        record.annotations["topology"] = topology


def record_is_linear(record, default=True):
    if "topology" not in record.annotations:
        return default
    else:
        return record.annotations["topology"] == "linear"


def reverse_complement(sequence):
    """Return the reverse-complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"GCCGTA"``.

    Uses BioPython for speed.
    """
    return complement(sequence)[::-1]


def sequence_to_biopython_record(
    sequence, id="<unknown id>", name="<unknown name>", features=()
):
    """Return a SeqRecord of the sequence, ready to be Genbanked."""
    return SeqRecord(
        Seq(sequence, alphabet=DNAAlphabet()),
        id=id,
        name=name,
        features=list(features),
    )


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


def load_record(
    filename,
    topology="auto",
    id="auto",
    upperize=True,
    default_topology="linear",
    max_name_length=20,
):
    if filename.lower().endswith(("gb", "gbk")):
        record = SeqIO.read(filename, "genbank")
    elif filename.lower().endswith(("fa", "fasta")):
        record = SeqIO.read(filename, "fasta")
    elif filename.lower().endswith(".dna"):
        record = snapgene_file_to_seqrecord(filename)
    else:
        raise ValueError("Unknown format for file: %s" % filename)
    if upperize:
        record = record.upper()
    if topology == "auto":
        set_record_topology(record, default_topology, pass_if_already_set=True)
    else:
        set_record_topology(record, topology)
    if id == "auto":
        id = record.id
        if id in [None, "", "<unknown id>", ".", " "]:
            id = os.path.splitext(os.path.basename(filename))[0]
            record.name = id.replace(" ", "_")[:max_name_length]
        record.id = id
    elif id is not None:
        record.id = id
        record.name = id.replace(" ", "_")[:max_name_length]

    return record


def annotate_record(
    seqrecord,
    location="full",
    feature_type="misc_feature",
    margin=0,
    **qualifiers
):
    """Add a feature to a Biopython SeqRecord.

    Parameters
    ----------

    seqrecord
      The biopython seqrecord to be annotated.

    location
      Either (start, end) or (start, end, strand). (strand defaults to +1)

    feature_type
      The type associated with the feature

    margin
      Number of extra bases added on each side of the given location.

    qualifiers
      Dictionnary that will be the Biopython feature's `qualifiers` attribute.
    """
    if location == "full":
        location = (margin, len(seqrecord) - margin)

    strand = location[2] if len(location) == 3 else 1
    seqrecord.features.append(
        SeqFeature(
            FeatureLocation(location[0], location[1], strand),
            qualifiers=qualifiers,
            type=feature_type,
        )
    )


def write_record(record, target, fmt="genbank"):
    """Write a record as genbank, fasta, etc. via Biopython, with fixes"""
    record = deepcopy(record)
    record.name = record.name[:20]
    if str(record.seq.alphabet.__class__.__name__) != "DNAAlphabet":
        record.seq.alphabet = DNAAlphabet()
    if hasattr(target, "open"):
        target = target.open("w")
    SeqIO.write(record, target, fmt)


def autoselect_enzyme(parts, enzymes=("BsmBI", "BsaI", "BbsI", "SapI")):
    """Finds the enzyme that the parts were probably meant to be assembled with

    Parameters
    ----------

    parts
      A list of SeqRecord files. They should have a "linear" attribute set to
      True or False, otherwise

    Returns
    --------
    The enzyme that has as near as possible as exactly 2 sites in the different
    constructs.
    """

    def enzyme_fit_score(enzyme_name):
        enzyme = Restriction.__dict__[enzyme_name]

        def number_of_sites(part):
            linear = record_is_linear(part, default=False)
            return len(enzyme.search(part.seq, linear=linear))

        return sum([abs(2 - number_of_sites(part)) for part in parts])

    return min(enzymes, key=enzyme_fit_score)