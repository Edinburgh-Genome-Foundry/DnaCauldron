import numpy as np
from copy import deepcopy
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
from Bio.SeqFeature import SeqFeature, FeatureLocation
from snapgene_reader import snapgene_file_to_seqrecord
from Bio.Seq import Seq

import os

def complement(dna_sequence):
    """Return the complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"TACGGC"``.

    Uses BioPython for speed.
    """
    return str(Seq(dna_sequence).complement())


def reverse_complement(sequence):
    """Return the reverse-complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"GCCGTA"``.

    Uses BioPython for speed.
    """
    return complement(sequence)[::-1]

def sequence_to_biopython_record(sequence, id='<unknown id>',
                                 name='<unknown name>', features=()):
    """Return a SeqRecord of the sequence, ready to be Genbanked."""
    return SeqRecord(Seq(sequence, alphabet=DNAAlphabet()),
                     id=id, name=name, features=list(features))

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


def load_record(filename, linear=True, id='auto', upperize=True):
    if filename.lower().endswith(("gb", "gbk")):
        record = SeqIO.read(filename, "genbank")
    elif filename.lower().endswith(('fa', 'fasta')):
        record = SeqIO.read(filename, "fasta")
    elif filename.lower().endswith('.dna'):
        record = snapgene_file_to_seqrecord(filename)
    else:
        raise ValueError('Unknown format for file: %s' % filename)
    if upperize:
        record = record.upper()
    record.linear = linear
    if id == 'auto':
        id = record.id
        if id in [None, '', "<unknown id>", '.', ' ']:
            id = os.path.splitext(os.path.basename(filename))[0]
            record.name = id.replace(" ", "_")[:20]
        record.id = id
    elif id is not None:
        record.id = id
        record.name = id.replace(" ", "_")[:20]

    return record


def load_genbank(filename, linear=True, name="unnamed"):
    """Load a genbank file

    Parameters
    ----------

    linear
      Set to True for linear constructs, False for circular constructs

    name
      The name of the record. Should be the name of the part if the record
      represents a part.

    """
    record = SeqIO.read(filename, "genbank")
    record.linear = linear
    record.id = name
    record.name = name.replace(" ", "_")[:20]
    return record


def annotate_record(seqrecord, location="full", feature_type="misc_feature",
                    margin=0, **qualifiers):
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
        location = (margin, len(seqrecord)-margin)

    strand = location[2] if len(location) == 3 else 1
    seqrecord.features.append(
        SeqFeature(
            FeatureLocation(location[0], location[1], strand),
            qualifiers=qualifiers,
            type=feature_type
        )
    )

def write_record(record, target, fmt='genbank'):
    """Write a record as genbank, fasta, etc. via Biopython, with fixes"""
    record = deepcopy(record)
    record.name = record.name[:20]
    if str(record.seq.alphabet.__class__.__name__) != 'DNAAlphabet':
        record.seq.alphabet = DNAAlphabet()
    if hasattr(target, 'open'):
        target = target.open('w')
    SeqIO.write(record, target, fmt)