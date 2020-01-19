
from copy import copy
from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

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
    sequence, id="<unknown id>", name="same_as_id", features=()
):
    """Return a SeqRecord of the sequence, ready to be Genbanked."""
    return SeqRecord(
        Seq(sequence, alphabet=DNAAlphabet()),
        id=id,
        name=id if name == "same_as_id" else name,
        features=list(features),
    )


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


def crop_record_with_saddling_features(record, start, end, filters=()):
    cropped = record[start:end]

    def is_saddling(f_start, f_end):
        return (f_start < start <= f_end) or (f_start <= end < f_end)

    saddling_features = [
        copy(f)
        for f in record.features
        if all([fl(f) for fl in filters])
        and f.location is not None
        and is_saddling(f.location.start, f.location.end)
    ]
    for f in saddling_features:
        f.location = FeatureLocation(
            start=max(f.location.start, start),
            end=min(f.location.end, end),
            strand=f.location.strand,
        )
        cropped.features.append(f)
    return cropped
