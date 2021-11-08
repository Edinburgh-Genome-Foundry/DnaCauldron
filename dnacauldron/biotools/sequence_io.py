import re
import os
from io import BytesIO, StringIO
from copy import deepcopy
import flametree
from snapgene_reader import snapgene_file_to_seqrecord
from Bio import SeqIO

try:
    # Biopython <1.78
    from Bio.Alphabet import DNAAlphabet

    has_dna_alphabet = True
except ImportError:
    # Biopython >=1.78
    has_dna_alphabet = False
from .record_operations import (
    set_record_topology,
    sequence_to_biopython_record,
)


def string_to_records(string):
    """Convert a string of a fasta, genbank... into a simple ATGC string.

    Can also be used to detect a format.
    """
    matches = re.match("([ATGC][ATGC]*)", string)
    if (matches is not None) and (matches.groups()[0] == string):
        return [sequence_to_biopython_record(string)], "ATGC"

    for fmt in ("fasta", "genbank"):
        try:
            stringio = StringIO(string)
            records = list(SeqIO.parse(stringio, fmt))
            if len(records) > 0:
                return (records, fmt)
        except Exception:
            pass
    try:
        record = snapgene_file_to_seqrecord(filecontent=StringIO(string))
        return [record]
    except Exception:
        pass
    raise ValueError("Invalid sequence format")


def load_record(
    filepath,
    topology="default_to_linear",
    id="auto",
    upperize=True,
    max_name_length=20,
):
    """Return a Biopython record read from a Fasta/Genbank/Snapgene file.

    Parameters
    ----------

    filepath
      Path to a Genbank, Fasta, or Snapgene (.dna) file.

    topology
      Can be "circular", "linear", "default_to_circular" (will default
      to circular if ``annotations['topology']`` is not already set) or
      "default_to_linear".

    id
      Sets the record.id. If "auto", the original record.id is used, and if
      none is set the name of the file (without extension) is used instead.

    upperize
      If true, the sequence will get upperized (recommended in this library,
      as the mix of upper and lower case can cause problems in Biopython's
      enzyme site search).

    max_name_length
      The name of the record will be truncated if too long to avoid Biopython
      exceptions being raised.
    """
    if filepath.lower().endswith(("gb", "gbk")):
        record = SeqIO.read(filepath, "genbank")
    elif filepath.lower().endswith(("fa", "fasta")):
        record = SeqIO.read(filepath, "fasta")
    elif filepath.lower().endswith(".dna"):
        record = snapgene_file_to_seqrecord(filepath)
    else:
        raise ValueError("Unknown format for file: %s" % filepath)
    if upperize:
        record = record.upper()
    set_record_topology(record, topology)
    if id == "auto":
        id = record.id
        if id in [None, "", "<unknown id>", ".", " "]:
            id = os.path.splitext(os.path.basename(filepath))[0]
            id = id.replace(" ", "_")[:max_name_length]
        record.id = id
    elif id is not None:
        record.id = id.replace(" ", "_")[:max_name_length]

    return record


def _load_records_from_zip_file(zip_file, use_file_names_as_ids=False):
    """Return all fasta/genbank/snapgene in a zip as biopython records.

    Each record gets a ``source_file`` attribute from the zip's file name
    without the .zip extension.

    Used via "load_records_from_files".
    """
    zip_file = flametree.file_tree(zip_file)
    records = []
    for f in zip_file._all_files:
        ext = f._extension.lower()
        if ext in ["gb", "gbk", "fa", "dna"]:
            try:
                new_records, fmt = string_to_records(f.read())
                if not isinstance(new_records, list):
                    new_records = [new_records]
            except Exception:
                content_stream = BytesIO(f.read("rb"))
                try:
                    record = snapgene_file_to_seqrecord(fileobject=content_stream)
                    new_records, _ = [record], "snapgene"
                except Exception:
                    raise ValueError("Format not recognized for file " + f._path)

            single_record = len(new_records) == 1
            for i, record in enumerate(new_records):
                name = record.id
                if name in [
                    None,
                    "",
                    "<unknown id>",
                    ".",
                    " ",
                    "<unknown name>",
                ]:
                    number = "" if single_record else ("%04d" % i)
                    name = f._name_no_extension.replace(" ", "_") + number
                record.id = name
                record.id = name
                record.file_name = f._name_no_extension
                if use_file_names_as_ids and single_record:
                    basename = os.path.basename(record.source_file)
                    basename_no_extension = os.path.splitext(basename)[0]
                    record.id = basename_no_extension
            for record in new_records:
                record.source_file = f._path
            records += new_records
    return records


def load_records_from_file(filepath):
    """Autodetect file format and load biopython records from it."""

    with open(filepath, "rb") as f:
        content = f.read()
    try:
        records, fmt = string_to_records(content.decode("utf-8"))
    except Exception:
        try:
            record = snapgene_file_to_seqrecord(fileobject=BytesIO(content))
            records, fmt = [record], "snapgene"
        except Exception:
            raise ValueError("Format not recognized for file " + filepath)
    if not isinstance(records, list):
        records = [records]
    for record in records:
        record.source_file = filepath
    return records, fmt


def load_records_from_files(files=None, folder=None, use_file_names_as_ids=False):
    """Automatically convert files or a folder's content to biopython records.

    Parameters
    ----------

    files
      A list of path to files. A ``folder`` can be provided instead.

    folder
      A path to a folder containing sequence files.

    use_file_names_as_ids
      If True, for every file containing a single record, the file name
      (without extension) will be set as the record's ID.
    """
    if files is not None:
        for file in files:
            if isinstance(file, str) and not os.path.exists(file):
                raise IOError("File %s not found" % file)

    if folder is not None:
        files = [f._path for f in flametree.file_tree(folder)._all_files]
    records = []
    for filepath in files:
        filename = os.path.basename(filepath)
        if filename.lower().endswith("zip"):
            records += _load_records_from_zip_file(
                filepath, use_file_names_as_ids=use_file_names_as_ids
            )
            continue
        recs, fmt = load_records_from_file(filepath)
        single_record = len(recs) == 1
        for i, record in enumerate(recs):
            name_no_extension = "".join(filename.split(".")[:-1])
            name = name_no_extension + ("" if single_record else ("%04d" % i))
            name = name.replace(" ", "_")
            UNKNOWN_IDS = [
                "None",
                "",
                "<unknown id>",
                ".",
                "EXPORTED",
                "<unknown name>",
                "Exported",
            ]

            if has_dna_alphabet:  # Biopython <1.78
                record.seq.alphabet = DNAAlphabet()
            record.annotations["molecule_type"] = "DNA"

            # Sorry for this parts, it took a lot of "whatever works".
            # keep your part names under 20c and pointless, and everything
            # will be good
            if str(record.id).strip() in UNKNOWN_IDS:
                record.id = name
            record.file_name = name_no_extension
            if use_file_names_as_ids and single_record:
                basename = os.path.basename(record.source_file)
                basename_no_extension = os.path.splitext(basename)[0]
                record.id = basename_no_extension
        records += recs
    return records


def write_record(record, target, fmt="genbank"):
    """Write a record as genbank, fasta, etc. via Biopython, with fixes."""
    record = deepcopy(record)
    record.id = record.id[:20]

    if has_dna_alphabet:  # Biopython <1.78
        if str(record.seq.alphabet.__class__.__name__) != "DNAAlphabet":
            record.seq.alphabet = DNAAlphabet()
    record.annotations["molecule_type"] = "DNA"

    if hasattr(target, "open"):
        target = target.open("w")
    SeqIO.write(record, target, fmt)
