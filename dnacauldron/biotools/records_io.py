import re
import os
from io import BytesIO, StringIO
from copy import deepcopy
import flametree
from snapgene_reader import snapgene_file_to_seqrecord
from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
from .records_operations import (
    set_record_topology,
    sequence_to_biopython_record,
)


def string_to_record(string):
    """Convert a string of a fasta, genbank... into a simple ATGC string.

    Can also be used to detect a format.
    """
    matches = re.match("([ATGC][ATGC]*)", string)
    if (matches is not None) and (matches.groups()[0] == string):
        return sequence_to_biopython_record(string), "ATGC"

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
        return record
    except Exception:
        pass
    raise ValueError("Invalid sequence format")


def load_records_from_zip_file(zip_file):
    """Return all fasta/genbank/snapgene in a zip as biopython records."""
    zip_file = flametree.file_tree(zip_file)
    records = []
    for f in zip_file._all_files:
        ext = f._extension.lower()
        if ext in ["gb", "gbk", "fa", "dna"]:
            try:
                new_records, fmt = string_to_record(f.read())
            except Exception:
                content_stream = BytesIO(f.read("rb"))
                try:
                    record = snapgene_file_to_seqrecord(
                        fileobject=content_stream
                    )
                    new_records, _ = [record], "snapgene"
                except Exception:
                    raise ValueError(
                        "Format not recognized for file " + f._path
                    )

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
            for record in new_records:
                record.source_file = f._path
            records += new_records
    return records


def load_records_from_file(filepath):
    """Autodetect file format and load biopython records from it."""

    with open(filepath, "rb") as f:
        content = f.read()
    try:
        records, fmt = string_to_record(content.decode("utf-8"))
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
            record.id = id.replace(" ", "_")[:max_name_length]
        record.id = id
    elif id is not None:
        record.id = id
        record.id = id.replace(" ", "_")[:max_name_length]

    return record


def load_records_from_files(
    files=None, folder=None, use_file_names_as_ids=False
):
    """Automatically convert files or a folder's content to biopython records.
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
            records += load_records_from_zip_file(filepath)
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
            record.seq.alphabet = DNAAlphabet()
            # Sorry for this parts, it took a lot of "whatever works".
            # keep your part names under 20c and pointless, and everything
            # will be good
            if str(record.id).strip() in UNKNOWN_IDS:
                record.id = name
            if str(record.id).strip() in UNKNOWN_IDS:
                record.id = name
            record.file_name = name_no_extension
        records += recs
    if use_file_names_as_ids:
        for record in records:
            basename = os.path.basename(record.source_file)
            basename_no_extension = os.path.splitext(basename)[0]
            record.id = basename_no_extension
    return records


def write_record(record, target, fmt="genbank"):
    """Write a record as genbank, fasta, etc. via Biopython, with fixes"""
    record = deepcopy(record)
    record.id = record.id[:20]
    if str(record.seq.alphabet.__class__.__name__) != "DNAAlphabet":
        record.seq.alphabet = DNAAlphabet()
    if hasattr(target, "open"):
        target = target.open("w")
    SeqIO.write(record, target, fmt)
