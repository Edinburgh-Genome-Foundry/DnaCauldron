from .biotools import (
    load_records_from_files,
    set_record_topology,
    sequence_to_biopython_record,
)
from fuzzywuzzy import process


class NotInRepositoryError(Exception):
    def __init__(self, parts, repository):
        self.parts = parts
        self.repository = repository

        # CREATE THE MESSAGE AND INITIALIZE THE EXCEPTION:

        suggestions = [
            self.create_part_suggestion_string(part_name)
            for part_name in parts
        ]
        suggestions = ", ".join(suggestions)
        message = "Parts not found in %s: %s" % (repository.name, suggestions)
        super().__init__(message)
    
    def create_part_suggestion_string(self, part_name):
        suggestions = self.repository.suggest_part_names(part_name)
        if len(suggestions) == 0:
            return part_name
        return "%s (did you mean %s ?)" % (part_name, " or ".join(suggestions))


class RepositoryDuplicateError(Exception):
    def __init__(self, parts, repository):
        self.parts = parts
        self.repository = repository
        parts_list = ", ".join(parts)
        if len(parts_list) > 150:
            parts_list = parts_list[:150] + "..."
        parts = "Part ID%s %s" % ("s" if len(parts) > 1 else "", parts_list)
        repo_name = (" in " + repository.name) if repository.name else ""
        message = parts + " duplicated in " + repo_name
        super().__init__(message)


class SequenceRepository:
    """Sequence repositories store and provide sequence records.

    The records are organized into collections, for instance "parts" to host
    parts, "constructs" for records created during assembly plan simulation,
    or any other collection name like "emma_connectors" to store EMMA
    connectors.

    The suggested initialization of a sequence repository is:

    >>> repository = SequenceRepository()
    >>> repository.import_records(files=['part.fa', 'records.zip', etc.])


    
    Parameters
    ----------

    collections
      A dict {'collection_name': {'record_id': record, ...}, ...} giving for
      each collection a dict of Biopython records.
    
    name
      The name of the repository as it may appear in error messages and other
      reports.
    """

    def __init__(self, collections=None, name="repo"):
        self.collections = collections or {}
        self.name = name

    def add_record(self, record, collection="parts"):
        """Add one record to a collection, using its record.id as key.
        
        The collection is created if it doesn't exist.

        The record can also be a pair (id, "ATGTGCC...").
        """
        if isinstance(record, (tuple, list)):
            _id, _sequence = record
            record = sequence_to_biopython_record(_sequence, id=_id)
        if self.contains_record(record.id):
            raise RepositoryDuplicateError([record.id], repository=self)
        if collection not in self.collections:
            self.collections[collection] = {}
        self.collections[collection][record.id] = record

    def add_records(self, records, collection="parts"):
        """Add """

        if len(records) == 0:
            return
        for record in records:
            self.add_record(record, collection=collection)

    def contains_record(self, record_id):
        """Return whether the repo has a record corresponding to the given id
        """
        collections = self.collections.values()
        return any(record_id in collection for collection in collections)

    def get_record(self, record_id):
        """Return the record from the repository from its ID."""
        for collection in self.collections.values():
            if record_id in collection:
                return collection[record_id]
        raise NotInRepositoryError([record_id], self)

    def get_records(self, record_ids):
        """Get a list of records from a list of record IDs."""
        records = []
        not_in_repository = []
        for name in record_ids:
            if self.contains_record(name):
                records.append(self.get_record(name))
            else:
                not_in_repository.append(name)
        if len(not_in_repository):
            raise NotInRepositoryError(not_in_repository, repository=self)
        return records

    def import_records(
        self,
        files=None,
        folder=None,
        collection="parts",
        use_file_names_as_ids=True,
        topology="default_to_linear",
    ):
        """Import records into the repository, from files and zips and folders.

        Parameters
        ----------

        files
          A list of file paths, either Genbank, Fasta, Snapgene (.dna), or zips
          containing any of these formats.
        
        folder
          Path to a folder which can be provided instead of ``files``.
        
        collection
          Name of the collection under which to import the new records.
        
        use_file_names_as_ids
          If True, the file name will be used as ID for any record obtained
          from a single-record file (fasta files with many records will still
          use the internal ID).
        
        topology
          Can be "circular", "linear", "default_to_circular" (will default
          to circular if ``annotations['topology']`` is not already set) or
          "default_to_linear".
        """
        if folder is not None:
            records = load_records_from_files(
                folder=folder, use_file_names_as_ids=use_file_names_as_ids
            )
        elif files is not None:
            records = load_records_from_files(
                files=files, use_file_names_as_ids=use_file_names_as_ids,
            )
        else:
            raise ValueError("Provide either ``files`` or ``folder``")
        for r in records:
            set_record_topology(r, topology)

        self.add_records(records, collection=collection)

    def get_part_names_by_collection(self, format="dict"):
        """Return a dictionnary or a string representing the repo's content.
        
        Format: "dict" or "string"
        """
        result = {
            collection_name: list(parts.keys())
            for collection_name, parts in self.collections.items()
        }
        if format == "dict":
            return result
        else:
            return "\n".join(
                "\n".join([name] + ["- " + part for part in sorted(parts)])
                for name, parts in result.items()
            )

    def get_all_part_names(self):
        """Return the list of all part names"""
        parts = [
            part
            for collection in self.collections.values()
            for part in collection
        ]
        return sorted(parts)
    
    def suggest_part_names(self, query, cutoff=90, limit=3):
        """Suggest part names in the repo close to the given query."""
        search = process.extract(query, self.get_all_part_names())
        return [
            name
            for (name, score) in sorted(search, key=lambda e: -e[1])
            if score >= cutoff
        ][:limit]
