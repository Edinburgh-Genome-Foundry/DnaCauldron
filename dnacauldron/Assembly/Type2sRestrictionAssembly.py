from ..AssemblyMix.Type2sRestrictionMix import Type2sRestrictionMix
from .AssemblyBase import AssemblyBase
from ..AssemblyMix.Filter import FragmentSetContainsPartsFilter

class Type2sRestrictionAssembly(AssemblyBase):

    spreadsheet_import_parameters = (
        "enzyme",
        "level",
        "expected_constructs",
        "connectors_collection",
    )

    spreadsheet_columns = [
        "enzyme",
    ]

    def __init__(
        self,
        name,
        parts,
        enzyme="auto",
        connectors_collection=None,
        expected_constructs=1,
        max_constructs=40,
        no_skipped_parts=True,
        level=None,
    ):
        self.name = name
        self.parts = parts
        self.enzyme = enzyme
        self.enzymes = [enzyme]
        self.connectors_collection = connectors_collection
        self.expected_constructs = expected_constructs
        self.max_constructs = max_constructs
        self.no_skipped_parts = no_skipped_parts
        self.level = level

    def get_extra_construct_data(self):
        return dict(enzyme=self.enzyme)

    def simulate(self, sequence_repository, annotate_parts_homologies=True):
        records = sequence_repository.get_records(self.parts)
        mix = Type2sRestrictionMix(parts=records, enzyme=self.enzyme)
        self.enzyme = mix.enzyme_name  # in case the enzyme needs autoselection
        connectors_records = self.get_connectors_records(sequence_repository)
        if len(connectors_records) != 0:
            mix.autoselect_connectors(connectors_records)
        filters = ()
        if self.no_skipped_parts:
            filters = (FragmentSetContainsPartsFilter(self.parts),)
        generator = mix.compute_circular_assemblies(
            annotate_parts_homologies=annotate_parts_homologies,
            fragments_sets_filters=filters,
        )
        circular_assemblies = sorted(
            [asm for (i, asm) in zip(range(self.max_constructs), generator)],
            key=lambda asm: str(asm.seq),
        )
        return mix, circular_assemblies
