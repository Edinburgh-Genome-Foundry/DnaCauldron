from .AssemblyBase import AssemblyBase
from ..AssemblyMix import RestrictionLigationMix


class BioBrickStandardAssembly(AssemblyBase):
    def __init__(
        self,
        name,
        parts,
        connectors_collection=None,
        expected_constructs=1,
        max_constructs=40,
        level=None,
    ):
        self.name = name
        self.parts = parts
        self.connectors_collection = connectors_collection
        self.expected_constructs = expected_constructs
        self.max_constructs = max_constructs
        self.level = level
        self.enzymes = ["EcoRI", "SpeI", "XbaI"]
        self.extra_construct_data = dict(enzymes=self.enzymes)

    def simulate(self, sequences_repository, annotate_parts_homologies=True):
        left_part, right_part = sequences_repository.get_records(self.parts)
        E, X, S = "EcoRI", "XbaI", "SpeI"
        mix_1 = RestrictionLigationMix(parts=[left_part], enzymes=[E, S])
        fragments = mix_1.fragments + mix_1.reverse_fragments
        insert = [f for f in fragments if str(f.seq.right_end) == "CTAG"][0]

        mix_2 = RestrictionLigationMix(parts=[right_part], enzymes=[E, X])
        fragments = mix_2.fragments + mix_2.reverse_fragments
        backbone = [f for f in fragments if str(f.seq.left_end) == "CTAG"][0]
        mix = RestrictionLigationMix(fragments=[insert, backbone])

        def fragments_set_filter(fragments):
            if len(fragments) != 2:
                return False
            f1, f2 = fragments
            uses_both_parts = (
                f1.original_part.id != f2.original_part.id
            )
            return uses_both_parts and not f1.is_reverse

        generator = mix.compute_circular_assemblies(
            annotate_parts_homologies=annotate_parts_homologies,
            fragments_sets_filters=(fragments_set_filter,),
        )
        circular_assemblies = sorted(
            [asm for (i, asm) in zip(range(self.max_constructs), generator)],
            key=lambda asm: str(asm.seq),
        )
        return mix, circular_assemblies
