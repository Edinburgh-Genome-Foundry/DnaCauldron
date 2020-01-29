# from copy import copy
# from ...biotools import set_record_topology
# from ...biotools import annotate_record
from ..Assembly import Assembly
from ...AssemblyMix import RestrictionLigationMix
from ..AssemblySimulation import AssemblySimulation
from ..AssemblyFlaw import AssemblyFlaw


class BioBrickStandardAssembly(Assembly):
    """Representation and simulation of the Biobrick 2-part assembly standard.

    Parameters
    ----------

    parts
      A list of parts names corresponding to records in a repository.
      There must be exactly 2 parts and they must be represented on a
      backbone (i.e. circular constructs), and the first part will be inserted
      in the backbone of the second part, upstream of the second part.

    name
      Name of the assembly as it will appear in reports.

    max_constructs
      None or a number of maximum assemblies to compute (avoids complete
      freeze for combinatorial assemblies with extremely many possibilities).

    expected_constructs
      Either a number or a string ``'any_number'``. If the number of constructs
      doesn't match this value, the assembly will be considered invalid in
      reports and summaries

    connectors_collection
      Name of a collection in the repository from which to get candidates for
      connector autocompletion.

    dependencies
      (do not use). Metadata indicating which assemblies depend on this
      assembly, or are depended on by it.
    """
    def __init__(
        self,
        parts,
        name="unnamed_assembly",
        connectors_collection=None,
        expected_constructs=1,
        max_constructs=40,
        dependencies=None,
    ):
        Assembly.__init__(
            self,
            name=name,
            parts=parts,
            max_constructs=max_constructs,
            dependencies=dependencies,
        )
        self.connectors_collection = connectors_collection
        self.expected_constructs = expected_constructs
        self.enzymes = ["EcoRI", "SpeI", "XbaI"]
        self.extra_construct_data = dict(enzymes=self.enzymes)

    def simulate(self, sequence_repository, annotate_parts_homologies=True):
        """Simulate the assembly, return an AssemblySimulation."""
        if len(self.parts) != 2:
            error = AssemblyFlaw(
                assembly=self,
                message="BioBrick Std Assembly requires 2 parts exactly.",
                suggestion="Check assembly plan",
                data={"parts": len(self.parts)},
            )
            AssemblySimulation(
                assembly=self,
                sequence_repository=sequence_repository,
                errors=[error],
            )

        left_part, right_part = sequence_repository.get_records(self.parts)
        errors, mixes = [], []

        # E+S MIX TO DIGEST THE LEFT PART (INSERT)

        mix_1 = RestrictionLigationMix(
            parts=[left_part],
            enzymes=["EcoRI", "SpeI"],
            name=left_part.id + "_E+S_restriction",
        )
        mixes.append(mix_1)
        inserts = [
            fragment
            for fragment in mix_1.fragments + mix_1.reverse_fragments
            if fragment.seq.ends_tuple() == ("AATT", "CTAG")
            and not fragment.is_reversed
        ]
        if len(inserts) != 1:
            prefix = "No" if len(inserts) == 0 else "Multiple"
            error = AssemblyFlaw(
                assembly=self,
                message=prefix + " AATT-CTAG fragments in E+S digestion.",
                suggestion=(
                    "Check EcoRI/SpeI sites in %s. Also beware that the E+S "
                    "part should appear in direct sense in its record."
                )
                % left_part.id,
                data={"parts": len(self.parts)},
            )
            errors.append(error)

        # E+X MIX TO DIGEST THE RIGHT PART (BACKBONE)

        mix_2 = RestrictionLigationMix(
            parts=[right_part],
            enzymes=["EcoRI", "XbaI"],
            name=right_part.id + "_E+X_restriction",
        )
        mixes.append(mix_2)
        backbones = [
            fragment
            for fragment in mix_2.fragments + mix_2.reverse_fragments
            if fragment.seq.ends_tuple() == ("CTAG", "AATT")
            and len(fragment) > 100
        ]
        if len(backbones) != 1:
            prefix = "No" if len(backbones) == 0 else "Multiple"
            error = AssemblyFlaw(
                assembly=self,
                message=prefix + " left-CTAG overhang found in E+X mix.",
                suggestion="Check EcoRI/XbaI sites in %s" % right_part.id,
                data={"parts": len(self.parts)},
            )
            errors.append(error)

        if len(errors):
            return AssemblySimulation(
                assembly=self,
                errors=errors,
                mixes=mixes,
                sequence_repository=sequence_repository,
            )

        # FINAL ASSEMBLY MIX

        insert = inserts[0]
        backbone = backbones[0]

        mix = RestrictionLigationMix(
            fragments=[insert, backbone], name="assembly_mix"
        )
        mixes.append(mix)
        final_record = insert.assemble_with(backbone).circularized()
        final_record.id = self.name
        final_record.fragments = [insert, backbone]
        return AssemblySimulation(
            assembly=self,
            construct_records=[final_record],
            mixes=mixes,
            errors=errors,
            sequence_repository=sequence_repository,
        )

        # def fragments_set_filter(fragments):
        #     if len(fragments) != 2:
        #         return False
        #     f1, f2 = fragments
        #     uses_both_parts = f1.original_part.id != f2.original_part.id
        #     return uses_both_parts and not f1.is_reversed

        # generator = mix.compute_circular_assemblies(
        #     annotate_parts_homologies=annotate_parts_homologies,
        #     fragments_sets_filters=(fragments_set_filter,),
        # )
        # construct_records = sorted(
        #     [asm for (i, asm) in zip(range(self.max_constructs), generator)],
        #     key=lambda asm: str(asm.seq),
        # )
        # self.attribute_ids_to_constructs(construct_records)
        # n_constructs = len(construct_records)
        # if n_constructs != 1:
        #     error = AssemblyFlaw(
        #         assembly=self,
        #         message="Found %d constructs instead of 1." % n_constructs,
        #         suggestion="Hmmmm",
        #     )
        #     errors.append(error)

        # return AssemblySimulation(
        #     assembly=self,
        #     construct_records=construct_records,
        #     mixes=mixes,
        #     errors=errors,
        #     sequence_repository=sequence_repository,
        # )
