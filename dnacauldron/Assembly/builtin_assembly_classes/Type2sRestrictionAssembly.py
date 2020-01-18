from ...biotools import autoselect_enzyme
from ...Filter import NoRestrictionSiteFilter
from ...AssemblyMix import RestrictionLigationMix, AssemblyMixError
from ..Assembly import Assembly
from ..AssemblySimulation import AssemblySimulation
from ..AssemblySimulationError import AssemblySimulationError


class Type2sRestrictionAssembly(Assembly):

    spreadsheet_import_parameters = (
        "enzyme",
        "level",
        "expected_constructs",
        "expect_no_unused_parts",
        "connectors_collection",
    )

    spreadsheet_columns = [
        "enzyme",
    ]

    def __init__(
        self,
        parts,
        name="type2s_assembly",
        enzyme="auto",
        connectors_collection=None,
        expected_constructs=1,
        expect_no_unused_parts=True,
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
        self.enzyme = enzyme
        self.connectors_collection = connectors_collection
        self.expected_constructs = expected_constructs
        self.expect_no_unused_parts = expect_no_unused_parts
        

    def get_extra_construct_data(self):
        return dict(enzymes=self.enzymes)

    def _detect_constructs_number_error(self, found):
        expected = self.expected_constructs
        if expected != found:
            return AssemblySimulationError(
                assembly=self,
                message="Wrong number of constructs",
                suggestion="Check assembly or parts design",
                data={"expected_": expected, "found": found},
            )

    def _detect_unused_parts_error(self, construct_records):
        used_parts = [
            fragment.original_part.id
            for construct_record in construct_records
            for fragment in construct_record.fragments
        ]
        unused_parts = set(self.parts).difference(set(used_parts))
        if len(unused_parts):
            if len(construct_records):
                suggestion = "Check assembly plan"
            else:
                suggestion = "Check parts design"
            return AssemblySimulationError(
                assembly=self,
                message="Some parts are unused",
                suggestion=suggestion,
                data={"unused_parts": sorted(unused_parts)},
            )

    def _detect_parts_connections_errors(self, assembly_mix):
        errors = []
        slots_graph = assembly_mix.slots_graph(with_overhangs=False)
        slots_degrees = {
            slot: slots_graph.degree(slot) for slot in slots_graph.nodes()
        }
        slots_parts = assembly_mix.compute_slots()

        # CHECK PARTS WITH A SINGLE OVERHANG

        deadend_parts = [
            part
            for slot, degree in slots_degrees.items()
            if degree == 1
            for part in slots_parts[slot]
        ]
        if len(deadend_parts):
            error = AssemblySimulationError(
                assembly=self,
                message="Part(s) with single-sided sticky end",
                suggestion="Check restriction sites or assembly plan",
                data={"parts": deadend_parts},
            )
            errors.append(error)

        # CHECK PARTS WITH MORE THAN 2 SLOT CONNECTIONS (FORK/CROSSROAD)

        fork_parts = [
            part
            for slot, degree in slots_degrees.items()
            if degree > 2
            for part in slots_parts[slot]
        ]
        if len(fork_parts):
            error = AssemblySimulationError(
                assembly=self,
                message="Warning: parts at graph forking positions",
                suggestion="Check restriction sites in part sequences",
                data={"parts": fork_parts},
            )
            errors.append(error)
        return errors

    @property
    def enzymes(self):
        return [self.enzyme]

    def simulate(self, sequence_repository, annotate_parts_homologies=True):

        # CREATE THE MIX

        records = sequence_repository.get_records(self.parts)
        if self.enzyme == "auto":
            self.enzyme = autoselect_enzyme(records)
        mix = RestrictionLigationMix(
            parts=records,
            enzymes=[self.enzyme],
            fragments_filters=[NoRestrictionSiteFilter(str(self.enzyme))],
            name=self.name + "_type2s_mix",
        )

        # ATTEMPT CONNECTOR AUTOSELECTION IF NECESSARY

        connectors_records = self.get_connectors_records(sequence_repository)
        if len(connectors_records) != 0:
            try:
                mix.autoselect_connectors(connectors_records)
            except AssemblyMixError as err:
                error = AssemblySimulationError(
                    assembly=self,
                    message="Failed to find suitable connectors",
                    suggestion="Check assembly plan or parts design",
                    mixes=(err.mix,),
                )
                return AssemblySimulation(
                    assembly=self,
                    construct_records=[],
                    mixes=(mix,),
                    errors=[error],
                    sequence_repository=sequence_repository,
                )

        # COMPUTE ALL CIRCULAR ASSEMBLIES

        generator = mix.compute_circular_assemblies(
            annotate_parts_homologies=annotate_parts_homologies
        )
        construct_records = sorted(
            [asm for (i, asm) in zip(range(self.max_constructs), generator)],
            key=lambda asm: str(asm.seq),
        )
        self.attribute_ids_to_constructs(construct_records)

        # FIND ALL ERRORS

        errors = []
        found = len(construct_records)
        constructs_number_error = self._detect_constructs_number_error(found)
        if constructs_number_error is not None:
            errors.append(constructs_number_error)
        if self.expect_no_unused_parts:
            unused_error = self._detect_unused_parts_error(construct_records)
            if unused_error is not None:
                errors.append(constructs_number_error)
        if errors != []:
            errors.extend(self._detect_parts_connections_errors(mix))

        return AssemblySimulation(
            assembly=self,
            construct_records=construct_records,
            mixes=(mix,),
            errors=errors,
            sequence_repository=sequence_repository,
        )
