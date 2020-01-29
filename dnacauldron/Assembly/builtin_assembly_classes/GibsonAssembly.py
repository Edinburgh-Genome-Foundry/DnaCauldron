from ...AssemblyMix import HomologousAssemblyMix, AssemblyMixError
from ...Fragment.HomologousFragment import HomologyChecker
from ..AssemblySimulation import AssemblySimulation
from ..AssemblyFlaw import AssemblyFlaw
from ..Assembly import Assembly

class GibsonAssembly(Assembly):
    """Representation and simulation of Gibson Assembly

    Parameters
    ----------

    parts
      A list of parts names corresponding to records in a repository.
      The parts will be considered as assembling together if they
      have end homologies, as checked by the homology_checker.
    
    homology_checker
      An HomologyChecker instance defining which homology sizes and melting
      temperatures are valid.

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
    
    expect_no_unused_parts
      If True and some parts are unused, this will be considered an invalid
      assembly in summaries and reports.

    dependencies
      (do not use). Metadata indicating which assemblies depend on this
      assembly, or are depended on by it.
    """

    spreadsheet_import_parameters = (
        "level",
        "expected_constructs",
        "expect_no_unused_parts",
        "connectors_collection",
    )

    def __init__(
        self,
        parts,
        homology_checker="default",
        name="homologous_assembly",
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
        self.connectors_collection = connectors_collection
        self.expected_constructs = expected_constructs
        self.expect_no_unused_parts = expect_no_unused_parts
        if homology_checker == "default":
            homology_checker = HomologyChecker()
        self.homology_checker = homology_checker

    def get_extra_construct_data(self):
        description = self.homology_checker.parameters_as_string()
        return dict(homology_condition=description)

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
            return AssemblyFlaw(
                assembly=self,
                message="Some parts are unused",
                suggestion=suggestion,
                data={"unused_parts": sorted(unused_parts)},
            )

    def _detect_parts_connections_errors(self, assembly_mix):
        warnings = []
        graph = assembly_mix.connections_graph
        fragment_degrees = {
            fragment: graph.degree(fragment) for fragment in graph.nodes()
        }

        # CHECK PARTS WITH A SINGLE OVERHANG

        deadend_parts = [
            assembly_mix.fragments_dict[fragment].original_part.id
            for fragment, degree in fragment_degrees.items()
            if degree == 1
        ]
        deadend_parts = sorted(set(deadend_parts))
        if len(deadend_parts):
            warning = AssemblyFlaw(
                assembly=self,
                message="Part(s) with single-sided sticky end",
                suggestion="Check restriction sites or assembly plan",
                data={"parts": deadend_parts},
            )
            warnings.append(warning)

        # CHECK PARTS WITH MORE THAN 2 SLOT CONNECTIONS (FORK/CROSSROAD)

        fork_parts = [
            assembly_mix.fragments_dict[fragment].original_part.id
            for fragment, degree in fragment_degrees.items()
            if degree > 2
        ]
        fork_parts = sorted(set(fork_parts))
        if len(fork_parts):
            warning = AssemblyFlaw(
                assembly=self,
                message="Warning: parts at graph forking positions",
                suggestion="Check restriction sites in part sequences",
                data={"parts": fork_parts},
            )
            warnings.append(warning)
        return warnings

    def simulate(self, sequence_repository, annotate_parts_homologies=True):
        """Simulate the Gibson Assembly, return an AssemblySimulation."""

        # CREATE THE MIX

        warnings = []

        records = sequence_repository.get_records(self.parts)
        mix = HomologousAssemblyMix(
            parts=records,
            homology_checker=self.homology_checker,
            name=self.name + "_homology_mix",
        )

        # ATTEMPT CONNECTOR AUTOSELECTION IF NECESSARY

        connectors_records = self._get_connectors_records(sequence_repository)
        if len(connectors_records) != 0:
            try:
                connectors = mix.autoselect_connectors(connectors_records)
                if len(connectors):
                    connectors_ids = [c.id for c in connectors]
                    connectors_string = ", ".join(connectors_ids)
                    warning = AssemblyFlaw(
                        assembly=self,
                        message="Added connectors %s" % connectors_string,
                        suggestion="",
                        data={"selected_connectors": connectors_ids},
                    )
                    warnings.append(warning)
            except AssemblyMixError as err:
                error = AssemblyFlaw(
                    assembly=self,
                    message="Failed to find suitable connectors",
                    suggestion="Check assembly plan or parts design",
                )
                return AssemblySimulation(
                    assembly=self,
                    construct_records=[],
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
        self._detect_constructs_number_error(found, errors)
        self._detect_max_constructs_reached(found, warnings)
        if self.expect_no_unused_parts:
            unused_error = self._detect_unused_parts_error(construct_records)
            if unused_error is not None:
                errors.append(unused_error)
        if errors != []:
            warnings.extend(self._detect_parts_connections_errors(mix))

        return AssemblySimulation(
            assembly=self,
            construct_records=construct_records,
            mixes=(mix,),
            errors=errors,
            warnings=warnings,
            sequence_repository=sequence_repository,
        )
