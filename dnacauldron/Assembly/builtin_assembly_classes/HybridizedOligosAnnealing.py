from copy import deepcopy
from ..AssemblyFlaw import AssemblyFlaw
from ..AssemblySimulation import AssemblySimulation
from ..Assembly import Assembly
from ...Fragment.HomologousFragment import HomologyChecker
from ...Fragment.StickyEndFragment import (
    StickyEnd,
    StickyEndSeq,
    StickyEndFragment,
)


class HybridizedOligosAnnealing(Assembly):
    def __init__(
        self,
        parts,
        name="annealing",
        homology_checker="default",
        max_constructs=1,
        dependencies=None,
    ):
        Assembly.__init__(
            self,
            name=name,
            parts=parts,
            max_constructs=max_constructs,
            dependencies=dependencies,
        )
        if homology_checker == "default":
            homology_checker = HomologyChecker()
        self.homology_checker = homology_checker

    def simulate(self, sequence_repository):
        p1, p2 = sequence_repository.get_records(self.parts)

        def internal_annealing(p1_str, p2rv):
            i = p1_str.index(p2rv)
            start, end = i, i + len(p2rv)
            seq = StickyEndSeq(
                p1_str[start:end],
                left_end=StickyEnd(p1_str[:start], strand=1),
                right_end=StickyEnd(p1_str[end:], strand=1),
            )
            result = StickyEndFragment(seq)
            result.id = result.original_part = self.name
            p1_copy = deepcopy(p1)
            p2_copy = deepcopy(p2)
            p1_copy.is_reversed = False
            p2_copy.is_reversed = True
            result.fragments = [p1_copy, p2_copy]
            result.annotations['topology'] = 'linear'
            return result

        products = []
        p1_fwd = str(p1.seq)
        p2_fwd = str(p2.seq)
        p2_rv = str(p2.seq.reverse_complement())
        p1_rv = str(p1.seq.reverse_complement())
        
        
        if p2_rv in p1_fwd:
            products.append(internal_annealing(p1_fwd, p2_rv))
        if p1_rv in p2_fwd:
            products.append(internal_annealing(p2_fwd, p1_rv))

        errors = []
        if len(products) == 0:
            msg = "No homology found between oligos %s and %s" % (p1.id, p2.id)
            errors.append(AssemblyFlaw(message=msg, assembly=self))

        return AssemblySimulation(
            assembly=self,
            sequence_repository=sequence_repository,
            construct_records=products,
            errors=errors
        )
