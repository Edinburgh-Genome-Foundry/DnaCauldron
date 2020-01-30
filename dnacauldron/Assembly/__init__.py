from .builtin_assembly_classes import (
    Type2sRestrictionAssembly,
    BioBrickStandardAssembly,
    BASICAssembly,
    GibsonAssembly,
    OligoPairAnnealing
)
from .AssemblyReportWriter import AssemblyPlotTranslator, AssemblyReportWriter
from .AssemblySimulation import AssemblySimulation
ASSEMBLY_CLASS_DICT = {
    "type2s_assembly": Type2sRestrictionAssembly,
    "gibson_assembly": GibsonAssembly,
    "BASIC_assembly": BASICAssembly,
    "biobrick_assembly": BioBrickStandardAssembly,
    "oligo_annealing": OligoPairAnnealing
}

__all__ = [
    "Type2sRestrictionAssembly",
    "BioBrickStandardAssembly",
    "BASICAssembly",
    "GibsonAssembly",
    "AssemblyPlotTranslator",
    "AssemblyReportWriter",
    "OligoPairAnnealing",
    "AssemblySimulation",
    "ASSEMBLY_CLASS_DICT"
]
