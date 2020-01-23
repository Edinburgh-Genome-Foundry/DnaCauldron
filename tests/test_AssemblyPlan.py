import pytest
import dnacauldron as dc


def test_detect_duplicate_assemblies():
    assemblies = [
        dc.GibsonAssembly(name="bla", parts=[]),
        dc.GibsonAssembly(name="ble", parts=[]),
        dc.GibsonAssembly(name="bli", parts=[]),
        dc.GibsonAssembly(name="ble", parts=[]),
    ]
    with pytest.raises(ValueError):
        dc.AssemblyPlan(assemblies)
