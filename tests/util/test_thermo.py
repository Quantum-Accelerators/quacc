from ase.build import molecule
from ase.calculators.emt import EMT

from quacc.util.thermo import ideal_gas


def test_ideal_gas():
    co2 = molecule("CO2")
    igt = ideal_gas(co2, [526, 526, 1480, 2565], spin_multiplicity=2)
    assert igt.geometry == "linear"
    assert igt.spin == 0.5

    co2 = molecule("CO2")
    co2.calc = EMT()
    co2.calc.results["magmom"] = 1.0
    igt = ideal_gas(co2, [526, 526, 1480, 2565])
    assert igt.spin == 0.5
