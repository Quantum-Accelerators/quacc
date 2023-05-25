from ase.build import molecule
from ase.calculators.emt import EMT

from quacc.util.thermo import ideal_gas


def test_ideal_gas():
    co2 = molecule("CO2")
    igt = ideal_gas(co2, spin_multiplicity=3)

    co2 = molecule("CO2")
    co2.calc = EMT()
    co2.calc.results["magmom"] = 1.0
    igt = ideal_gas(co2)
