from __future__ import annotations

import pytest
from ase.build import molecule
from ase.calculators.emt import EMT
from ase.units import invcm

from quacc.runners.thermo import ThermoRunner


def test_run_ideal_gas():
    co2 = molecule("CO2")
    igt = ThermoRunner(co2, [526, 526, 1480, 2565]).run_ideal_gas(spin_multiplicity=2)
    assert igt.geometry == "linear"
    assert igt.spin == 0.5

    co2 = molecule("CO2")
    co2.calc = EMT()
    co2.calc.results["magmom"] = 1.0
    igt = ThermoRunner(co2, [526, 526, 1480, 2565]).run_ideal_gas()
    assert igt.spin == 0.5
    assert igt.get_ZPE_correction() == pytest.approx(2548.5 * invcm)

    co2 = molecule("CO2")
    co2.calc = EMT()
    co2.calc.results["magmom"] = 1.0
    igt = ThermoRunner(co2, [-12, 526, 526, 1480, 2565]).run_ideal_gas()
    assert igt.get_ZPE_correction() == pytest.approx(2548.5 * invcm)


def test_run_harmonic_thermo():
    co2 = molecule("CO2")
    ht = ThermoRunner(co2, [526, 526, 1480, 2565]).run_harmonic_thermo()
    assert ht.get_ZPE_correction() == pytest.approx(2548.5 * invcm)

    co2 = molecule("CO2")
    ht = ThermoRunner(co2, [-12, 526, 526, 1480, 2565]).run_harmonic_thermo()
    assert ht.get_ZPE_correction() == pytest.approx(2548.5 * invcm)