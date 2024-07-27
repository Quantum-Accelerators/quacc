from __future__ import annotations

import logging

import numpy as np
import pytest
from ase.build import molecule
from ase.calculators.emt import EMT
from ase.units import invcm
from monty.json import MontyDecoder, jsanitize

from quacc.schemas.thermo import ThermoSummarize

LOGGER = logging.getLogger(__name__)
LOGGER.propagate = True


def test_run_ideal_gas(tmp_path):
    co2 = molecule("CO2")
    igt = ThermoSummarize(
        co2, [526, 526, 1480, 2565], directory=tmp_path
    )._make_ideal_gas()
    assert igt.geometry == "linear"
    assert igt.spin == 0
    assert igt.get_ZPE_correction() == pytest.approx(2548.5 * invcm)

    co2 = molecule("CO2")
    igt = ThermoSummarize(
        co2, [526, 526, 1480, 2565], directory=tmp_path
    )._make_ideal_gas(spin_multiplicity=2)
    assert igt.geometry == "linear"
    assert igt.spin == 0.5
    assert igt.get_ZPE_correction() == pytest.approx(2548.5 * invcm)

    co2 = molecule("CO2")
    igt = ThermoSummarize(
        co2, [-12, 526, 526, 1480, 2565], directory=tmp_path
    )._make_ideal_gas(spin_multiplicity=2)
    assert igt.geometry == "linear"
    assert igt.spin == 0.5
    assert igt.get_ZPE_correction() == pytest.approx(2548.5 * invcm)


def test_summarize_ideal_gas_thermo(tmp_path, caplog):
    # Make sure metadata is made
    atoms = molecule("N2")
    results = ThermoSummarize(atoms, [0.34 / invcm], directory=tmp_path).ideal_gas()
    assert results["natoms"] == len(atoms)
    assert results["atoms"] == atoms
    assert results["parameters_thermo"]["vib_energies"] == [0.34]
    assert results["parameters_thermo"]["vib_freqs"] == [0.34 / invcm]
    assert results["results"]["energy"] == 0
    assert "pymatgen_version" in results["builder_meta"]

    # Make sure right number of vib energies are reported
    atoms = molecule("N2")
    results = ThermoSummarize(
        atoms, [0.0 / invcm, 0.34 / invcm], energy=-1, directory=tmp_path
    ).ideal_gas()
    assert results["natoms"] == len(atoms)
    assert results["atoms"] == atoms
    assert results["parameters_thermo"]["vib_energies"] == [0.34]
    assert results["parameters_thermo"]["vib_freqs"] == [0.34 / invcm]
    assert results["results"]["energy"] == -1

    # # # Make sure info tags are handled appropriately
    atoms = molecule("N2")
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    atoms.calc = EMT()
    results = ThermoSummarize(
        atoms, [0.0 / invcm, 0.34 / invcm], energy=-1, directory=tmp_path
    ).ideal_gas()
    assert results["atoms"].info.get("test_dict", None) == {"hi": "there", "foo": "bar"}

    # Make sure spin works right
    atoms = molecule("CH3")
    vib_energies = [
        9.551077150221621e-06j,
        3.1825877476455407e-06j,
        2.7223332245579342e-06j,
        (0.03857802457526743 + 0j),
        (0.038762952842240087 + 0j),
        (0.03876411386029907 + 0j),
        (0.07135067701372912 + 0j),
        (0.1699785717790056 + 0j),
        (0.1700229358789492 + 0j),
        (0.3768719400148424 + 0j),
        (0.38803854931751625 + 0j),
        (0.3880868821616261 + 0j),
    ]
    with caplog.at_level(logging.INFO):
        results = ThermoSummarize(
            atoms, np.array(vib_energies) / invcm, energy=-10.0, directory=tmp_path
        ).ideal_gas(temperature=1000.0, pressure=20.0)
    assert (
        "No multiplicity provided. Automatically detecting a spin multiplicity of 2 from the Atoms object"
        in caplog.text
    )

    assert results["natoms"] == len(atoms)
    assert results["atoms"] == atoms
    assert len(results["parameters_thermo"]["vib_energies"]) == 6
    assert results["parameters_thermo"]["vib_energies"][0] == vib_energies[-6]
    assert results["parameters_thermo"]["vib_energies"][-1] == vib_energies[-1]
    assert results["results"]["energy"] == -10.0
    assert results["results"]["enthalpy"] == pytest.approx(-8.749341973959462)
    assert results["results"]["entropy"] == pytest.approx(0.0023506788982171896)
    assert results["results"]["gibbs_energy"] == pytest.approx(-11.100020872176652)
    assert results["parameters_thermo"]["temperature"] == 1000.0
    assert results["parameters_thermo"]["pressure"] == 20.0
    assert results["parameters_thermo"]["sigma"] == 6
    assert results["parameters_thermo"]["spin_multiplicity"] == 2

    # Test custom spin
    results = ThermoSummarize(
        atoms,
        np.array(vib_energies) / invcm,
        energy=-10.0,
        charge_and_multiplicity=(0, 2),
        directory=tmp_path,
    ).ideal_gas(temperature=1000.0, pressure=20.0)
    assert results["results"]["entropy"] == pytest.approx(0.0023506788982171896)
    assert results["parameters_thermo"]["spin_multiplicity"] == 2

    # Test custom spin
    results = ThermoSummarize(
        atoms,
        np.array(vib_energies) / invcm,
        energy=-10.0,
        charge_and_multiplicity=(0, 4),
        directory=tmp_path,
    ).ideal_gas(temperature=1000.0, pressure=20.0)
    assert results["results"]["entropy"] == pytest.approx(0.0024104096804891486)
    assert results["parameters_thermo"]["spin_multiplicity"] == 4

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)


def test_run_harmonic_thermo(tmp_path):
    co2 = molecule("CO2")
    ht = ThermoSummarize(
        co2, [526, 526, 1480, 2565], directory=tmp_path
    )._make_harmonic_thermo()
    assert ht.get_ZPE_correction() == pytest.approx(2548.5 * invcm)

    co2 = molecule("CO2")
    ht = ThermoSummarize(
        co2, [-12, 526, 526, 1480, 2565], directory=tmp_path
    )._make_harmonic_thermo()
    assert ht.get_ZPE_correction() == pytest.approx(2548.5 * invcm)


def test_summarize_harmonic_thermo(tmp_path):
    atoms = molecule("H2")

    # Make sure metadata is made
    results = ThermoSummarize(atoms, [0.34 / invcm], directory=tmp_path).harmonic()
    assert results["parameters_thermo"]["vib_energies"] == [0.34]
    assert results["parameters_thermo"]["vib_freqs"] == [0.34 / invcm]
    assert results["results"]["energy"] == 0
    assert results["results"]["helmholtz_energy"] == pytest.approx(0.16999995401497991)
    assert results["results"]["internal_energy"] == pytest.approx(0.1700006085385999)
    assert results["results"]["entropy"] == pytest.approx(2.1952829783392438e-09)
    assert results["results"]["zpe"] == pytest.approx(0.17)

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)
