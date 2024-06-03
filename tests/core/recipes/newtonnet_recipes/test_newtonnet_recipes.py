from __future__ import annotations

import pytest

pytest.importorskip("sella")
pytest.importorskip("newtonnet")
from pathlib import Path

import numpy as np
from ase.build import molecule

from quacc import SETTINGS
from quacc.recipes.newtonnet.core import freq_job, relax_job, static_job
from quacc.recipes.newtonnet.ts import irc_job, quasi_irc_job, ts_job

DEFAULT_SETTINGS = SETTINGS.model_copy()


def setup_module():
    current_file_path = Path(__file__).parent

    SETTINGS.NEWTONNET_CONFIG_PATH = current_file_path / "config0.yml"
    SETTINGS.NEWTONNET_MODEL_PATH = current_file_path / "best_model_state.tar"
    SETTINGS.CHECK_CONVERGENCE = False


def teardown_module():
    SETTINGS.NEWTONNET_CONFIG_PATH = DEFAULT_SETTINGS.NEWTONNET_CONFIG_PATH
    SETTINGS.NEWTONNET_MODEL_PATH = DEFAULT_SETTINGS.NEWTONNET_MODEL_PATH
    SETTINGS.CHECK_CONVERGENCE = DEFAULT_SETTINGS.CHECK_CONVERGENCE


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O")
    output = static_job(atoms)
    assert output["spin_multiplicity"] == 1
    assert output["natoms"] == len(atoms)
    assert output["results"]["energy"] == pytest.approx(-9.515200426406743)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())


def test_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O")
    output = relax_job(atoms)
    assert output["spin_multiplicity"] == 1
    assert output["natoms"] == len(atoms)
    assert output["results"]["energy"] == pytest.approx(-9.517354818364769)
    assert not np.array_equal(output["atoms"].get_positions(), atoms.get_positions())
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01


def test_freq_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = molecule("H2O")
    output = freq_job(atoms)
    assert output["atoms"] == molecule("H2O")
    assert len(output["results"]["vib_freqs_raw"]) == 9
    assert len(output["results"]["vib_freqs"]) == 3
    assert output["results"]["vib_freqs_raw"][-1] == pytest.approx(4090.37777396351)
    assert output["results"]["vib_freqs"][0] == pytest.approx(1814.0941260498644)
    assert output["results"]["vib_freqs"][-1] == pytest.approx(4090.37777396351)
    assert output["results"]["n_imag"] == 0
    assert output["results"]["imag_vib_freqs"] == []

    assert output["symmetry"]["point_group"] == "C2v"
    assert output["symmetry"]["rotation_number"] == 2
    assert output["symmetry"]["linear"] is False
    assert output["results"]["energy"] == pytest.approx(-9.515200426406743)
    assert output["results"]["enthalpy"] == pytest.approx(-8.807932688921495)
    assert output["results"]["entropy"] == pytest.approx(0.0019582788098945945)
    assert output["results"]["gibbs_energy"] == pytest.approx(-9.391793516091568)

    atoms = molecule("CH3")
    output = freq_job(atoms, temperature=1000, pressure=20)
    assert output["atoms"] == molecule("CH3")
    assert len(output["results"]["vib_freqs_raw"]) == 12
    assert len(output["results"]["vib_freqs"]) == 6
    assert output["results"]["vib_energies_raw"][0] == pytest.approx(
        -0.09441402482739979
    )
    assert output["results"]["vib_energies_raw"][-1] == pytest.approx(
        0.3925829460532815
    )
    assert output["results"]["vib_energies"][0] == pytest.approx(-0.09441402482739979)
    assert output["results"]["vib_energies"][-1] == pytest.approx(0.3925829460532815)
    assert output["results"]["n_imag"] == 1
    assert output["results"]["imag_vib_freqs"] == pytest.approx([-761.5004719152678])


def test_ts_job_with_default_args(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Define test inputs
    atoms = molecule("H2O")

    # Call the function
    output = ts_job(atoms)

    # Perform assertions on the result
    assert isinstance(output, dict)

    assert "freq_job" in output
    assert output["results"]["energy"] == pytest.approx(-6.796914263061945)
    assert output["freq_job"]["results"]["imag_vib_freqs"][0] == pytest.approx(
        -2426.7398321816004
    )


def test_ts_job_with_custom_hessian(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Define test inputs
    atoms = molecule("H2O")
    use_custom_hessian = True
    opt_params = {"max_steps": 4}
    # Call the function
    output = ts_job(atoms, use_custom_hessian=use_custom_hessian, opt_params=opt_params)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "results" in output["freq_job"]


def test_irc_job_with_default_args(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Define test inputs
    atoms = molecule("H2O")

    # Call the function
    output = irc_job(atoms)

    # Perform assertions on the result
    assert isinstance(output, dict)

    assert output["results"]["energy"] == pytest.approx(-9.517354091813969)
    assert output["freq_job"]["results"]["energy"] == pytest.approx(-9.517354091813969)


def test_irc_job_with_custom_fmax(tmp_path, monkeypatch):
    from ase.build import molecule

    from quacc.recipes.newtonnet.ts import irc_job

    monkeypatch.chdir(tmp_path)
    # Define test inputs
    atoms = molecule("H2O")
    fmax = 0.001

    # Call the function
    output = irc_job(atoms, opt_params={"fmax": fmax})

    # Perform assertions on the result
    assert isinstance(output, dict)

    assert output["results"]["energy"] == pytest.approx(-9.517354091813969)
    assert output["freq_job"]["results"]["energy"] == pytest.approx(-9.517354091813969)


def test_irc_job_with_custom_max_steps(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Define test inputs
    atoms = molecule("H2O")
    max_steps = 500

    # Call the function
    output = irc_job(atoms, opt_params={"max_steps": max_steps})

    # Perform assertions on the result
    assert isinstance(output, dict)

    assert output["results"]["energy"] == pytest.approx(-9.517354091813969)
    assert output["freq_job"]["results"]["energy"] == pytest.approx(-9.517354091813969)


def test_irc_job_with_custom_temperature_and_pressure(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Define test inputs
    atoms = molecule("H2O")
    temperature = 500.0
    pressure = 10.0

    # Call the function
    output = irc_job(
        atoms, freq_job_kwargs={"temperature": temperature, "pressure": pressure}
    )

    # Perform assertions on the result
    assert isinstance(output, dict)

    assert output["results"]["energy"] == pytest.approx(-9.517354091813969)
    assert output["freq_job"]["results"]["energy"] == pytest.approx(-9.517354091813969)


def test_irc_job_with_custom_opt_params(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Define test inputs
    atoms = molecule("H2O")
    opt_params = {"run_kwargs": {"direction": "reverse"}}

    # Call the function
    output = irc_job(atoms, opt_params=opt_params)

    # Perform assertions on the result
    assert isinstance(output, dict)

    assert output["results"]["energy"] == pytest.approx(-9.517354965639784)
    assert output["freq_job"]["results"]["energy"] == pytest.approx(-9.517354965639784)


def test_quasi_irc_job_with_default_args(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Define test inputs
    atoms = molecule("H2O")

    # Call the function
    output = quasi_irc_job(atoms)

    # Perform assertions on the result
    assert isinstance(output, dict)

    assert output["irc_job"]["results"]["energy"] == pytest.approx(-9.517354091813969)
    assert output["results"]["energy"] == pytest.approx(-9.517354091813969)
    assert output["freq_job"]["results"]["energy"] == pytest.approx(-9.517354091813969)


def test_quasi_irc_job_with_custom_direction(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Define test inputs
    atoms = molecule("H2O")
    direction = "reverse"

    # Call the function
    output = quasi_irc_job(atoms, direction=direction)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc_job" in output

    assert output["irc_job"]["results"]["energy"] == pytest.approx(-9.517354965639784)
    assert output["irc_job"]["results"]["energy"] == pytest.approx(-9.517354965639784)
    assert output["results"]["energy"] == pytest.approx(-9.517354965639784)
    assert output["freq_job"]["results"]["energy"] == pytest.approx(-9.517354965639784)


def test_quasi_irc_job_with_custom_temperature_and_pressure(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Define test inputs
    atoms = molecule("H2O")
    temperature = 500.0
    pressure = 10.0

    # Call the function
    output = quasi_irc_job(
        atoms, freq_job_kwargs={"temperature": temperature, "pressure": pressure}
    )

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc_job" in output

    assert output["irc_job"]["results"]["energy"] == pytest.approx(-9.517354091813969)
    assert output["results"]["energy"] == pytest.approx(-9.517354091813969)
    assert output["freq_job"]["results"]["energy"] == pytest.approx(-9.517354091813969)


def test_quasi_irc_job_with_custom_irc_swaps(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Define test inputs
    atoms = molecule("H2O")
    irc_job_kwargs = {"run_kwargs": {"direction": "reverse"}}

    # Call the function
    output = quasi_irc_job(atoms, irc_job_kwargs=irc_job_kwargs)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc_job" in output

    assert output["irc_job"]["results"]["energy"] == pytest.approx(-9.517354965639784)
    assert output["results"]["energy"] == pytest.approx(-9.517354965639784)
    assert output["freq_job"]["results"]["energy"] == pytest.approx(-9.517354965639784)
