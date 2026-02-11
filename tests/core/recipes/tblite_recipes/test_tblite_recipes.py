from __future__ import annotations

import pytest

pytest.importorskip("tblite.ase")
from copy import deepcopy

import numpy as np
from ase.build import bulk, molecule
from numpy.testing import assert_array_equal

from quacc import change_settings
from quacc.recipes.tblite.core import freq_job, relax_job, static_job


def test_static_job_v1(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O")
    output = static_job(atoms)
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert output["parameters"]["method"] == "GFN2-xTB"
    assert output["results"]["energy"] == pytest.approx(-137.96777594361672)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())
    assert_array_equal(output["atoms"].get_initial_magnetic_moments(), [0.0, 0.0, 0.0])


def test_static_job_v2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O")
    output = static_job(atoms, method="GFN1-xTB")
    assert output["parameters"]["method"] == "GFN1-xTB"
    assert output["results"]["energy"] == pytest.approx(-156.96750578831137)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())


def test_static_job_v3(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H")
    atoms.set_initial_magnetic_moments([1.0])
    output = static_job(atoms, method="GFN1-xTB")
    assert output["parameters"]["method"] == "GFN1-xTB"
    assert output["results"]["energy"] == pytest.approx(-10.92345239113973)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())
    assert_array_equal(output["atoms"].get_initial_magnetic_moments(), [1.0])


def test_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O")
    output = relax_job(atoms)
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert output["parameters"]["method"] == "GFN2-xTB"
    assert output["results"]["energy"] == pytest.approx(-137.97654191396492)
    assert not np.array_equal(output["atoms"].get_positions(), atoms.get_positions())
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01


def test_relax_job_cell(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    atoms.set_initial_magnetic_moments([1.0])
    output = relax_job(atoms, method="GFN1-xTB", relax_cell=True)
    assert output["parameters"]["method"] == "GFN1-xTB"
    assert output["trajectory_results"][-1]["energy"] == pytest.approx(
        -130.4917433263215
    )
    assert "spin_multiplicity" not in output


def test_freq_job_v1(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O")
    output = freq_job(atoms)
    assert output["atoms"] == molecule("H2O")
    assert len(output["results"]["vib_freqs_raw"]) == 9
    assert len(output["results"]["vib_freqs"]) == 3
    assert output["results"]["vib_freqs_raw"][0] == pytest.approx(-0.10864429415434408)
    assert output["results"]["vib_freqs_raw"][-1] == pytest.approx(3526.9940431752034)
    assert output["results"]["vib_freqs"][0] == pytest.approx(1586.623114694335)
    assert output["results"]["vib_freqs"][-1] == pytest.approx(3526.9940431752034)
    assert output["results"]["n_imag"] == 0
    assert output["results"]["imag_vib_freqs"] == []

    assert output["molecule_metadata"]["symmetry"]["point_group"] == "C2v"
    assert output["molecule_metadata"]["symmetry"]["rotation_number"] == 2
    assert output["molecule_metadata"]["symmetry"]["linear"] is False
    assert len(output["parameters_thermo"]["vib_freqs"]) == 3
    assert output["results"]["vib_freqs"][0] == pytest.approx(1586.623114694335)
    assert output["parameters_thermo"]["vib_freqs"][-1] == pytest.approx(
        3526.9940431752034
    )
    assert output["results"]["energy"] == 0.0
    assert output["results"]["enthalpy"] == pytest.approx(0.6375973622705722)
    assert output["results"]["entropy"] == pytest.approx(0.0019584992229988523)
    assert output["results"]["gibbs_energy"] == pytest.approx(0.05367081893346437)


def test_freq_job_v2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H")
    atoms.set_initial_magnetic_moments([1.0])
    initial_atoms = deepcopy(atoms)
    output = freq_job(atoms, energy=-1.0)
    assert output["atoms"] == initial_atoms
    assert len(output["results"]["vib_freqs_raw"]) == 3
    assert len(output["results"]["vib_freqs"]) == 0
    assert output["results"]["vib_freqs_raw"][0] == 0
    assert output["results"]["vib_freqs_raw"][-1] == 0
    assert output["results"]["vib_freqs"] == []
    assert output["results"]["n_imag"] == 0
    assert output["results"]["imag_vib_freqs"] == []

    assert output["molecule_metadata"]["symmetry"]["linear"] is False
    assert len(output["parameters_thermo"]["vib_freqs"]) == 0
    assert output["results"]["energy"] == -1.0
    assert output["results"]["enthalpy"] == pytest.approx(-0.9357685739989672)
    assert output["results"]["entropy"] == pytest.approx(0.0011889660575166018)
    assert output["results"]["gibbs_energy"] == pytest.approx(-1.290258804047542)


def test_freq_job_v3(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("CH3")
    atoms.set_initial_magnetic_moments([1.0, 0.0, 0.0, 0.0])
    initial_atoms = deepcopy(atoms)
    output = freq_job(atoms, energy=-10.0, temperature=1000, pressure=20)
    assert output["atoms"] == initial_atoms
    assert len(output["results"]["vib_freqs_raw"]) == 12
    assert len(output["results"]["vib_freqs"]) == 6
    assert output["results"]["vib_energies_raw"][0] == pytest.approx(
        -9.551076713062095e-06
    )
    assert output["results"]["vib_energies_raw"][-1] == pytest.approx(
        0.3880868821616259
    )
    assert output["results"]["vib_energies"][0] == pytest.approx(0.0713506770137291)
    assert output["results"]["vib_energies"][-1] == pytest.approx(0.3880868821616259)
    assert output["results"]["n_imag"] == 0
    assert output["results"]["imag_vib_freqs"] == []

    assert output["parameters_thermo"]["temperature"] == 1000.0
    assert output["parameters_thermo"]["pressure"] == 20.0
    assert output["parameters_thermo"]["sigma"] == 6
    assert output["parameters_thermo"]["spin_multiplicity"] == 2
    assert output["molecule_metadata"]["symmetry"]["linear"] is False
    assert output["molecule_metadata"]["symmetry"]["rotation_number"] == 6
    assert len(output["parameters_thermo"]["vib_freqs"]) == 6
    assert output["results"]["energy"] == -10.0
    assert output["results"]["enthalpy"] == pytest.approx(-8.749341973959462)
    assert output["results"]["entropy"] == pytest.approx(0.0023506788982171896)
    assert output["results"]["gibbs_energy"] == pytest.approx(-11.100020872176652)
    assert "nid" in output
    assert "dir_name" in output
    assert_array_equal(
        output["atoms"].get_initial_magnetic_moments(), [1.0, 0.0, 0.0, 0.0]
    )


def test_unique_workdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    with change_settings({"CREATE_UNIQUE_DIR": True}):
        test_static_job_v1(tmp_path, monkeypatch)
        test_relax_job(tmp_path, monkeypatch)
