from __future__ import annotations

import numpy as np
import pytest
from ase.build import bulk, molecule
from maggma.stores import MemoryStore

from quacc import change_settings
from quacc.recipes.lj.core import freq_job, relax_job, static_job


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    with change_settings({"STORE": MemoryStore()}):
        atoms = molecule("H2O")

        output = static_job(atoms)
        assert output["molecule_metadata"]["natoms"] == len(atoms)
        assert output["parameters"]["epsilon"] == 1.0
        assert output["parameters"]["sigma"] == 1.0
        assert output["parameters"]["rc"] == 3
        assert output["parameters"]["ro"] == 0.66 * 3
        assert output["results"]["energy"] == pytest.approx(1.772068860679255)

        output = static_job(atoms, epsilon=2.0, rc=0.5)
        assert output["molecule_metadata"]["natoms"] == len(atoms)
        assert output["parameters"]["epsilon"] == 2.0
        assert output["parameters"]["sigma"] == 1.0
        assert output["parameters"]["rc"] == 0.5
        assert output["parameters"]["ro"] == 0.66 * 0.5
        assert output["results"]["energy"] == pytest.approx(0.0)


def test_static_job_v2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Pt")
    atoms[0].symbol = "Au"
    assert static_job(atoms)


def test_static_job_v3(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Pt")
    atoms.pbc = False
    atoms[0].symbol = "Au"
    assert static_job(atoms)


def test_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O")

    output = relax_job(atoms)
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert output["parameters"]["epsilon"] == 1.0
    assert output["parameters"]["sigma"] == 1.0
    assert output["parameters"]["rc"] == 3
    assert output["parameters"]["ro"] == 0.66 * 3
    assert output["results"]["energy"] == pytest.approx(-2.983561029599189)
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01

    output = relax_job(atoms, fmax=0.03, epsilon=2.0, rc=0.5)
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert output["parameters"]["epsilon"] == 2.0
    assert output["parameters"]["sigma"] == 1.0
    assert output["parameters"]["rc"] == 0.5
    assert output["parameters"]["ro"] == 0.66 * 0.5
    assert output["results"]["energy"] == pytest.approx(0.0)
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.03


def test_freq_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O")

    output = freq_job(relax_job(atoms)["atoms"])
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert output["parameters"]["epsilon"] == 1.0
    assert output["parameters"]["sigma"] == 1.0
    assert output["parameters"]["rc"] == 3
    assert output["parameters"]["ro"] == 0.66 * 3
    assert len(output["results"]["vib_freqs_raw"]) == 3 * len(atoms)
    assert len(output["results"]["vib_freqs"]) == 3 * len(atoms) - 6
    assert len(output["parameters_thermo"]["vib_freqs"]) == 3 * len(atoms) - 6
    assert output["parameters_thermo"]["n_imag"] == 0
