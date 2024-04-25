from __future__ import annotations

import os
from pathlib import Path

import pytest
from ase.build import molecule
from ase.io import read

from quacc.recipes.orca.core import (
    ase_quasi_irc_perturb_job,
    ase_relax_job,
    freq_job,
    relax_job,
    static_job,
)

FILE_DIR = Path(__file__).parent


@pytest.fixture()
def test_atoms():
    return read(FILE_DIR / "xyz" / "ts_test.xyz")


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = static_job(atoms, charge=0, spin_multiplicity=1, nprocs=1)
    assert output["natoms"] == len(atoms)
    assert (
        output["parameters"]["orcasimpleinput"]
        == "def2-tzvp engrad normalprint wb97x-d3bj xyzfile"
    )
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["spin_multiplicity"] == 1
    assert output["charge"] == 0
    assert output.get("attributes")


@pytest.mark.skipif(os.name == "nt", reason="mpirun not available on Windows")
def test_static_job_parallel(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = static_job(
        atoms,
        charge=-2,
        spin_multiplicity=3,
        orcasimpleinput=["def2-svp", "#def2-tzvp"],
        orcablocks=["%scf maxiter 300 end"],
        nprocs=2,
    )
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == -2
    assert output["parameters"]["mult"] == 3
    assert (
        output["parameters"]["orcasimpleinput"]
        == "def2-svp engrad normalprint wb97x-d3bj xyzfile"
    )
    assert "%scf maxiter 300 end" in output["parameters"]["orcablocks"]
    assert output.get("attributes")


def test_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = relax_job(atoms, charge=0, spin_multiplicity=1)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert (
        output["parameters"]["orcasimpleinput"]
        == "def2-tzvp normalprint opt wb97x-d3bj xyzfile"
    )
    assert output["trajectory"][0] != output["trajectory"][-1]

    output = relax_job(
        atoms,
        charge=-2,
        spin_multiplicity=3,
        orcasimpleinput=["def2-svp", "#def2-tzvp", "#wb97x-d3bj", "hf"],
        orcablocks=["%scf maxiter 300 end"],
        nprocs=2,
    )
    assert output["natoms"] == len(atoms)
    assert (
        output["parameters"]["orcasimpleinput"] == "def2-svp hf normalprint opt xyzfile"
    )
    assert (
        output["parameters"]["orcablocks"] == "%pal nprocs 2 end\n%scf maxiter 300 end"
    )
    assert output.get("trajectory")
    assert len(output["trajectory"]) > 1
    assert output.get("attributes")


def test_relax_freq_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = relax_job(
        atoms,
        xc="hf",
        basis="def2-svp",
        charge=0,
        spin_multiplicity=1,
        orcasimpleinput=["#normalprint"],
        run_freq=True,
    )
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["parameters"]["orcasimpleinput"] == "def2-svp freq hf opt xyzfile"
    assert output["trajectory"][0] != output["trajectory"][-1]


def test_ase_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = ase_relax_job(atoms, opt_params={"fmax": 0.1}, nprocs=1)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert (
        output["parameters"]["orcasimpleinput"]
        == "def2-tzvp engrad normalprint wb97x-d3bj xyzfile"
    )
    assert output["fmax"] == 0.1
    assert output.get("trajectory")
    assert output.get("trajectory_results")
    assert output.get("attributes")


def test_ase_relax_job_store(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O")
    output = ase_relax_job(atoms, opt_params={"store_intermediate_results": True})
    nsteps = len(output["trajectory"])
    for i in range(nsteps):
        assert f"step{i}" in os.listdir(output["dir_name"])
        assert "orca.xyz.gz" in os.listdir(Path(output["dir_name"], f"step{i}"))
    assert len(output["steps"]) == nsteps
    assert "attributes" in output["steps"][0]


def test_freq_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = freq_job(
        atoms,
        xc="hf",
        basis="def2-svp",
        charge=0,
        spin_multiplicity=1,
        orcasimpleinput=["#normalprint"],
    )
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["parameters"]["orcasimpleinput"] == "def2-svp freq hf xyzfile"
    assert output.get("attributes")


def test_ase_quasi_irc_perturb_job(test_atoms, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    mode = [
        [-0.164243, 0.289973, 0.026016],
        [0.111755, -0.021282, -0.002864],
        [0.011866, -0.071662, -0.040428],
        [-0.087153, 0.038657, -0.038008],
        [-0.017777, 0.013758, 0.001151],
        [0.030409, -0.18767, 0.028847],
        [0.750577, -0.378458, 0.179569],
        [0.042058, 0.035122, 0.024727],
        [-0.008014, -0.000956, -0.009153],
        [-0.052972, -0.176662, -0.077928],
        [0.037381, 0.036482, 0.028861],
        [0.043279, 0.040925, 0.022537],
        [0.035434, 0.032613, 0.019516],
        [-0.002674, -0.03398, 0.011123],
        [-0.006118, -0.009193, -0.122432],
        [0.014124, -0.035613, 0.097518],
    ]

    output = ase_quasi_irc_perturb_job(
        test_atoms,
        mode,
        charge=0,
        spin_multiplicity=1,
        perturb_magnitude=0.01,
        direction="reverse",
        xc="hf",
        basis="def2-svp",
        orcasimpleinput=["#normalprint"],
    )
    assert output["natoms"] == len(test_atoms)
    assert output["atoms"] != test_atoms
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["parameters"]["orcasimpleinput"] == "def2-svp engrad hf xyzfile"
    assert output.get("attributes")
