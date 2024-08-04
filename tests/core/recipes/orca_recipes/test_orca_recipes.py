from __future__ import annotations

import os
from pathlib import Path

import pytest
from ase.build import molecule

from quacc.recipes.orca.core import (
    ase_quasi_irc_job,
    ase_relax_job,
    freq_job,
    relax_job,
    static_job,
)

FILE_DIR = Path(__file__).parent


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


@pytest.mark.skipif(os.name == "nt", reason="mpirun not available on Windows")
def test_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = relax_job(atoms, charge=0, spin_multiplicity=1, nprocs=2)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert (
        output["parameters"]["orcasimpleinput"]
        == "def2-tzvp normalprint opt wb97x-d3bj xyzfile"
    )

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


@pytest.mark.skipif(os.name == "nt", reason="mpirun not available on Windows")
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
        nprocs=2,
    )
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["parameters"]["orcasimpleinput"] == "def2-svp freq hf opt xyzfile"


@pytest.mark.skipif(os.name == "nt", reason="mpirun not available on Windows")
def test_ase_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = ase_relax_job(atoms, opt_params={"fmax": 0.1}, nprocs=2)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert (
        output["parameters"]["orcasimpleinput"]
        == "def2-tzvp engrad normalprint wb97x-d3bj xyzfile"
    )
    assert output["parameters_opt"]["fmax"] == 0.1
    assert output.get("trajectory")
    assert output.get("trajectory_results")


@pytest.mark.skipif(os.name == "nt", reason="mpirun not available on Windows")
def test_ase_relax_job_store(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O")
    output = ase_relax_job(
        atoms, opt_params={"store_intermediate_results": True}, nprocs=2
    )
    nsteps = len(output["trajectory"])
    for i in range(nsteps):
        assert f"step{i}" in os.listdir(output["dir_name"])
        assert "orca.xyz.gz" in os.listdir(Path(output["dir_name"], f"step{i}"))


@pytest.mark.skipif(os.name == "nt", reason="mpirun not available on Windows")
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
        nprocs=2,
    )
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["parameters"]["orcasimpleinput"] == "def2-svp freq hf xyzfile"


@pytest.mark.skipif(os.name == "nt", reason="mpirun not available on Windows")
def test_ase_quasi_irc_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    mode = [[0.0, 0.0, 0.1], [0.0, 0.1, 0.0]]

    output = ase_quasi_irc_job(
        atoms,
        mode,
        charge=0,
        spin_multiplicity=1,
        perturb_magnitude=0.01,
        direction="reverse",
        xc="hf",
        basis="def2-svp",
        orcasimpleinput=["#normalprint"],
        nprocs=2,
    )
    assert output["natoms"] == len(atoms)
    assert output["atoms"] != atoms
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["parameters"]["orcasimpleinput"] == "def2-svp engrad hf xyzfile"
