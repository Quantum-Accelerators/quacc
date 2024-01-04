import os

import pytest
from ase.build import molecule

from quacc import Remove
from quacc.recipes.orca.core import ase_relax_job, relax_job, static_job


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = static_job(atoms, charge=0, spin_multiplicity=1, nprocs=1)
    assert output["natoms"] == len(atoms)
    assert (
        output["parameters"]["orcasimpleinput"]
        == "def2-tzvp normalprint slowconv sp wb97x-d3bj xyzfile"
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
        == "def2-svp normalprint slowconv sp wb97x-d3bj xyzfile"
    )
    assert "%scf maxiter 300 end" in output["parameters"]["orcablocks"]
    assert output.get("attributes")


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
        == "wb97x-d3bj def2-tzvp opt slowconv normalprint xyzfile"
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
        output["parameters"]["orcasimpleinput"]
        == "def2-svp hf normalprint opt slowconv xyzfile"
    )
    assert "%scf maxiter 300 end" in output["parameters"]["orcablocks"]
    assert output.get("trajectory")
    assert len(output["trajectory"]) > 1
    assert output.get("attributes")


def test_ase_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = ase_relax_job(atoms, opt_params={"fmax": 0.1}, nprocs=1)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert (
        output["parameters"]["orcasimpleinput"]
        == "def2-tzvp engrad normalprint slowconv wb97x-d3bj xyzfile"
    )
    assert output["fmax"] == 0.1
    assert output.get("trajectory")
    assert output.get("trajectory_results")
    assert output.get("attributes")
