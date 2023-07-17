import multiprocessing
import os
from pathlib import Path
from shutil import copy

import pytest
from ase.build import molecule
from numpy.testing import assert_allclose

from quacc.recipes.orca.core import relax_job, static_job

FILE_DIR = Path(__file__).resolve().parent
ORCA_DIR = os.path.join(FILE_DIR, "orca_run")
BAD_ORCA_DIR = os.path.join(FILE_DIR, "orca_failed_run")


def prep_files():
    for f in os.listdir(ORCA_DIR):
        copy(os.path.join(ORCA_DIR, f), f)


def prep_files_bad():
    for f in os.listdir(BAD_ORCA_DIR):
        copy(os.path.join(BAD_ORCA_DIR, f), f)


def test_static_job(monkeypatch, tmpdir):
    tmpdir.chdir()
    prep_files()

    atoms = molecule("H2")

    output = static_job(atoms)
    assert output["natoms"] == len(atoms)
    assert (
        output["parameters"]["orcasimpleinput"]
        == "wb97x-d3bj def2-tzvp sp slowconv normalprint xyzfile"
    )
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["spin_multiplicity"] == 1
    assert output["charge"] == 0

    output = static_job(
        atoms,
        charge=-2,
        multiplicity=3,
        input_swaps={"def2-svp": True, "def2-tzvp": None},
        block_swaps={"%scf maxiter 300 end": True},
    )
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == -2
    assert output["parameters"]["mult"] == 3
    assert (
        output["parameters"]["orcasimpleinput"]
        == "wb97x-d3bj sp slowconv normalprint xyzfile def2-svp"
    )
    assert "%scf maxiter 300 end" in output["parameters"]["orcablocks"]

    atoms = molecule("H2")
    monkeypatch.setenv("mpirun", "test")
    output = static_job(atoms)
    nprocs = multiprocessing.cpu_count()
    assert f"%pal nprocs {nprocs} end" in output["parameters"]["orcablocks"]


def test_relax_job(monkeypatch, tmpdir):
    tmpdir.chdir()
    prep_files()

    atoms = molecule("H2")

    output = relax_job(atoms)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert (
        output["parameters"]["orcasimpleinput"]
        == "wb97x-d3bj def2-tzvp opt slowconv normalprint xyzfile"
    )

    output = relax_job(
        atoms,
        charge=-2,
        multiplicity=3,
        input_swaps={
            "hf": True,
            "wb97x-d3bj": None,
            "def2-svp": True,
            "def2-tzvp": None,
        },
        block_swaps={"%scf maxiter 300 end": True},
    )
    assert output["natoms"] == len(atoms)
    assert (
        output["parameters"]["orcasimpleinput"]
        == "opt slowconv normalprint xyzfile hf def2-svp"
    )
    assert "%scf maxiter 300 end" in output["parameters"]["orcablocks"]
    assert "trajectory" in output["attributes"]
    assert len(output["attributes"]["trajectory"]) > 1
    assert (
        output["attributes"]["trajectory"][0] != output["attributes"]["trajectory"][-1]
    )
    assert_allclose(
        output["attributes"]["trajectory"][-1]["atoms"].get_positions(),
        output["atoms"].get_positions(),
        rtol=1e-5,
    )

    atoms = molecule("H2")
    monkeypatch.setenv("mpirun", "test")
    output = relax_job(atoms)
    nprocs = multiprocessing.cpu_count()
    assert f"%pal nprocs {nprocs} end" in output["parameters"]["orcablocks"]


def test_bad_relax_job(tmpdir):
    tmpdir.chdir()
    prep_files_bad()

    atoms = molecule("H2")
    with pytest.raises(ValueError):
        relax_job(atoms)
