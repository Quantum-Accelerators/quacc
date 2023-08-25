import multiprocessing

import pytest
from ase.build import molecule

from quacc import SETTINGS
from quacc.recipes.orca.core import relax_job, static_job


@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_static_job(monkeypatch, tmpdir):
    tmpdir.chdir()

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


@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_relax_job(monkeypatch, tmpdir):
    tmpdir.chdir()

    atoms = molecule("H2")

    output = relax_job(atoms)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert (
        output["parameters"]["orcasimpleinput"]
        == "wb97x-d3bj def2-tzvp opt slowconv normalprint xyzfile"
    )
    assert (
        output["attributes"]["trajectory"][0] != output["attributes"]["trajectory"][-1]
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

    atoms = molecule("H2")
    monkeypatch.setenv("mpirun", "test")
    output = relax_job(atoms)
    nprocs = multiprocessing.cpu_count()
    assert f"%pal nprocs {nprocs} end" in output["parameters"]["orcablocks"]
