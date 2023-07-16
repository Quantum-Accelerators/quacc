from shutil import which

import pytest
from ase.build import molecule
from numpy.testing import assert_allclose

from quacc import SETTINGS
from quacc.recipes.orca.core import relax_job, static_job

has_orca = bool(which(SETTINGS.ORCA_CMD))

DEFAULT_SETTINGS = SETTINGS.copy()


def setup_function():
    SETTINGS.CREATE_UNIQUE_WORKDIR = False


def teardown_function():
    SETTINGS.CREATE_UNIQUE_WORKDIR = DEFAULT_SETTINGS.CREATE_UNIQUE_WORKDIR


@pytest.mark.skipif(has_orca is False, reason="ORCA not installed")
def test_static_job(tmpdir):
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


@pytest.mark.skipif(has_orca is False, reason="ORCA not installed")
def test_relax_job(tmpdir):
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


def test_unique_workdir(tmpdir):
    SETTINGS.CREATE_UNIQUE_WORKDIR = True
    test_static_job(tmpdir)
    test_relax_job(tmpdir)

