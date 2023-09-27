from shutil import which

import pytest
from ase.build import molecule
from numpy.testing import assert_allclose

from quacc import SETTINGS

has_orca = bool(which("orca"))

pytestmark = pytest.mark.skipif(
    not has_orca or SETTINGS.WORKFLOW_ENGINE != "local",
    reason="Need ORCA and Need to use local as workflow manager to run this test.",
)


def test_static_job(tmpdir):
    from quacc.recipes.orca.core import static_job

    tmpdir.chdir()

    atoms = molecule("H2")

    output = static_job(atoms, 0, 1)
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
        -2,
        3,
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


def test_relax_job(tmpdir):
    from quacc.recipes.orca.core import relax_job

    tmpdir.chdir()

    atoms = molecule("H2")

    output = relax_job(atoms, 0, 1)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert (
        output["parameters"]["orcasimpleinput"]
        == "wb97x-d3bj def2-tzvp opt slowconv normalprint xyzfile"
    )

    output = relax_job(
        atoms,
        -2,
        3,
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
